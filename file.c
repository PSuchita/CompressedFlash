/*
 *  linux/fs/ext2/file.c
 *
 * Copyright (C) 1992, 1993, 1994, 1995
 * Remy Card (card@masi.ibp.fr)
 * Laboratoire MASI - Institut Blaise Pascal
 * Universite Pierre et Marie Curie (Paris VI)
 *
 *  from
 *
 *  linux/fs/minix/file.c
 *
 *  Copyright (C) 1991, 1992  Linus Torvalds
 *
 *  ext2 fs regular file handling primitives
 *
 *  64-bit file support on 64-bit platforms by Jakub Jelinek
 * 	(jj@sunsite.ms.mff.cuni.cz)
 */

#ifdef CONFIG_EXT2_COMPRESS
#include <linux/fs.h>
#include <linux/ext2_fs_c.h>
#include <linux/buffer_head.h>
#include <asm/uaccess.h>
#include <linux/kmod.h>
#include <linux/slab.h>
#include <linux/swap.h>
#include <linux/pagemap.h>
#include <linux/quotaops.h>
#include <linux/writeback.h>
#else
#include <linux/time.h>
#include <linux/pagemap.h>
#include <linux/quotaops.h>
#include "ext2.h"
#endif


#include "xattr.h"
#include "acl.h"

/*
 * Called when filp is released. This happens when all file descriptors
 * for a single struct file are closed. Note that different open() calls
 * for the same file yield different struct file structures.
 */

/*
 * pjm 1998-01-09: I would note that this is different from `when no
 * process has the inode open'.
 */
static int ext2_release_file (struct inode * inode, struct file * filp)
{
#ifdef CONFIG_EXT2_COMPRESS
	/*
	 * Now's as good a time as any to clean up wrt compression.
	 * Previously (before 2.1.4x) we waited until
	 * ext2_put_inode(), but now the dcache sometimes delays that
	 * call until umount time.
	 */
	//printk(KERN_DEBUG "ext2_release_file: pid=%d, i_ino=%lu, i_count=%d\n", current->pid, inode->i_ino,  atomic_read(&inode->i_count));
	
	if (S_ISREG (inode->i_mode)
	    && inode->i_nlink
	    && (EXT2_I(inode)->i_compr_flags & EXT2_CLEANUP_FL)) {
#ifdef EXT2_COMPR_REPORT_PUT
		printk(KERN_DEBUG "ext2_release_file: pid=%d, i_ino=%lu, i_count=%d\n", current->pid, inode->i_ino,  atomic_read(&inode->i_count));
#endif
		/*
		 * todo: See how the return code of
		 * ext2_release_file() is used, and decide whether it
		 * might be appropriate to pass any errors to
		 * caller.
		 */
		//dump_stack();
		(void) ext2_cleanup_compressed_inode (inode);
	}
	
#endif
	if (filp->f_mode & FMODE_WRITE) {
		mutex_lock(&EXT2_I(inode)->truncate_mutex);
		ext2_discard_reservation(inode);
		mutex_unlock(&EXT2_I(inode)->truncate_mutex);
	}
	return 0;
}

int ext2_fsync(struct file *file, loff_t start, loff_t end, int datasync)
{
	int ret;
	struct super_block *sb = file->f_mapping->host->i_sb;
	struct address_space *mapping = sb->s_bdev->bd_inode->i_mapping;

	ret = generic_file_fsync(file, start, end, datasync);
	if (ret == -EIO || test_and_clear_bit(AS_EIO, &mapping->flags)) {
		/* We don't really know where the IO error happened... */
		ext2_error(sb, __func__,
			   "detected IO error when writing metadata buffers");
		ret = -EIO;
	}
	return ret;
}

#ifdef CONFIG_EXT2_COMPRESS
struct page_cluster {
	struct page *	page;
	loff_t		pos;
	unsigned	bytes;
	unsigned long	offset;
	unsigned char	in_range;
	const char *	buf;
};

#define PAGE_IN_RANGE 1
#define PAGE_KMAPPED  2


/**
 * generic_osync_inode - flush all dirty data for a given inode to disk
 * @inode: inode to write
 * @mapping: the address_space that should be flushed
 * @what:  what to write and wait upon
 *
 * This can be called by file_write functions for files which have the
 * O_SYNC flag set, to flush dirty writes to disk.
 *
 * @what is a bitmask, specifying which part of the inode's data should be
 * written and waited upon.
 *
 *    OSYNC_DATA:     i_mapping's dirty data
 *    OSYNC_METADATA: the buffers at i_mapping->private_list
 *    OSYNC_INODE:    the inode itself
 */

/* mw: see generic_osync_inode() in kernel<2.6.30 for orginal method.
       basically we want all of it:  OSYNC_DATA and OSYNC_METADATA  and OSYNC_INODE */
int ex_generic_osync_inode(struct inode *inode, struct address_space *mapping) //, int what)
{
        int err = 0;
        int need_write_inode_now = 0;
        int err2;

        err = filemap_fdatawrite(mapping);

        err2 = sync_mapping_buffers(mapping);
        if (!err)
 		err = err2;
 
        err2 = filemap_fdatawait(mapping);
        if (!err)
 		err = err2;
 
        /* check if data is dirty */
        spin_lock(&inode->i_lock);
        if (inode->i_state & I_DIRTY)
        	need_write_inode_now = 1;
        spin_unlock(&inode->i_lock);

        if (need_write_inode_now) {
                err2 = write_inode_now(inode, 1);
                if (!err)
                        err = err2;
        }
        else
                inode_sync_wait(inode);

        return err;
}


/*
 * Write to a file through the page cache.
 *
 * We currently put everything into the page cache prior to writing it.
 * This is not a problem when writing full pages. With partial pages,
 * however, we first have to read the data into the cache, then
 * dirty the page, and finally schedule it for writing. Alternatively, we
 * could write-through just the portion of data that would go into that
 * page, but that would kill performance for applications that write data
 * line by line, and it's prone to race conditions.
 *
 * Note that this routine doesn't try to keep track of dirty pages. Each
 * file system has to do this all by itself, unfortunately.
 *                                                    okir@monad.swb.de
 */
ssize_t
ext2_file_write(struct file *file,const char *buf,size_t count,loff_t *ppos)
{
	struct address_space *mapping = file->f_dentry->d_inode->i_mapping;
	struct inode	*inode = mapping->host;
	unsigned long	limit = current->signal->rlim[RLIMIT_FSIZE].rlim_cur, written, last_index;	   /* last page index */
	loff_t	pos;
	long	status;
	int		err;
	unsigned	bytes;
	u32		comprblk_mask=0;
	struct ext2_inode_info *ei = EXT2_I(inode);

	if (!(ei->i_flags & (EXT2_COMPR_FL|EXT2_COMPRBLK_FL)) 
#undef DUD //mw: I think this is a buggy bug-fix
#ifdef DUD
			|| (count < inode->i_sb->s_blocksize) 
#endif
	)		
	{
		return do_sync_write(file, buf, count, ppos);
	}

	if ((ssize_t) count < 0)
		return -EINVAL;

	if (!access_ok(VERIFY_READ, buf, count))
		return -EFAULT;

#ifdef EXT2_COMPR_REPORT_MUTEX
    printk(KERN_DEBUG "EXT2_FILE_WRITE_LOCK of PID %u @ inode:%lu\n", current->pid, inode->i_ino );
#endif
	mutex_lock(&inode->i_mutex);
	/* mw:	down_read(&inode->i_alloc_sem); // as used by ocsf2 TLL 02/21/07 
		was removed with kernel 3.1 */
	atomic_inc(&inode->i_dio_count);

	pos = *ppos;
	err = -EINVAL;
	if (pos < 0)
		goto out;

	written = 0;

	/* FIXME: this is for backwards compatibility with 2.4 */
	if (!S_ISBLK(inode->i_mode) && file->f_flags & O_APPEND)
	{
		pos = inode->i_size;
	}

	/*
	 * Check whether we've reached the file size limit.
	 */
	err = -EFBIG;

	if (limit != RLIM_INFINITY) {
		if (pos >= limit) {
			send_sig(SIGXFSZ, current, 0);
			goto out;
		}
		if (pos > 0xFFFFFFFFULL || count > limit - (u32)pos) {
			/* send_sig(SIGXFSZ, current, 0); */
			count = limit - (u32)pos;
		}
	}

	/*
	 *      LFS rule
	 */
	if ( pos + count > MAX_NON_LFS && !(file->f_flags&O_LARGEFILE)) {
		if (pos >= MAX_NON_LFS) {
			send_sig(SIGXFSZ, current, 0);
			goto out;
		}
		if (count > MAX_NON_LFS - (u32)pos) {
			/* send_sig(SIGXFSZ, current, 0); */
			count = MAX_NON_LFS - (u32)pos;
		}
	}

	/*
	 *	Are we about to exceed the fs block limit ?
	 *
	 *	If we have written data it becomes a short write
	 *	If we have exceeded without writing data we send
	 *	a signal and give them an EFBIG.
	 *
	 *	Linus frestrict idea will clean these up nicely..
	 */
	if (!S_ISBLK(inode->i_mode)) {
		if (pos >= inode->i_sb->s_maxbytes) {
			if (count || pos > inode->i_sb->s_maxbytes) {
				send_sig(SIGXFSZ, current, 0);
				err = -EFBIG;
				goto out;
			}
			/* zero-length writes at ->s_maxbytes are OK */
		}

		if (pos + count > inode->i_sb->s_maxbytes)
			count = inode->i_sb->s_maxbytes - pos;
	} else { 
		if (bdev_read_only(inode->i_sb->s_bdev)) {
			err = -EPERM;
			goto out;
		}
		if (pos >= inode->i_size) {
			if (count || pos > inode->i_size) {
				err = -ENOSPC;
				goto out;
			}
		}

		if (pos + count > inode->i_size)
		{
			count = inode->i_size - pos;			
		}
	}

	err = 0;
	if (count == 0)
		goto out;

	status  = 0;

	if (file->f_flags & O_DIRECT)
	{
		err = -EINVAL;
		goto out;
	}
	/*
	 *	We must still check for EXT2_ECOMPR_FL, as it may have been
	 *	set after we got the write permission to this file.
	 */
	if ((ei->i_flags & (EXT2_ECOMPR_FL | EXT2_NOCOMPR_FL))   == (EXT2_ECOMPR_FL | 0)) 
	{
		err = -EXT2_ECOMPR;
		goto out;
	}

	should_remove_suid(file->f_dentry);
	inode->i_ctime = inode->i_mtime = CURRENT_TIME;
	mark_inode_dirty_sync(inode);

	if ((pos+count) > inode->i_size)
		last_index = (pos+count-1) >> PAGE_CACHE_SHIFT;
	else
		last_index = (inode->i_size-1) >> PAGE_CACHE_SHIFT;

	comprblk_mask = ei->i_flags | ~EXT2_COMPRBLK_FL;

	//mw: now do it cluster-wise
	do {
		//unsigned long index, offset, clusters_page_index0, 
		unsigned long index, nextClusterFirstByte, cluster_compressed=0;
		u32  cluster=0;
		status = -ENOMEM;	/* we'll assign it later anyway */

#ifdef EXT2_COMPRESS_WHEN_CLU	
		ei->i_flags |= EXT2_COMPRBLK_FL;
		assert( (file->f_flags & O_DIRECT) == 0);
		assert(mapping_mapped(inode->i_mapping) == 0);
#endif			

		index = pos >> PAGE_CACHE_SHIFT;	/*mw: pageindex (start)*/
		cluster = ext2_page_to_cluster(inode, index);		

		/*
		 * We decompress the cluster if needed, and write
		 * the data as normal.  The cluster will be
		 * compressed again when the inode is cleaned up.
		 */
		if ((comprblk_mask == ~(u32)0)
		    && !(ei->i_flags & EXT2_NOCOMPR_FL)) {
		      /* AUFFÄLLIG 2*/
			/* assert (block == pos >> inode->i_sb->s_blocksize_bits); */

			cluster_compressed = ext2_cluster_is_compressed_fn(inode, cluster);
			if (cluster_compressed < 0) {
				if (! written)
					written = cluster_compressed;
				break;
			}
		}

		if (cluster_compressed > 0) {
			/* Here, decompression take place  */
			cluster_compressed = ext2_decompress_cluster(inode, cluster);
			if (cluster_compressed < 0) {
				if (! written) {
					written = cluster_compressed;
				}
				break;
			}
		}
			
		nextClusterFirstByte = (ext2_cluster_page0(inode, cluster+1) * PAGE_CACHE_SIZE);	
		bytes = nextClusterFirstByte - pos;   /*mw: bytes todo in this cluster*/
		if (bytes > count) {
			bytes = count; /*mw: if end of data*/
		}
	
#ifdef EXT2_COMPR_DEBUG
		//assert we stay inside the cluster!
		{
			int endpos;
			int endindex;
			int endcluster;
			unsigned long thisClusterFirstByte;
			int relstart, relend, startblock, endblock;
			
			thisClusterFirstByte = (ext2_cluster_page0(inode, cluster) * PAGE_CACHE_SIZE);
			
			relstart =   pos - thisClusterFirstByte;
			relend   =  bytes + relstart;
			
			startblock = relstart >> 10; 
			endblock = relend >> 10;
			
					
			endpos = pos + bytes;
			//printk("do_sync_write cluster %d: inode:%lu, \t start:%i(%i), end:%i(%i), \t ccount:%d \t tcount:%d\n", cluster , inode->i_ino, relstart, startblock, relend , endblock,  (int)bytes, count);
	 		endindex = (endpos-1) >> PAGE_CACHE_SHIFT;	/*mw: pageindex (start)*/
			endcluster = ext2_page_to_cluster(inode, endindex);
			assert(cluster == endcluster); 
		}
#endif
	
		//mw: must unlock here, do_sync_write() will aquire the mutex again
		mutex_unlock(&inode->i_mutex);	
		
		//mw: this is pretty clever: we use the generic method now :-)
		//printk("do_sync_write cluster %d, mapped:%i\n", cluster, mapping_mapped(inode->i_mapping));
		//status = do_sync_write_nolock(file, buf,  bytes, &pos); //without locking mutex
		status = do_sync_write(file, buf,  bytes, &pos);  //with locking mutex
		assert(status>=0);
		
		mutex_lock(&inode->i_mutex);
				
		written += status;
		count -= status;
		buf += status;

#ifdef EXT2_COMPRESS_WHEN_CLU
		assert (ei->i_flags & EXT2_COMPRBLK_FL);
		if ((ei->i_flags & EXT2_COMPR_FL)
		    && (ext2_offset_is_clu_boundary(inode, pos)) ) {
			
			if  (mapping_mapped(inode->i_mapping) == 0 ) 
			/*
			 * Pierre Peiffer: For file mapped (via mmap, I mean),
			 * compression will occure when releasing the file.
			 * We must, in this case, avoid the pages (possibly
			 * mapped by a process) to be compressed under them.
			 */
			{
				int error;
				assert(mapping_mapped(inode->i_mapping) == 0);
				error = ext2_compress_cluster(inode, cluster);
				/*if (ext2_cluster_is_compressed_fn(inode, cluster))
					ext2_decompress_cluster(inode, cluster);*/
				assert(mapping_mapped(inode->i_mapping) == 0);
				/*
				 * Actually, raising write_error may be a
				 * mistake.  For example,
				 * ext2_cleanup_compressed_cluster() doesn't
				 * usually return any errors to user.  todo:
				 * Have a look at ext2_compress_cluster, and
				 * check whether its errors are such that they
				 * should be returned to user.  Some of the
				 * will be, of course, but it might be
				 * possible for it to return without
				 * change.
				 */
				if (error > 0)
					comprblk_mask = ~(u32)0;
			} else {
#ifdef EXT2_COMPR_REPORT
				char bdn[BDEVNAME_SIZE];
				bdevname(inode->i_sb->s_bdev, bdn);
#endif

				trace_e2c("ext2_file_write: (dev. %s): "
				    "ino=%ld, cluster=%d: file mapped, does "
				    "not compress cluster\n",
				    bdn, inode->i_ino, cluster);
				ei->i_flags |= EXT2_DIRTY_FL;
				ei->i_compr_flags |= EXT2_CLEANUP_FL;
			}
		}
#endif

	} while (count);
	*ppos = pos;

	/*
	 * For now, when the user asks for O_SYNC, we'll actually
	 * provide O_DSYNC.
	 */
	if (status >= 0) {
		if ((file->f_flags & O_SYNC) || IS_SYNC(inode)) {
			/*if (ei->i_compr_flags & EXT2_OSYNC_INODE) {
				osync_already = 1;
			} else {
				osync_already = 0;
				ei->i_compr_flags |= EXT2_OSYNC_INODE;
			}*/
			/* Should 2nd arg be inode->i_mapping? */
			status = ex_generic_osync_inode(inode, file->f_mapping
				/*, OSYNC_METADATA|OSYNC_DATA*/);
			/*if (osync_already == 0) {
				ei->i_compr_flags &= ~EXT2_OSYNC_INODE;
			}*/
		}
	}

	err = written ? written : status;

# ifdef EXT2_COMPRESS_WHEN_CLU
	//mw: ext2_compress_cluster() might remove EXT2_COMPRBLK_FL
	//if the file does not compress at all. this is NO error: remove next line?
	//assert (ei->i_flags & EXT2_COMPRBLK_FL);
	
	ei->i_flags &= comprblk_mask;
	if ( (ei->i_flags & EXT2_COMPR_FL)
	    && (!ext2_offset_is_clu_boundary(inode, pos)) )
	{
		ei->i_flags |= EXT2_DIRTY_FL;
		ei->i_compr_flags |= EXT2_CLEANUP_FL;
	}

# else
	if (ei->i_flags & EXT2_COMPR_FL) {
		ei->i_flags |= EXT2_DIRTY_FL;
		ei->i_compr_flags |= EXT2_CLEANUP_FL;
	}
# endif
out:

#ifdef EXT2_COMPR_REPORT_MUTEX
    printk(KERN_DEBUG "EXT2_FILE_WRITE_UNLOCK of PID %u @ inode:%lu\n", current->pid, inode->i_ino);
#endif
	/* mw:	up_read(&inode->i_alloc_sem); // as used by ocsf2 TLL 02/21/07 
		was removed with kernel 3.1 */
	inode_dio_done(inode);
	mutex_unlock(&inode->i_mutex); 
	return err;	
}

/*
 * Called when an inode is about to be open.
 * We use this to disallow opening RW large files on 32bit systems if
 * the caller didn't specify O_LARGEFILE.  On 64bit systems we force
 * on this flag in sys_open.
 * Prevent opening compressed file with O_DIRECT.
 */
static int ext2_file_open(struct inode * inode, struct file * filp)
{
	if ((filp->f_flags & O_DIRECT) && (EXT2_I(inode)->i_flags &
	    (EXT2_COMPR_FL|EXT2_COMPRBLK_FL)))
		return -EINVAL;
	if (!(filp->f_flags & O_LARGEFILE) && inode->i_size > MAX_NON_LFS)
		return -EFBIG;
			
  	return 0;
 }
#endif /* CONFIG_EXT2_COMPRESS*/

/*
 * We have mostly NULL's here: the current defaults are ok for
 * the ext2 filesystem.
 */
const struct file_operations ext2_file_operations = {
	.llseek		= generic_file_llseek,
	.read		= do_sync_read,
#ifdef CONFIG_EXT2_COMPRESS
        .write          = ext2_file_write,
#else
	.write		= do_sync_write,
#endif

	.aio_read	= generic_file_aio_read,
	.aio_write	= generic_file_aio_write,
	.unlocked_ioctl = ext2_ioctl,
#ifdef CONFIG_COMPAT
	.compat_ioctl	= ext2_compat_ioctl,
#endif
	.mmap		= generic_file_mmap,
#ifdef CONFIG_EXT2_COMPRESS
    	.open         	= ext2_file_open,
#else
	.open		= dquot_file_open,
#endif
	.release	= ext2_release_file,
	.fsync		= ext2_fsync,
	.splice_read	= generic_file_splice_read,
	.splice_write	= generic_file_splice_write,
};

#ifdef CONFIG_EXT2_FS_XIP
const struct file_operations ext2_xip_file_operations = {
	.llseek		= generic_file_llseek,
	.read		= xip_file_read,
	.write		= xip_file_write,
	.unlocked_ioctl = ext2_ioctl,
#ifdef CONFIG_COMPAT
	.compat_ioctl	= ext2_compat_ioctl,
#endif
	.mmap		= xip_file_mmap,
	.open		= dquot_file_open,
	.release	= ext2_release_file,
	.fsync		= ext2_fsync,
};
#endif

const struct inode_operations ext2_file_inode_operations = {
#ifdef CONFIG_EXT2_FS_XATTR
	.setxattr	= generic_setxattr,
	.getxattr	= generic_getxattr,
	.listxattr	= ext2_listxattr,
	.removexattr	= generic_removexattr,
#endif
	.setattr	= ext2_setattr,
	.get_acl	= ext2_get_acl,
	.fiemap		= ext2_fiemap,
};
