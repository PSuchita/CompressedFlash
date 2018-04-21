/*
 * linux/fs/ext2/ioctl.c
 *
 * Copyright (C) 1993, 1994, 1995
 * Remy Card (card@masi.ibp.fr)
 * Laboratoire MASI - Institut Blaise Pascal
 * Universite Pierre et Marie Curie (Paris VI)
 */

#ifdef CONFIG_EXT2_COMPRESS
#include <linux/fs.h>
#include <linux/ext2_fs_c.h>
#include <linux/kmod.h>
#include <linux/stat.h>
#else
#include "ext2.h"
#endif
#include <linux/capability.h>
#include <linux/time.h>
#include <linux/sched.h>
#include <linux/compat.h>
#include <linux/mount.h>
#include <asm/current.h>
#include <asm/uaccess.h>


#ifdef CONFIG_EXT2_COMPRESS

#ifndef MIN
# define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifdef CONFIG_GZ_HACK
static int check_name(struct inode *ino)
{
  struct dentry *dentry = list_entry(ino->i_dentry.next, struct dentry, d_alias);
  if (dentry)
    if (

           (dentry->d_name.len >= 4) &&
                 (((dentry->d_name.name[dentry->d_name.len - 2] == 'g')
                   && (dentry->d_name.name[dentry->d_name.len - 1] == 'z')
                   && ((dentry->d_name.name[dentry->d_name.len - 3] == '.')
                       || (dentry->d_name.name[dentry->d_name.len - 4] == '.')))

                  || ((dentry->d_name.name[dentry->d_name.len - 3] == 't')
                      && (dentry->d_name.name[dentry->d_name.len - 2] == 'g')
                      && (dentry->d_name.name[dentry->d_name.len - 1] == 'z')
                      && (dentry->d_name.name[dentry->d_name.len - 4] == '.')
                      && (dentry->d_name.len >= 5))

                  || ((dentry->d_name.name[dentry->d_name.len - 3] == 'p')
                      && (dentry->d_name.name[dentry->d_name.len - 2] == 'n')
                      && (dentry->d_name.name[dentry->d_name.len - 1] == 'g')
                      && (dentry->d_name.name[dentry->d_name.len - 4] == '.')
                      && (dentry->d_name.len >= 5))

                  || ((dentry->d_name.name[dentry->d_name.len - 3] == 'j')
                      && (dentry->d_name.name[dentry->d_name.len - 2] == 'p')
                      && (dentry->d_name.name[dentry->d_name.len - 1] == 'g')
                      && (dentry->d_name.name[dentry->d_name.len - 4] == '.')
                      && (dentry->d_name.len >= 5))

                  || ((dentry->d_name.name[dentry->d_name.len - 3] == 'b')
                      && (dentry->d_name.name[dentry->d_name.len - 2] == 'z')
                      && (dentry->d_name.name[dentry->d_name.len - 1] == '2')
                      && (dentry->d_name.name[dentry->d_name.len - 4] == '.')
                      && (dentry->d_name.len >= 5))

                  || ((dentry->d_name.name[dentry->d_name.len - 3] == 'm')
                      && (dentry->d_name.name[dentry->d_name.len - 2] == 'n')
                      && (dentry->d_name.name[dentry->d_name.len - 1] == 'g')
                      && (dentry->d_name.name[dentry->d_name.len - 4] == '.')
                      && (dentry->d_name.len >= 5))
                  )
       ) {
        return 1;
    }
  return 0;
}
#endif
#endif



long ext2_ioctl(struct file *filp, unsigned int cmd, unsigned long arg)
{
	struct inode *inode = filp->f_dentry->d_inode;
	struct ext2_inode_info *ei = EXT2_I(inode);
	unsigned int flags;
	unsigned short rsv_window_size;
	int ret;
#ifdef CONFIG_EXT2_COMPRESS
	unsigned long datum;
	int err;
#endif

	ext2_debug ("cmd = %u, arg = %lu\n", cmd, arg);

	switch (cmd) {
	case EXT2_IOC_GETFLAGS:
		ext2_get_inode_flags(ei);
		flags = ei->i_flags & EXT2_FL_USER_VISIBLE;
		return put_user(flags, (int __user *) arg);
	case EXT2_IOC_SETFLAGS: {
		unsigned int oldflags;

		ret = mnt_want_write_file(filp);
		if (ret)
			return ret;

		if (!inode_owner_or_capable(inode)) {
			ret = -EACCES;
			goto setflags_out;
		}

		if (get_user(flags, (int __user *) arg)) {
			ret = -EFAULT;
			goto setflags_out;
		}

		flags = ext2_mask_flags(inode->i_mode, flags);

		mutex_lock(&inode->i_mutex);
		/* Is it quota file? Do not allow user to mess with it */
		if (IS_NOQUOTA(inode)) {
			mutex_unlock(&inode->i_mutex);
			ret = -EPERM;
			goto setflags_out;
		}
		oldflags = ei->i_flags;

		/*
		 * The IMMUTABLE and APPEND_ONLY flags can only be changed by
		 * the relevant capability.
		 *
		 * This test looks nicer. Thanks to Pauline Middelink
		 */
		if ((flags ^ oldflags) & (EXT2_APPEND_FL | EXT2_IMMUTABLE_FL)) {
			if (!capable(CAP_LINUX_IMMUTABLE)) {
				mutex_unlock(&inode->i_mutex);
				ret = -EPERM;
				goto setflags_out;
			}
		}

		flags = flags & EXT2_FL_USER_MODIFIABLE;
#ifdef CONFIG_EXT2_COMPRESS
		if (S_ISREG (inode->i_mode) || S_ISDIR (inode->i_mode)) {		
			
			/* pjm 1998-01-14: In previous versions of
			     e2compr, the kernel forbade raising
			     EXT2_ECOMPR_FL from userspace.  I can't
			     think of any purpose for forbidding this,
			     and I find it useful to raise
			     EXT2_ECOMPR_FL for testing purposes, so
			     I've removed the forbidding code. */
			if (S_ISREG (inode->i_mode)
			    && (EXT2_NOCOMPR_FL
				& (flags ^ ei->i_flags))) {   // mw hint:  ^ is a (excluisive OR)
				/* NOCOMPR_FL can only be changed if
				   nobody else has the file opened.  */
				/* pjm 1998-02-16: inode->i_count is
				   useless to us because only dentries
				   use inodes now.  Unfortunately,
				   there isn't an easy way of finding
				   the equivalent.  We'd have to go
				   through all dentries using the
				   inode, and sum their d_count
				   values.  Rather than do that, I'd
				   rather get rid of the exclusion
				   constraint.  todo. */				
				//printk("i_count: %i\n", atomic_read(&inode->i_count));
				//if (atomic_read(&inode->i_count) > 1)
				//if (0)
				if (ext2_get_dcount(inode) > 1)
				{
					mutex_unlock(&inode->i_mutex); /*mw*/
					return -ETXTBSY;
				}
				else {
					/* pjm 970429: Discarding
					     cached pages is not very
					     clean, but should work. */
					/* pjm 980114: Not quite.  We
					     should also sync any
					     mappings to buffers first.
					     This isn't very important,
					     as none of the current
					     e2compr programs can
					     trigger this, but todo. */
					invalidate_remote_inode (inode);
				}
			}

			if (EXT2_COMPR_FL
			    & (flags ^ ei->i_flags)) {
				if (flags & EXT2_COMPR_FL) {
					if (ei->i_flags & EXT2_COMPRBLK_FL) {
						/* There shouldn't actually be any
						   compressed blocks, AFAIK.  However,
						   this is still possible because sometimes
						   COMPRBLK gets raised just to stop 
						   us changing cluster size at the wrong
						   time.

						   todo: Call a function that just
						   checks that there are not compressed
						   clusters, and print a warning if any are
						   found. */
					} else {
						int bits = MIN(EXT2_DEFAULT_LOG2_CLU_NBLOCKS,
							       (EXT2_LOG2_MAX_CLUSTER_BYTES
								- inode->i_sb->s_blocksize_bits));

						ei->i_log2_clu_nblocks = bits;
						ei->i_clu_nblocks = 1 << bits;
					}
					ei->i_compr_method = EXT2_DEFAULT_COMPR_METHOD;
					if (S_ISREG (inode->i_mode)) {
						//compress	
#ifdef CONFIG_GZ_HACK
						/*  mw: check for .gz-files and similar
						 *  I think this is the most clever place for
						 * rejecting files. They remain regular, uncompressed
						 * files and though can be read bypassing all 
						 * compression stuff (= fast) :-). And it seems to save 
						 * space... somehow */						
						if (check_name (inode))
						{
							//printk("non-compressable file extension\n");
							mutex_unlock(&inode->i_mutex);
							return 0;
						}
#endif
						//set flags to trigger compression later on
						flags |= EXT2_DIRTY_FL;
						ei->i_compr_flags |= EXT2_CLEANUP_FL;
					}
				} else if (S_ISREG (inode->i_mode)) {
					if (ei->i_flags & EXT2_COMPRBLK_FL) {
						int err;
						
						if (ext2_get_dcount(inode) > 1){
							mutex_unlock(&inode->i_mutex); //mw
							return -ETXTBSY;
						}
						err = ext2_decompress_inode(inode);
						if (err)
						{
				            		mutex_unlock(&inode->i_mutex); //mw
							return err;
						}
					}
					ei->i_flags &= ~EXT2_DIRTY_FL;
					ei->i_compr_flags &= ~EXT2_CLEANUP_FL;
				}
			}
		}
#endif
		flags |= oldflags & ~EXT2_FL_USER_MODIFIABLE;
#ifdef CONFIG_EXT2_COMPRESS
		/* bug fix: scrub 'B' flag from uncompressed files TLL 02/28/07 */
		if (!(flags & EXT2_COMPR_FL) && (flags & EXT2_COMPRBLK_FL) )
		{ 
				flags &= ~EXT2_COMPRBLK_FL;
		}
#endif
		ei->i_flags = flags;

		ext2_set_inode_flags(inode);
		inode->i_ctime = CURRENT_TIME_SEC;
		mutex_unlock(&inode->i_mutex);

		mark_inode_dirty(inode);
setflags_out:
		mnt_drop_write_file(filp);
		return ret;
	}
	case EXT2_IOC_GETVERSION:
		return put_user(inode->i_generation, (int __user *) arg);
	case EXT2_IOC_SETVERSION: {
		__u32 generation;

		if (!inode_owner_or_capable(inode))
			return -EPERM;
		ret = mnt_want_write_file(filp);
		if (ret)
			return ret;
		if (get_user(generation, (int __user *) arg)) {
			ret = -EFAULT;
			goto setversion_out;
		}

		mutex_lock(&inode->i_mutex);
		inode->i_ctime = CURRENT_TIME_SEC;
		inode->i_generation = generation;
		mutex_unlock(&inode->i_mutex);

		mark_inode_dirty(inode);
setversion_out:
		mnt_drop_write_file(filp);
		return ret;
	}
	case EXT2_IOC_GETRSVSZ:
		if (test_opt(inode->i_sb, RESERVATION)
			&& S_ISREG(inode->i_mode)
			&& ei->i_block_alloc_info) {
			rsv_window_size = ei->i_block_alloc_info->rsv_window_node.rsv_goal_size;
			return put_user(rsv_window_size, (int __user *)arg);
		}
		return -ENOTTY;
	case EXT2_IOC_SETRSVSZ: {

		if (!test_opt(inode->i_sb, RESERVATION) ||!S_ISREG(inode->i_mode))
			return -ENOTTY;

		if (!inode_owner_or_capable(inode))
			return -EACCES;

		if (get_user(rsv_window_size, (int __user *)arg))
			return -EFAULT;

		ret = mnt_want_write_file(filp);
		if (ret)
			return ret;

		if (rsv_window_size > EXT2_MAX_RESERVE_BLOCKS)
			rsv_window_size = EXT2_MAX_RESERVE_BLOCKS;

		/*
		 * need to allocate reservation structure for this inode
		 * before set the window size
		 */
		/*
		 * XXX What lock should protect the rsv_goal_size?
		 * Accessed in ext2_get_block only.  ext3 uses i_truncate.
		 */
		mutex_lock(&ei->truncate_mutex);
		if (!ei->i_block_alloc_info)
			ext2_init_block_alloc_info(inode);

		if (ei->i_block_alloc_info){
			struct ext2_reserve_window_node *rsv = &ei->i_block_alloc_info->rsv_window_node;
			rsv->rsv_goal_size = rsv_window_size;
		}
		mutex_unlock(&ei->truncate_mutex);
		mnt_drop_write_file(filp);
		return 0;
	}
#ifdef CONFIG_EXT2_COMPRESS
	case EXT2_IOC_GETCOMPRMETHOD:	/* Result means nothing if COMPR_FL is not set */
 		return put_user (ei->i_compr_method, (long *) arg);
	case EXT2_IOC_SETCOMPRMETHOD:
		if ((current_fsuid() != inode->i_uid) && !capable(CAP_FOWNER))
			return -EPERM;
		if (IS_RDONLY (inode))
			return -EROFS;
		if (get_user (datum, (long*) arg))
			return -EFAULT;
		if (!S_ISREG (inode->i_mode) && !S_ISDIR (inode->i_mode)) 
			return -ENOSYS;
		/* todo: Allow the below, but set initial value of
		   i_compr_meth at read_inode() time (using default if
		   !/) instead of +c time.  Same for cluster
		   size. */
		if ((unsigned) datum >= EXT2_N_METHODS)
			return -EINVAL;
		if (ei->i_compr_method != datum) {
			if ((ei->i_compr_method == EXT2_NEVER_METH)
			    && (ei->i_flags & EXT2_COMPR_FL))
				return -EPERM;
			/* If the previous method was `defer' then
			   take a look at all uncompressed clusters
			   and try to compress them.  (pjm 1997-04-16) */
			if ((ei->i_compr_method == EXT2_DEFER_METH)
			    && S_ISREG (inode->i_mode)) {
				ei->i_flags |= EXT2_DIRTY_FL;
				ei->i_compr_flags |= EXT2_CLEANUP_FL;
			}
			if ((datum == EXT2_NEVER_METH)
			    && S_ISREG (inode->i_mode)) {
				//printk("SETCOMPR\n");
				if ((ei->i_flags & EXT2_COMPRBLK_FL))
				{
					/*mw*/
					mutex_lock(&inode->i_mutex);
					if (ext2_get_dcount(inode) > 1){
						mutex_unlock(&inode->i_mutex); /*mw*/
						return -ETXTBSY;
					}
					err = ext2_decompress_inode(inode);
					mutex_unlock(&inode->i_mutex);
					if ( err < 0)
						return err;
				}
				ei->i_flags &= ~EXT2_DIRTY_FL;
				ei->i_compr_flags &= ~EXT2_CLEANUP_FL;
			}
			ei->i_compr_method = datum;
			inode->i_ctime = CURRENT_TIME;
			mark_inode_dirty(inode);
		}
#ifdef CONFIG_KMOD
		if (!ext2_algorithm_table[ext2_method_table[datum].alg].avail) {
			char str[32];

			sprintf(str, "ext2-compr-%s", ext2_algorithm_table[ext2_method_table[datum].alg].name);
			request_module(str);
		}
#endif
		datum = ((datum < EXT2_N_METHODS)
			 && (ext2_algorithm_table[ext2_method_table[datum].alg].avail));
		return put_user(datum, (long *)arg);

	case EXT2_IOC_GETCLUSTERBIT:
		if (get_user (datum, (long*) arg))
			return -EFAULT;
		if (!S_ISREG (inode->i_mode))
			return -ENOSYS;
		/* We don't do `down(&inode->i_sem)' here because
		   there's no way for userspace to do the
		   corresponding up().  Userspace must rely on
		   EXT2_NOCOMPR_FL if it needs to lock. */
		err = ext2_cluster_is_compressed (inode, datum);
		if (err < 0)
			return err;
		return put_user ((err ? 1 : 0),
				 (long *) arg);

	case EXT2_IOC_RECOGNIZE_COMPRESSED:
		if (get_user (datum, (long*) arg))
			return -EFAULT;
		if (!S_ISREG (inode->i_mode))
			return -ENOSYS;
		if (IS_RDONLY (inode))
			return -EROFS;
		return ext2_recognize_compressed (inode, datum);
		
	case EXT2_IOC_GETCLUSTERSIZE:
		/* Result means nothing if COMPR_FL is not set (until
                   SETCLUSTERSIZE w/o COMPR_FL is implemented;
                   todo). */
		if (!S_ISREG (inode->i_mode)
		    && !S_ISDIR (inode->i_mode)) 
			return -ENOSYS;
		return put_user (ei->i_clu_nblocks, (long *) arg);

	case EXT2_IOC_GETFIRSTCLUSTERSIZE:
		/* Result means nothing if COMPR_FL is not set (until
                   SETCLUSTERSIZE w/o COMPR_FL is implemented;
                   todo). */
		if (!S_ISREG (inode->i_mode)
		    && !S_ISDIR (inode->i_mode)) 
			return -ENOSYS;
		return put_user (ext2_first_cluster_nblocks(inode), (long *) arg);

	case EXT2_IOC_SETCLUSTERSIZE:
		if ((current_fsuid() != inode->i_uid) && !capable(CAP_FOWNER))
			return -EPERM;
		if (IS_RDONLY (inode))
			return -EROFS;
		if (get_user (datum, (long *) arg))
			return -EFAULT;
		if (!S_ISREG (inode->i_mode)
		    && !S_ISDIR (inode->i_mode)) 
			return -ENOSYS;

		/* These are the only possible cluster sizes.  The
		   cluster size must be a power of two so that
		   clusters don't straddle address (aka indirect)
		   blocks.  At the moment, the upper limit is constrained
		   by how much memory is allocated for de/compression.
		   Also, the gzip algorithms have some optimisations
		   that assume tht the input is no more than 32KB,
		   and in compress.c we would need to zero more bits
		   of head->holemap.  (In previous releases, the file
		   format was limited to 32 blocks and under 64KB.) */
// #if EXT2_MAX_CLUSTER_BLOCKS > 32 || EXT2_MAX_CLUSTER_NBYTES > 32768
// # error "This code not updated for cluster size yet."
// #endif
		switch (datum) {
		case (1 << 2): datum = 2; break;
		case (1 << 3): datum = 3; break;
		case (1 << 4): datum = 4; break;
		case (1 << 5): datum = 5; break;
		default: return -EINVAL;
		}

		assert (ei->i_clu_nblocks == (1 << ei->i_log2_clu_nblocks));
		if (datum == ei->i_log2_clu_nblocks)
			return 0;

		if (ei->i_flags & EXT2_ECOMPR_FL)
			return -EPERM;
		if (!(ei->i_flags & EXT2_COMPR_FL))
			return -ENOSYS;

		/* We currently lack a mechanism to change the cluster
		   size if there are already some compressed clusters.
		   The compression must be done in userspace
		   (e.g. with the e2compress program) instead.  */
		if (ei->i_flags & EXT2_COMPRBLK_FL)
			return -ENOSYS;

		if (datum + inode->i_sb->s_blocksize_bits
		    > EXT2_LOG2_MAX_CLUSTER_BYTES)
			return -EINVAL;

		ei->i_log2_clu_nblocks = datum;
		ei->i_clu_nblocks = 1 << datum;
		inode->i_ctime = CURRENT_TIME;
		mark_inode_dirty(inode);
		return 0;

	case EXT2_IOC_GETCOMPRRATIO:		
		if (!S_ISREG (inode->i_mode)) 
			return -ENOSYS;
		if (ei->i_flags & EXT2_ECOMPR_FL)
			return -EPERM;
		if ((long) (datum = ext2_count_blocks (inode)) < 0)
			return datum;
		if ((err = put_user ((long) datum, (long*) arg)))
			return err;
		return put_user ((long) inode->i_blocks, (long*) arg + 1);
		
		
#endif
	default:
		return -ENOTTY;
	}
}

#ifdef CONFIG_COMPAT
long ext2_compat_ioctl(struct file *file, unsigned int cmd, unsigned long arg)
{
	/* These are just misnamed, they actually get/put from/to user an int */
	switch (cmd) {
	case EXT2_IOC32_GETFLAGS:
		cmd = EXT2_IOC_GETFLAGS;
		break;
	case EXT2_IOC32_SETFLAGS:
		cmd = EXT2_IOC_SETFLAGS;
		break;
	case EXT2_IOC32_GETVERSION:
		cmd = EXT2_IOC_GETVERSION;
		break;
	case EXT2_IOC32_SETVERSION:
		cmd = EXT2_IOC_SETVERSION;
		break;
	default:
		return -ENOIOCTLCMD;
	}
	return ext2_ioctl(file, cmd, (unsigned long) compat_ptr(arg));
}
#endif
