/*
 *  linux/fs/ext2/compress.c
 *
 *  Copyright (C) 1995  Antoine Dumesnil de Maricourt (dumesnil@etca.fr) 
 *      (transparent compression code)
 */

/*
 *  Copyright (C) 2001 Alcatel Business Systems - R&D Illkirch FRANCE
 *
 *  	Transparent compression code for 2.4 kernel.
 *
 *  Denis Richard (denis.richard@sxb.bsf.alcatel.fr)
 *  Pierre Peiffer (pierre.peiffer@sxb.bsf.alcatel.fr)
 *
 *  Adapted from patch e2compr-0.4.39-patch-2.2.18 .
 */

#include <asm/segment.h>
#include <linux/errno.h>
#include <linux/fs.h>
#include <linux/ext2_fs.h>
#include <linux/ext2_fs_c.h>
#include <linux/fcntl.h>
#include <linux/sched.h>
#include <linux/stat.h>
#include <linux/buffer_head.h>
#include <linux/string.h>
#include <linux/kernel.h>
#include <linux/quotaops.h>
#include <linux/kmod.h>
#include <linux/vmalloc.h>
#include <linux/swap.h>
#include <linux/slab.h>
#include <linux/pagemap.h>
#include <linux/writeback.h>
#include <linux/rmap.h>
#include <linux/swap.h>
#include <linux/mm.h>
#include <linux/module.h>
#include <linux/slab.h>
#include <linux/kernel_stat.h>
#include <linux/swap.h>
#include <linux/pagemap.h>
#include <linux/init.h>
#include <linux/highmem.h>
#include <linux/vmstat.h>
#include <linux/file.h>
#include <linux/writeback.h>
#include <linux/blkdev.h>
#include <linux/buffer_head.h>
#include <linux/mm_inline.h>
#include <linux/pagevec.h>
#include <linux/backing-dev.h>
#include <linux/rmap.h>
#include <linux/topology.h>
#include <linux/cpu.h>
#include <linux/cpuset.h>
#include <linux/notifier.h>
#include <linux/rwsem.h>
#include <linux/delay.h>
#include <linux/kthread.h>
#include <linux/freezer.h>
#include <asm/tlbflush.h>
#include <asm/div64.h>
#include <linux/swapops.h>
#include <linux/percpu.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))

#ifdef CONFIG_HIGHMEM
#define restore_b_data_himem(bh)   assert(page_address(bh->b_page));  bh->b_data = page_address(bh->b_page) + bh_offset(bh)



int ext2_kmap_cluster_pages(struct page *page, struct page *pg[],
			    struct page *epg[])
{
    int i = 0;

    for (i = 0; i < EXT2_MAX_CLUSTER_PAGES; i++) {
	if (!pg[i])
	    break;
	if (epg && epg[i])
	    kmap(epg[i]);
	else
	    kmap(pg[i]);
    }

    if (page)
	kmap(page);
    return 0;
}


int ext2_kunmap_cluster_pages(struct page *page, struct page *pg[],
			      struct page *epg[])
{
    int i = 0;

    for (i = 0; i < EXT2_MAX_CLUSTER_PAGES; i++) {
	if (!pg[i])
	    break;
	if (epg && epg[i])
	    kunmap(epg[i]);
	else
	    kunmap(pg[i]);
    }

    if (page)
	kunmap(page);
    return 0;
}
#else //no high-mem:
#define restore_b_data_himem(bh)	;
#endif


/*none compression dummy functions*/
size_t ext2_iNONE (int action) { return 0; }
size_t ext2_wNONE (__u8 *ibuf, __u8 *obuf, void *wa, size_t ilen, size_t olen, int xarg) { return 0; }
size_t ext2_rNONE (__u8 *ibuf, __u8 *obuf, void *wa, size_t ilen, size_t olen, int xarg) { return 0; }

/*
 *    Algorithm and method tables
 */
struct ext2_algorithm ext2_algorithm_table[] = {
    /* Note: all algorithms must have the `name' field filled in.
       This is used to autoload algorithm modules (ext2-compr-%s), and
       in kernel printk. */
    /* N.B. Do not renumber these algorithms!  (To do so is to change
       the binary format.)  It's OK for `none' and `undef' to be
       renumbered, though. */

    /* Fields:
       name; available; routines for:
       init,  compress,   decompress. */
    {"lzv1", 0, ext2_iNONE, ext2_wNONE, ext2_rNONE},
    {"lzrw3a", 0, ext2_iNONE, ext2_wNONE, ext2_rNONE},
    {"gzip", 1, ext2_iZLIB, ext2_wZLIB, ext2_rZLIB},	//Andreas: workaround
    {"bzip2", 0, ext2_iNONE, ext2_wNONE, ext2_rNONE},
    {"lzo", 0, ext2_iNONE, ext2_wNONE, ext2_rNONE},
    {"none", 1, ext2_iNONE, ext2_wNONE, ext2_rNONE},

    /* This "algorithm" is for unused entries in the method table.
       It differs from EXT2_NONE_ALG in that it is considered
       unavailable, whereas `none' is always available. */
    {"undef", 0, ext2_iNONE, ext2_wNONE, ext2_rNONE},

};

/* Note: EXT2_N_ALGORITHMS can't be increased beyond 16 without
   changing the width of the s_algorithms_used field in the in-memory
   superblock.  The on-disk s_algorithms_used field is 32 bits long.
   (This is in a state of flux.  Currently (1998-02-05) there is no
   distinction: we always use the s_es copy. */

/* The size of this table must be 32 to prevent Oopsen from
   invalid data.  We index this from 5 bits of i_flags, so
   the size is (1 << 5) == 32. */
struct ext2_method ext2_method_table[32] = {
    /* Fields: algorithm id, algorithm argument. */
    {EXT2_LZV1_ALG, 0},
    {EXT2_NONE_ALG, 0},		/* 1: auto */
    {EXT2_NONE_ALG, 0},		/* 2: defer */
    {EXT2_NONE_ALG, 0},		/* 3: never */
    {EXT2_BZIP2_ALG, 0},	/* 4: bzip2 */
    {EXT2_UNDEF_ALG, 0},
    {EXT2_UNDEF_ALG, 0},
    {EXT2_UNDEF_ALG, 0},
    {EXT2_LZRW3A_ALG, 0},	/* 8: lzrw3a */
    {EXT2_UNDEF_ALG, 0},
    {EXT2_LZO_ALG, 0},		/* 10: lzo1x_1 */
    {EXT2_UNDEF_ALG, 0},
    {EXT2_UNDEF_ALG, 0},
    {EXT2_UNDEF_ALG, 0},
    {EXT2_UNDEF_ALG, 0},
    {EXT2_UNDEF_ALG, 0},
    {EXT2_GZIP_ALG, 1},		/* 16 */
    {EXT2_GZIP_ALG, 2},
    {EXT2_GZIP_ALG, 3},
    {EXT2_GZIP_ALG, 4},
    {EXT2_GZIP_ALG, 5},
    {EXT2_GZIP_ALG, 6},
    {EXT2_GZIP_ALG, 7},
    {EXT2_GZIP_ALG, 8},
    {EXT2_GZIP_ALG, 9},
    {EXT2_UNDEF_ALG, 0},
    {EXT2_UNDEF_ALG, 0},
    {EXT2_UNDEF_ALG, 0},
    {EXT2_UNDEF_ALG, 0},
    {EXT2_UNDEF_ALG, 0},
    {EXT2_UNDEF_ALG, 0},
    {EXT2_UNDEF_ALG, 0}
};


static void ext2_mark_algorithm_use(struct inode *inode, unsigned alg)
{
    struct ext2_sb_info *sbi = EXT2_SB(inode->i_sb);

    /* Hopefully, lock_super() isn't needed here, as we don't
       block in the critical region.  True? */
    assert(alg < EXT2_N_ALGORITHMS);
    if (sbi->s_es->s_feature_incompat
	& cpu_to_le32(EXT2_FEATURE_INCOMPAT_COMPRESSION)) {
	sbi->s_es->s_algorithm_usage_bitmap |= cpu_to_le32(1 << alg);
    } else {
	struct ext2_super_block *es = sbi->s_es;

	es->s_algorithm_usage_bitmap = cpu_to_le32(1 << alg);
	es->s_feature_incompat
	    |= cpu_to_le32(EXT2_FEATURE_INCOMPAT_COMPRESSION);
	if (es->s_rev_level < EXT2_DYNAMIC_REV) {
	    /* Raise the filesystem revision level to
	       EXT2_DYNAMIC_REV so that s_feature_incompat
	       is honoured (except in ancient kernels /
	       e2fsprogs).  We must also initialize two
	       other dynamic-rev fields.  The remaining
	       fields are assumed to be already correct
	       (e.g. still zeroed). */
	    es->s_rev_level = cpu_to_le32(EXT2_DYNAMIC_REV);
	    es->s_first_ino = cpu_to_le32(EXT2_GOOD_OLD_FIRST_INO);
	    es->s_inode_size = cpu_to_le16(EXT2_GOOD_OLD_INODE_SIZE);
	}
    }
    mark_buffer_dirty(sbi->s_sbh);
}


/* Displays an error message if algorithm ,alg` is not marked in use,
   and then marks it in use. */
static void ext2_ensure_algorithm_use(struct inode *inode, unsigned alg)
{
    assert(alg < EXT2_N_ALGORITHMS);

    if (!(EXT2_SB(inode->i_sb)->s_es->s_algorithm_usage_bitmap
	  & cpu_to_le32(1 << alg))) {
	ext2_msg(inode->i_sb, "algorithm usage bitmap algorithm %s not marked used in inode %lu",
		     ext2_algorithm_table[alg].name, inode->i_ino);
	ext2_mark_algorithm_use(inode, alg);
    }
}


/*mw: out of cache bug fix 5-16-07 */
static void create_empty_buffers_e2c(struct page *page,
				     unsigned long blocksize,
				     unsigned long b_state,
				     struct inode *inode)
{
    struct buffer_head *bh, *head, *tail;

    head = alloc_page_buffers(page, blocksize, 1);
    bh = head;
    do {
	bh->b_state |= b_state;
	tail = bh;
	bh->b_bdev = NULL;	//mw: make it like 2.4
	bh->b_blocknr = 0;	//mw: make it like 2.4
	bh->b_end_io = NULL;	//mw: make it like 2.4
	bh = bh->b_this_page;
    } while (bh);
    tail->b_this_page = head;
    spin_lock(&inode->i_mapping->private_lock);
    if (PageUptodate(page) || PageDirty(page)) {
	bh = head;
	do {
	    if (PageDirty(page))
		set_buffer_dirty(bh);
	    if (PageUptodate(page))
		set_buffer_uptodate(bh);
	    bh = bh->b_this_page;
	} while (bh != head);
    }
    attach_page_buffers(page, head);
    spin_unlock(&inode->i_mapping->private_lock);
}

int ext2_get_cluster_pages(struct inode *inode, u32 cluster,
			   struct page *pg[], struct page *page, int compr)
{
    int nbpg, npg, i;
    u32 page0;	/* = position within file (not position within fs). */
    u32 idx = 0;

    /*mw */
    for (i = 0; i < EXT2_MAX_CLUSTER_PAGES; i++)
	pg[i] = NULL;

    page0 = ext2_cluster_page0(inode, cluster);	
    nbpg = ext2_cluster_npages(inode, cluster);	

    if (compr && (((page0 + nbpg) << PAGE_CACHE_SHIFT) > inode->i_size))	
	nbpg = ((inode->i_size - 1) >> PAGE_CACHE_SHIFT) - page0 + 1;	
#ifdef  EXT2_COMPR_REPORT
    trace_e2c("ext2_get_cluster_pages: page0=%d, nbpg=%d page=%ld\n",
	      page0, nbpg, ((page != NULL) ? page->index : 0));
#endif
    for (npg = 0; npg < nbpg; npg++) {
	if ((page == NULL) || ((page0 + npg) != page->index)) {
		//pg[npg] = __grab_cache_page(inode->i_mapping,  page0+npg); /* &cached_page, &lru_pvec);*/
		pg[npg] = grab_cache_page_write_begin(inode->i_mapping, page0+npg,  0); 
	    if (!pg[npg])
		goto error;
	} else {
	    pg[npg] = page;
	}
	if (!page_has_buffers(pg[npg])) {
	    ClearPageUptodate(pg[npg]);	
	    ClearPageDirty(pg[npg]);	
	    create_empty_buffers_e2c(pg[npg], inode->i_sb->s_blocksize, 0, inode);	
	    if (unlikely(!page_has_buffers(pg[npg])))
		trace_e2c("ext2_get_cluster_pages: NOMEM!\n");
	    assert(!PageUptodate(pg[npg]));
	    assert(!PageDirty(pg[npg]));
	}
    }
    //set remaining pages to NULL
    for (idx = npg; idx < EXT2_MAX_CLUSTER_PAGES; idx++)
	pg[idx] = NULL;

    return (npg);

  error:
    while (--npg >= 0) {
	if ((page == NULL) || ((page0 + npg) != page->index)) {
	    unlock_page(pg[npg]);
	    page_cache_release(pg[npg]);
	}
	pg[npg] = NULL;
    }
    trace_e2c("ext2_get_cluster_pages: error no page\n");
    return (-ENOMEM);
}


int ext2_get_cluster_extra_pages(struct inode *inode, u32 cluster,
				 struct page *pg[], struct page *epg[])
{
    struct page *page;
    int nbpg, npg, i;

    for (i = 0; i < EXT2_MAX_CLUSTER_PAGES; i++)
	epg[i] = NULL;

    nbpg = ext2_cluster_npages(inode, cluster);
    for (npg = 0; npg < nbpg; npg++) {
	if (pg[npg] == NULL)
	    break;
	if (PageUptodate(pg[npg])) {
	    //page = page_cache_alloc(inode->i_mapping);   
	    //mw: has gfp-mask of adress-space:  gfp_t mapping_gfp_mask(struct address_space * mapping)
	    //    don't trigger. shrink_dcache_memory which might call ext2_cleanup_compressed_inode with the SAME mutex.
	    page = __page_cache_alloc(GFP_NOFS);

	    if (!page) {
		goto error;
	    }
	    ClearPageError(page);
	    ClearPageReferenced(page);
	    ClearPageUptodate(page);
	    ClearPageDirty(page);
	    lock_page(page);
	    page->index = pg[npg]->index;

	    if (!page_has_buffers(page)) {
		create_empty_buffers_e2c(page, inode->i_sb->s_blocksize, 0,
					 inode);
		/*mw : only the "extra_pages" for decompression need create_empty_buffers_unlocked, because
		 *     they have no mapping-context and they must not have one. Otherwise they get need a page->index
		 *     which belongs always to an address_space object (e.g.: inode). But I think this is not intented here.
		 *     we just need thei buffers for a short time of decompression */
		if (unlikely(!page_has_buffers(page)))
		    return printk("Error: NOMEM!\n");
	    }

	    epg[npg] = page;
#ifdef  EXT2_COMPR_REPORT
	    trace_e2c
		("ext2_get_cluster_extra_pages: allocated page idx=%ld\n",
		 pg[npg]->index);
#endif
	} else {
	    epg[npg] = NULL;
	}
    }
    return (npg);
  error:
    while (--npg >= 0)
	if (epg[npg]) {
	    ClearPageDirty(epg[npg]);
	    ClearPageUptodate(epg[npg]);
	    try_to_free_buffers(epg[npg]);
	    unlock_page(epg[npg]);
	    assert(page_count(epg[npg]) == 1);
	    page_cache_release(epg[npg]);
	}
    trace_e2c("ext2_get_cluster_extra_pages: error no page\n");
    return (-ENOMEM);

}

/* Read every block in the cluster.  The blocks are stored in the bh
   array, which must be big enough.

   Return the number of block contained in the cluster, or -errno if an
   error occured.  The buffers should be released by the caller
   (unless an error occurred).
 
   The inode must be locked, otherwise it is possible that we return
   some out of date blocks.
 
   Called by :
 
         ext2_decompress_cluster()      [i_sem]
         ext2_compress_cluster()        [i_sem]
         ext2_readpage()      		[i_sem] */


int ext2_get_cluster_blocks(struct inode *inode, u32 cluster,
			    struct buffer_head *bh[], struct page *pg[],
			    struct page *epg[], int compr)
{
    struct buffer_head *br[EXT2_MAX_CLUSTER_BLOCKS];
    int nreq, nbh = 0, npg, i;
    u32 clu_nblocks;
    int err;
    const int blocks = PAGE_CACHE_SIZE >> inode->i_sb->s_blocksize_bits;

    /*mw */
    for (i = 0; i < EXT2_MAX_CLUSTER_BLOCKS; i++)
	bh[i] = NULL;

    assert(atomic_read(&inode->i_mutex.count) <= 0);	/* i.e. mutex_lock */

    /*
     *  Request full cluster.
     */
    {
	u32 endblk;
	u32 block; /* = position within file (not position within fs). */
	u32 nbpg;
	u32 page0; /* = position within file (not position within fs). */
	u32 idx;

	block = ext2_cluster_block0(inode, cluster);
	clu_nblocks = ext2_cluster_nblocks(inode, cluster);
	/* impl: Don't shorten endblk for i_size.  The
	   remaining blocks should be NULL anyway, except in
	   the case when called from ext2_decompress_cluster
	   from ext2_truncate, in which case i_size is short
	   and we _want_ to get all of the blocks. */
	endblk = block + clu_nblocks;

	page0 = ext2_cluster_page0(inode, cluster);
	nbpg = ext2_cluster_npages(inode, cluster);

	if (compr
	    && (((page0 + nbpg) << PAGE_CACHE_SHIFT) > inode->i_size)) {
	    nbpg = ((inode->i_size - 1) >> PAGE_CACHE_SHIFT) - page0 + 1;
	    endblk =
		block +
		(nbpg <<
		 (PAGE_CACHE_SHIFT - inode->i_sb->s_blocksize_bits));
	}

	idx = page0 << (PAGE_CACHE_SHIFT - inode->i_sb->s_blocksize_bits);
#ifdef  EXT2_COMPR_REPORT
	trace_e2c("ext2_get_cluster_blocks: page0=%d, nbpg=%d\n", page0,
		  nbpg);
#endif
	for (npg = 0; npg < nbpg; npg++) {
	    struct buffer_head *buffer;

	    if ((epg != NULL) && (epg[npg] != NULL))
		buffer = page_buffers(epg[npg]);
	    else
		buffer = page_buffers(pg[npg]);
	    for (i = 0; i < blocks && (block + nbh) < endblk;
		 buffer = buffer->b_this_page, i++) {
		if (idx == (block + nbh)) {
		    bh[nbh] = buffer;
		    nbh++;
		}
		idx++;
	    }
	}
#ifdef  EXT2_COMPR_REPORT
	trace_e2c
	    ("ext2_get_cluster_blocks: get every pages and %d buffers\n",
	     nbh);
#endif

	for (nbh = 0, nreq = 0; block < endblk; nbh++) {
	    assert(bh[nbh] != NULL);
	    bh[nbh]->b_blocknr = 0;
	    clear_bit(BH_Mapped, &bh[nbh]->b_state);

	    //mw: does not work with 2.6 and holes!!!
	    //err=ext2_get_block(inode, block++, bh[nbh], (PageDirty(bh[nbh]->b_page) ? 1 : 0));  
	    err = ext2_get_block(inode, block++, bh[nbh], 0);
	    /* mw: 0: we dont' create non existing blocks here
	     *     let's do it just before the writeback, when we know, which blocks we really need...*/
	    //err=ext2_get_block(inode, block++, bh[nbh], (buffer_dirty(bh[nbh]) ? 1 : 0));

	    /* mw: bdev-bug-fix: for files which got compressed and now consume less buffers
	     * ext2_get_block returns 0, for a empty-block. As these buffer were used before
	     * the bh[nbh]->b_bdev might be != NULL or just invalid. So we set them explicitly
	     * to NULL. */
	    //printk("Get Block cluster %i: (%#x):%i Blk-NR:%lu(%lu)[%lu-%lu] Bdev:%#x(%#x), PGDirty:%i, mapped:%i, PID: %lu\n", cluster, bh[nbh], nbh, block, 

	    //if we are not mapped, then the blocknr will be wrong
	    //we set a bdev here the we will write to some "random" block
	    if (!buffer_mapped(bh[nbh])) {
		bh[nbh]->b_bdev = NULL;	/* don't write wrongly mapped blocks !!! */
		/* mw: you encounter null pointer oops you MUST
		 *         map your buffer using ext2_get_block()*/
	    }

	    if (bh[nbh]->b_blocknr != 0) {
		if (!buffer_uptodate(bh[nbh])
		    /* TODO: Do we need this
		       `!buffer_locked' test? */
		    && !buffer_locked(bh[nbh])
		    && !PageDirty(bh[nbh]->b_page))
		    br[nreq++] = bh[nbh];
	    } else if ((err != 0)
		       && (err != -EFBIG))
		/* impl: for some unknown reason,
		   ext2_getblk() returns -EFBIG if
		   !create and there's a hole. ==> not right any more in 2.4 */
		goto error;
	}
	for (i = nbh; i < EXT2_MAX_CLUSTER_BLOCKS; i++) {
	    bh[i] = NULL;
	}
    }
#ifdef  EXT2_COMPR_REPORT_CPR
    trace_e2c("ext2_get_cluster_blocks: nreq=%d for cluster=%d\n", nreq,
	      cluster);
#endif

    //read all blocks, which are not null-blocks
    if (nreq > 0)
	ll_rw_block(READ, nreq, br);

    /*
     *  Adjust nbh if we have some null blocks at end of cluster.
     */
    while ((nbh != 0) && (bh[nbh - 1]->b_blocknr == 0))
	nbh--;

    /*
     *  Wait for blocks.
     */
    err = -EIO;
    CHECK_NOT_ATOMIC
    for (i = 0; i < nbh; i++)
	if ((!PageDirty(bh[i]->b_page)) && (bh[i]->b_blocknr != 0)) {  
	    wait_on_buffer(bh[i]);
	    if (!buffer_uptodate(bh[i])) {	/* Read error ??? */
		trace_e2c
		    ("ext2_get_cluster_blocks: wait_on_buffer error (blocknr=%ld)\n",
		     bh[i]->b_blocknr);
		goto error;
	    }
	}
    assert(nbh <= EXT2_MAX_CLU_NBLOCKS);

    return nbh;

  error:
    printk("ERROR: ext2_get_cluster_blocks()\n");
    return err;
}


/* Iterations over block in the inode are done with a generic
   iteration key mechanism.  We need one method to convert a block
   number into a new key, one method to iterate (i.e., increment the
   key) and one method to free the key.  The code could be shared with
   truncate.c, as this mechanism is very general.
 
   This code assumes tht nobody else can read or write the file
   between ext2_get_key() and ext2_free_key(), so callers need to have
   i_sem (which they all do anyway). */

/* TODO: Get all of the bkey routines to return -errno instead of
   true/false. */
/* TODO: The bkey routines currently assume tht address blocks are
   allocated even if all contained addresses are NULL, but this is not
   true.  Make sure tht we differentiate between NULL block and error,
   and then fix up ext2_set_key_blkaddr() and anything else (including
   the pack/unpack routines). */
struct ext2_bkey {
    int level;
    u32 block;
    struct inode *inode;
    int off[4];
    u32 *ptr[4];
    struct buffer_head *ibh[4];
};


/*
 *    Method to convert a block number into a key.
 *
 *    Returns 1 on success, 0 on failure.  You may safely, but need
 *    not, free the key even if ext2_get_key() fails. 
 */
static int ext2_get_key(struct ext2_bkey *key, struct inode *inode,
			u32 block)
{
    int x, level;
    int addr_per_block = EXT2_ADDR_PER_BLOCK(inode->i_sb);

    assert(atomic_read(&inode->i_mutex.count) <= 0);

    /*
     *  The first step can be viewed as translating the
     *    original block number in a special base (powers
     *    of addr_per_block).
     */

    key->block = block;

    key->off[0] = key->off[1] = key->off[2] = key->off[3] = 0;
    key->ibh[0] = key->ibh[1] = key->ibh[2] = key->ibh[3] = NULL;
    key->ptr[0] = key->ptr[1] = key->ptr[2] = key->ptr[3] = NULL;

    if (block >= EXT2_NDIR_BLOCKS) {
	block -= EXT2_NDIR_BLOCKS;

	if (block >= addr_per_block) {
	    block -= addr_per_block;

	    if (block >= addr_per_block * addr_per_block) {
		block -= addr_per_block * addr_per_block;

		key->off[0] = EXT2_TIND_BLOCK;
		key->off[1] = (block / (addr_per_block * addr_per_block));
		key->off[2] =
		    (block % (addr_per_block * addr_per_block)) /
		    addr_per_block;
		key->off[3] = (block % addr_per_block);
		level = 3;
	    } else {
		key->off[0] = EXT2_DIND_BLOCK;
		key->off[1] = block / addr_per_block;
		key->off[2] = block % addr_per_block;
		level = 2;
	    }
	} else {
	    key->off[0] = EXT2_IND_BLOCK;
	    key->off[1] = block;
	    level = 1;
	}
    } else {
	key->off[0] = block;
	level = 0;
    }

    /*
     *  In the second step, we load the needed buffers.
     */

    key->level = level;
    key->inode = inode;

    key->ptr[0] = (u32 *) (&(EXT2_I(inode)->i_data));

    for (x = 1; x <= level; x++) {
	u32 *ptr;

	ptr = key->ptr[x - 1];
	if (ptr == NULL)
	    break;
/* Paul Whittaker tweak 19 Feb 2005 */
	block = le32_to_cpu(ptr[key->off[x - 1]]);
	if (block == 0)
	    continue;		// TLL 05/01/07
	if (x - 1 != 0)
	    block = le32_to_cpu(block);
	if ((key->ibh[x] = __bread(inode->i_sb->s_bdev,
				   block, inode->i_sb->s_blocksize))
	    == NULL)
	    goto error;
	key->ptr[x] = (u32 *) (key->ibh[x]->b_data);
    }

    return 1;
  error:
    for (; x != 0; x--)
	if (key->ibh[x] != NULL)
	    brelse(key->ibh[x]);
    return 0;
}


/*
 *    Find the block for a given key.  Return 0 if there
 *      is no block for this key.
 */
static inline u32 ext2_get_key_blkaddr(struct ext2_bkey *key)
{
    assert(key->inode);
    assert(atomic_read(&(key->inode)->i_mutex.count) <= 0);

/* Paul Whittaker tweak 19 Feb 2005 */
    if (key->ptr[key->level] == NULL)
	return 0;
    return le32_to_cpu(key->ptr[key->level][key->off[key->level]]);
}


/*
 *    Change the block for a given key.  Return 0 on success,
 *      -errno on failure.
 */
static inline int ext2_set_key_blkaddr(struct ext2_bkey *key, u32 blkaddr)
{
    char bdn[BDEVNAME_SIZE];
    assert(key->inode);
    assert(atomic_read(&(key->inode)->i_mutex.count) <= 0);

    if (key->ptr[key->level] == NULL) {
	/* The reason that this "can't happen" is that this
	   routine is only used to shuffle block numbers or by
	   free_cluster_blocks.  Cluster sizes are such that
	   clusters can't straddle address blocks.  So the
	   indirect block address can't be zero.  AFAIK, ptr
	   can only be NULL on error or on null indirect block
	   address.  Hmm, come to think of it, I think there
	   are still some callers that don't check for errors
	   from ext2_get_key(), so this still can happen until
	   those are fixed up. */
	printk(KERN_ERR
	       "ext2_set_key_blkaddr: can't happen: NULL parent.  "
	       "dev=%s, ino=%lu, level=%u.\n",
	       bdevname(key->inode->i_sb->s_bdev, bdn),
	       key->inode->i_ino, key->level);
	return -ENOSYS;
    }
    /* Paul Whittaker tweak 19 Feb 2005 */
    key->ptr[key->level][key->off[key->level]] = le32_to_cpu(blkaddr);
    if (key->level > 0)
	mark_buffer_dirty(key->ibh[key->level]);
    return 0;
}


/*
 *    Increment the key.  Returns 0 if we go beyond the limits,
 *      1 otherwise.
 *
 *    Precondition: -key->off[level] <= incr < addr_per_block.
 */
static int ext2_next_key(struct ext2_bkey *key, int incr)
{
    int addr_per_block = EXT2_ADDR_PER_BLOCK(key->inode->i_sb);
    int x, level = key->level;
    u32 tmp;

    assert(key->inode);
    assert(atomic_read(&(key->inode)->i_mutex.count) <= 0);


    /*
     *  Increment the key. This is done in two step: first
     *    adjust the off array, then reload buffers that should
     *    be reloaded (we assume level > 0).
     */

    assert(key->off[level] >= -incr);
    assert(incr < addr_per_block);
    key->block += incr;
    key->off[level] += incr;

    /*
     *  First step: should be thought as the propagation
     *    of a carry.
     */

    if (level == 0) {
	if (key->off[0] >= EXT2_NDIR_BLOCKS) {
	    key->off[1] = key->off[0] - EXT2_NDIR_BLOCKS;
	    key->off[0] = EXT2_IND_BLOCK;
	    level = 1;
	}
	x = 0;
    } else {
	for (x = level; x > 0; x--) {
	    if (key->off[x] >= addr_per_block) {
		key->off[x] -= addr_per_block;
		key->off[x - 1]++;

		if (x == 1) {
		    if (++level < 4) {
			key->off[level] = key->off[level - 1];
			key->off[level - 1] = 0;
		    } else
			return 0;
		}
	    } else
		break;
	}
    }

    /*
     *  Second step: reload the buffers that have changed.
     */

    key->level = level;

    CHECK_NOT_ATOMIC
    while (x++ < level) {
	if (key->ibh[x] != NULL) {
	    if (IS_SYNC(key->inode) && buffer_dirty(key->ibh[x])) {
		//mw:
		assert(buffer_mapped(key->ibh[x])
		       && (key->ibh[x]->b_bdev != NULL));
		ll_rw_block(WRITE, 1, &(key->ibh[x]));
		wait_on_buffer(key->ibh[x]);
	    }
	    brelse(key->ibh[x]);
	}
/* Paul Whittaker tweak 19 Feb 2005 */
	if ((key->ptr[x - 1] != NULL)
	    && ((tmp = le32_to_cpu(key->ptr[x - 1][key->off[x - 1]])) !=
		0)) {
	    if ((key->ibh[x] =
		 __bread(key->inode->i_sb->s_bdev, tmp,
			 key->inode->i_sb->s_blocksize))
		!= NULL)
		key->ptr[x] = (u32 *) (key->ibh[x]->b_data);
	    else
		key->ptr[x] = NULL;
	} else {
	    key->ibh[x] = NULL;
	    key->ptr[x] = NULL;
	}
    }

    return 1;
}


/* Method to free the key: just release buffers.

   Returns 0 on success, -errno on error.
*/

static int ext2_free_key(struct ext2_bkey *key)
{
    int x, n;
    struct buffer_head *bh[4];

    assert(key->inode);
    assert(atomic_read(&(key->inode)->i_mutex.count) <= 0);


    for (x = 0, n = 0; x <= key->level; x++) {
	if (key->ibh[x] != NULL) {
	    if (IS_SYNC(key->inode) && buffer_dirty(key->ibh[x]))
		bh[n++] = key->ibh[x];
	    else
		brelse(key->ibh[x]);
	}
    }

    if (n > 0) {
	int ncopy = n;
	while (ncopy-- > 0) {
	    assert(buffer_mapped(bh[ncopy])
		   && (bh[ncopy]->b_bdev != NULL));
	}

	ll_rw_block(WRITE, n, bh);

	CHECK_NOT_ATOMIC

	while (n-- > 0) {
	    wait_on_buffer(bh[n]);
	    /* TODO: Check for error. */
	    brelse(bh[n]);
	}
    }
    return 0;
}


/* Returns positive if specified cluster is compressed,
   zero if not,
   -errno if an error occurred.

   If you need the result to be accurate, then down i_sem before
   calling this, and don't raise i_sem until after you've used the
   result. */
int ext2_cluster_is_compressed_fn(struct inode *inode, unsigned cluster)
{
    unsigned block = (ext2_cluster_block0(inode, cluster)
		      + ext2_cluster_nblocks(inode, cluster)
		      - 1);
    struct ext2_bkey key;
    int result;

    assert(atomic_read(&inode->i_mutex.count) <= 0);

    /* impl: Not all callers of ext2_cluster_is_compressed_fn() have
       i_sem down.  Of course it is impossible to guarantee
       up-to-date information for such callers (someone may
       compress or decompress between when we check and when they
       use the information), so hopefully it won't matter if the
       information we return is slightly inaccurate (e.g. because
       someone is de/compressing the cluster while we check). */
    if (!ext2_get_key(&key, inode, block))
	return -EIO;

    result = (ext2_get_key_blkaddr(&key) == EXT2_COMPRESSED_BLKADDR);
    ext2_free_key(&key);
    return result;
}


/* Support for the GETCOMPRRATIO ioctl() call.  We calculate how many
   blocks the file would hold if it weren't compressed.  This requires
   reading the cluster head for every compressed cluster.

   Returns either -EAGAIN or the number of blocks that the file would
   take up if uncompressed.  */
int ext2_count_blocks(struct inode *inode)
{
    struct buffer_head *head_bh;
    int count;
    int cluster;
    struct ext2_bkey key;
    u32 end_blknr;

    if (!(EXT2_I(inode)->i_flags & EXT2_COMPRBLK_FL))
	return inode->i_blocks;

    mutex_lock(&inode->i_mutex);
    end_blknr = ROUNDUP_RSHIFT(inode->i_size,
			       inode->i_sb->s_blocksize_bits);

    /* inode->i_blocks is stored in units of 512-byte blocks.  It's
       more convenient for us to work in units of s_blocksize. */
    {
	u32 shift = inode->i_sb->s_blocksize_bits - 9;

	count = inode->i_blocks;
	if (count & ((1 << shift) - 1))
	    ext2_msg(inode->i_sb,
			 "ext2_count_blocks",
			 "i_blocks not multiple of blocksize");
	count >>= shift;
    }

    cluster = 0;
    if (!ext2_get_key(&key, inode, 0)) {
	count = -EIO;
	goto out;
    }
    while (key.block < end_blknr) {
	u32 head_blkaddr = ext2_get_key_blkaddr(&key);

	/* bug fix: init head_bh for each iteration TLL 2/21/07 */
	head_bh = NULL;
	if (head_blkaddr == EXT2_COMPRESSED_BLKADDR) {
	    count = -EXT2_ECOMPR;
	    break;
	}
	if (!ext2_next_key(&key, ext2_cluster_nblocks(inode, cluster) - 1))
	    break;
	if (ext2_get_key_blkaddr(&key) == EXT2_COMPRESSED_BLKADDR) {
	    struct ext2_cluster_head *head;

	    if (head_blkaddr == 0) {
		count = -EXT2_ECOMPR;
		break;
	    }
	    head_bh = __getblk(inode->i_sb->s_bdev,
			       head_blkaddr, inode->i_sb->s_blocksize);
	    if (head_bh == NULL) {
		/* Hmm, EAGAIN or EIO? */
		count = -EAGAIN;
		break;
	    }
	    if (!buffer_uptodate(head_bh))
		ll_rw_block(READ, 1, &head_bh);

	    CHECK_NOT_ATOMIC

	    wait_on_buffer(head_bh);

#ifdef CONFIG_HIGHMEM
	    if (!page_address(head_bh->b_page)) {
		BUG();
	    }
#endif

	    head = (struct ext2_cluster_head *) head_bh->b_data;
	    /* remove clen > ulen test TLL 2/21/07 */
	    if ((head->magic != cpu_to_le16(EXT2_COMPRESS_MAGIC_04X))
		|| (le32_to_cpu(head->ulen) > EXT2_MAX_CLUSTER_BYTES)
		|| (head->holemap_nbytes > 4)) {
		count = -EXT2_ECOMPR;
		break;
	    }
	    assert(sizeof(struct ext2_cluster_head) == 16);
	    count += (ROUNDUP_RSHIFT(le32_to_cpu(head->ulen),
				     inode->i_sb->s_blocksize_bits)
		      - ROUNDUP_RSHIFT((le32_to_cpu(head->clen)
					+ sizeof(struct ext2_cluster_head)
					+ head->holemap_nbytes),
				       inode->i_sb->s_blocksize_bits));
	    brelse(head_bh);
	    head_bh = NULL;
	}

	if (!ext2_next_key(&key, 1))
	    break;
	cluster++;
    }
    ext2_free_key(&key);
    if (head_bh != NULL)
	brelse(head_bh);
  out:
    mutex_unlock(&inode->i_mutex);
    if (count == -EXT2_ECOMPR) {
	ext2_msg(inode->i_sb,
		     "ext2_count_blocks",
		     "invalid compressed cluster %u of inode %lu",
		     cluster, inode->i_ino);
	EXT2_I(inode)->i_flags |= EXT2_ECOMPR_FL;
    }

    /* The count should be in units of 512 (i.e. 1 << 9) bytes. */
    if (count >= 0)
	count <<= inode->i_sb->s_blocksize_bits - 9;
    return count;
}


/* Decompress some blocks previously obtained from a cluster.
   Decompressed data is stored in ext2_rd_wa.u.  Buffer heads in the bh
   array are packed together at the begining of the array.  The ulen
   argument is an indication of how many bytes the caller wants to
   obtain, excluding holes.  (This can be less than head->ulen, as in the
   case of readpage.)  No hole processing is done; we don't even look at
   head->holemap.

   Note the semantic difference between this and
   (): the latter decompresses a cluster _and
   stores it as such_, whereas ext2_decompress_blocks() just
   decompresses the contents of the blocks into ext2_rd_wa.u.

   The working area is supposed to be available and locked.
 
   Returns a negative value on failure, the number of bytes
   decompressed otherwise.
 
   Called by :
 
         ext2_decompress_cluster ()    [sem down]
         ext2_readpage () [sem down, but only ifndef EXT2_LOCK_BUFFERS] */
	
/* TODO: ext2_decompress_blocks() scribbles in ext2_rd_wa.c.
   Check callers to make sure this isn't a problem. */

/* mw: caller must already have done: "get_cpu_var(ext2_rd_wa)" */
size_t
ext2_decompress_blocks(struct inode * inode,
		       struct buffer_head ** bh,
		       int nblk, size_t ulen, u32 cluster)
{
    struct ext2_cluster_head *head;
    int count, src_ix, x;
    unsigned char *dst;
    unsigned meth, alg;
    char bdn[BDEVNAME_SIZE];

#ifdef EXT2_COMPR_DEBUG
    //mw: 30.04.2012: seems to fail... ? assert(in_atomic());
    assert(atomic_read(&inode->i_mutex.count) <= 0);	/* i.e. mutex_lock */
#endif

    /*
       We pack the buffer together before (and must take care
       not to duplicate the buffer heads in the array).

       pjm 1998-01-09: Starting from e2compr-0.4.0, they should
       already be packed together in the blkaddr array.  TODO:
       Insert appropriate assert() statements checking tht this is
       the case.  TODO: Check that callers have bh[] packed. */
#ifdef  EXT2_COMPR_REPORT
    trace_e2c("ext2_decompress_blocks: nblk=%d\n", nblk);
#endif
    for (src_ix = 0, x = 0; src_ix < nblk; src_ix++) {
	if (bh[src_ix] == NULL)
	    printk("no_bheader()\n");
	if ((bh[src_ix] != NULL) && (bh[src_ix]->b_blocknr != 0)) {

	    if (x < src_ix) {
		ext2_msg(inode->i_sb, "bad buffer table",
			     "inode = %lu", inode->i_ino);
		goto error;
	    }
	    x++;
	}
    }

    nblk = x;
#ifdef  EXT2_COMPR_REPORT_CPR
    trace_e2c("ext2_decompress_blocks (2): nblk=%d\n", nblk);
#endif
    if (nblk == 0) {
	ext2_msg(inode->i_sb, "no block in cluster", "inode = %lu",
		     inode->i_ino);
	goto error;
    }

    restore_b_data_himem(bh[0]);
    head = (struct ext2_cluster_head *) (bh[0]->b_data);

    /*
     *  Do some consistency checks.
     */

    if (head->magic != cpu_to_le16(EXT2_COMPRESS_MAGIC_04X)) {
	ext2_msg(inode->i_sb,
		     "bad magic number",
		     "inode = %lu, magic = %#04x",
		     inode->i_ino, le16_to_cpu(head->magic));
	goto error;
    }
#if EXT2_GRAIN_SIZE & (EXT2_GRAIN_SIZE - 1)
# error "This code assumes EXT2_GRAIN_SIZE to be a power of two."
#endif
    /* The macro also assumes that _a > 0, _b > 0. */
#define ROUNDUP_GE(_a, _b, _d) (   (  ((_a) - 1) \
				    | ((_d) - 1)) \
				>= (  ((_b) - 1) \
				    | ((_d) - 1)))

    //mw: following 3 just for debugging!!!
    assert(!((le32_to_cpu(head->ulen) > EXT2_MAX_CLUSTER_BYTES)));
    assert(!((head->clen == 0)));
    assert(!(ROUNDUP_GE(le32_to_cpu(head->clen)
	+ head->holemap_nbytes + sizeof(struct ext2_cluster_head), 
	le32_to_cpu(head->ulen), EXT2_GRAIN_SIZE)));

    if ((le32_to_cpu(head->ulen) > EXT2_MAX_CLUSTER_BYTES)
	|| (head->clen == 0)
	|| ROUNDUP_GE(le32_to_cpu(head->clen)
		      + head->holemap_nbytes
		      + sizeof(struct ext2_cluster_head),
		      le32_to_cpu(head->ulen), EXT2_GRAIN_SIZE)) {
	ext2_msg(inode->i_sb,
		     "invalid cluster len",
		     "inode = %lu, len = %u:%u",
		     inode->i_ino,
		     le32_to_cpu(head->clen), le32_to_cpu(head->ulen));
	goto error;
    }
#undef ROUNDUP_GE

    /* TODO: Test for `nblk != 1 + ...' instead of the current
       one-sided test.  However, first look at callers, and make
       sure that they handle the situation properly (e.g. freeing
       unneeded blocks) and tht they always pass a correct
       value for nblk. */
    if (nblk <= ((le32_to_cpu(head->clen)
		  + head->holemap_nbytes + sizeof(struct ext2_cluster_head)
		  - 1)
		 / bh[0]->b_size)) {
	int i;
	ext2_msg(inode->i_sb,
		     "missing blocks",
		     "inode = %lu, blocks = %d/%u",
		     inode->i_ino, nblk, ((le32_to_cpu(head->clen)
					   + head->holemap_nbytes
					   + sizeof(struct ext2_cluster_head)
					   - 1)
					  / bh[0]->b_size) + 1);
	printk("i_size=%d\n", (int) inode->i_size);
	for (i = 0; i < 12; i++)
	    printk("i_data[%d]=%d\n", i, EXT2_I(inode)->i_data[i]);
	    printk("cluster_head (sizeof head=%u):\n\tmagic=0x%4x\n\tmethod=%d\n\t  \
	     holemap_nbytes=%d\n\tulen=%d\n\tclen=%d\n\tbh->b_size=%zu\n",
	     sizeof(struct ext2_cluster_head), head->magic,
	     (int) head->method, (int) head->holemap_nbytes, head->ulen,
	     head->clen, bh[0]->b_size);
	goto error;
    }

    /* I moved it here in case we need to load a module that
     * needs more heap that is currently allocated.
     * In such case "init_module" for that algorithm forces
     * re-allocation of ext2_wa. It should be safe here b/c the
     * first reference to ext2_wa comes just after and we have
     * locked ext2_wa before.
     *
     * FIXME: Totally separate working areas for reading and writing.
     *      Jan R.
     */
    meth = head->method;	/* only a byte, so no swabbing needed. */
    if (meth >= EXT2_N_METHODS) {
	ext2_msg(inode->i_sb,
		     "Ass: illegal method id",
		     "inode = %lu, id = %u", inode->i_ino, meth);
	dump_stack();
	goto error;
    }
    alg = ext2_method_table[meth].alg;

    /*
     *  Adjust the length if too many bytes are requested.
     *
     *    TODO: Traiter les bitmaps ici, et non plus au niveau de    
     *          l'appelant. Faire un petit cache en memorisant le    
     *          numero du dernier noeud decompresse et du dernier    
     *          cluster. Le pb, c'est qu'on ne peut pas savoir si    
     *          les blocs ont ete liberes et realloue entre temps    
     *          -> il faut etre prevenu pour invalider le buffer.    
     *
     *          pjm fixme tr: Take care of the bitmaps here,
     *          instead of by the caller as we currently do.  Keep
     *          a small cache that holds the number of the
     *          previous <inode, cluster> to have been
     *          decompressed.  The problem is that we have no way
     *          of knowing whether the blocks have been freed and
     *          reallocated in the meantime / since last time ->
     *          we must be informed so that we can invalidate the
     *          buffer.  */
    if (ulen > le32_to_cpu(head->ulen)) {
	memset(__get_cpu_var(ext2_rd_wa)->u + le32_to_cpu(head->ulen), 0, ulen - le32_to_cpu(head->ulen));
	ulen = le32_to_cpu(head->ulen);

	assert((bh[0]->b_size & (bh[nblk - 1]->b_size - 1)) == 0);
	if (((le32_to_cpu(head->clen)
	      + head->holemap_nbytes + sizeof(struct ext2_cluster_head)
	      - 1)
	     | (bh[0]->b_size - 1))
	    >= ((ulen - 1) | (bh[0]->b_size - 1))) {
	    printk(KERN_WARNING
		   "ext2_decompress_blocks: "
		   "ulen (=%zu) or clen (=%u) wrong "
		   "in dev %s, inode %lu.\n",
		   ulen, le32_to_cpu(head->clen),
		   bdevname(inode->i_sb->s_bdev, bdn), inode->i_ino);
	    goto error;
	}
    }

    /*
     *  Now, decompress data.
     */
    /* TODO: Is this (ulen == 0) possible? */
    if (ulen == 0)
	return 0;

    for (x = 0, dst = __get_cpu_var(ext2_rd_wa)->c; x < nblk; dst += bh[x++]->b_size) {
	restore_b_data_himem(bh[x]);
	memcpy(dst, bh[x]->b_data, bh[x]->b_size);
    }


    if (!ext2_algorithm_table[alg].avail) {
	ext2_msg(inode->i_sb,
		     "ext2_decompress_blocks",
		     "algorithm `%s' not available for inode %lu",
		     ext2_algorithm_table[alg].name, inode->i_ino);
	ext2_mark_algorithm_use(inode, alg);
	goto error;
    }


#ifdef EXT2_COMPR_DEBUG
    {
    	struct ext2_cluster_head *wa1head = (struct ext2_cluster_head *) __get_cpu_var(ext2_rd_wa)->c;
    	unsigned clen = le32_to_cpu(wa1head->clen);
    	if (wa1head->checksum !=
    		cpu_to_le32(ext2_adler32
			(le32_to_cpu(*(u32 *) __get_cpu_var(ext2_rd_wa)->c),
			 __get_cpu_var(ext2_rd_wa)->c + 8,
			 (sizeof(struct ext2_cluster_head) - 8 +
			  head->holemap_nbytes + clen))))
    	{
    			head->checksum = cpu_to_le32(0);
    			ext2_msg(inode->i_sb, "ext2_decompress_blocks: corrupted compressed data ",
				 "in inode %lu", inode->i_ino);
    			//goto error; 
			//mw: we try to go on. if data is corrupt we will get an compression error anyway.
    	}
    }
#endif

    count = ext2_algorithm_table[alg].decompress(__get_cpu_var(ext2_rd_wa)->c +
					     sizeof(struct
						    ext2_cluster_head) +
					     head->holemap_nbytes,
					     __get_cpu_var(ext2_rd_wa)->u,
					     __get_cpu_var(ext2_rd_wa)->heap,
					     le32_to_cpu(head->clen), ulen,
					     ext2_method_table[meth].xarg);

    /* If we got fewer than ulen bytes, there is a problem, since
       we corrected the ulen value before decompressing.  Note
       that it's OK for count to exceed ulen, because ulen can be
       less than head->ulen. */
    if ((count < ulen) || (count != le32_to_cpu(head->ulen))) {
	ext2_msg(inode->i_sb, 
		"ext2_decompress_blocks: corrupted compressed data ", "inode = %lu, count = %u of %zu (%u/%u)",
		     inode->i_ino, count, ulen, le32_to_cpu(head->clen), le32_to_cpu(head->ulen));
	goto error;
    }
    ext2_ensure_algorithm_use(inode, alg);
    return count;

  error:

    /* Raise the ECOMPR flag for this file.  What this means is
       that the file cannot be written to, and can only be read if
       the user raises the NOCOMPR flag.

       pjm 1997-01-16: I've changed it so that files with ECOMPR
       still have read permission, so user can still read the rest
       of the file but get an I/O error (errno = EXT2_ECOMPR) when
       they try to access anything from this cluster. */

    EXT2_I(inode)->i_flags |= EXT2_ECOMPR_FL;

    inode->i_ctime = CURRENT_TIME;
    mark_inode_dirty_sync(inode);
    /* pjm 1998-02-21: We used to do `memset(ext2_rd_wa.u, 0, ulen)'
       here because once upon a time the user could sometimes see
       buf contents.  I believe that this can never happen any
       more. */
    return -EXT2_ECOMPR;
}


/* ext2_calc_free_ix: Calculates the position of the C_NBLK'th non-hole
   block; equals C_NBLK plus the number of holes in the first CALC_FREE_IX()
   block positions of the cluster.

   pre: 1 =< c_nblk < EXT2_MAX_CLUSTER_BLOCKS,
        Number of 1 bits in ,ubitmap` > ,c_nblk`.
   post: c_nblk =< calc_free_ix() < EXT2_MAX_CLUSTER_BLOCKS

   Called by:
       ext2_decompress_cluster()
       ext2_file_write()

   TODO: Have ext2_compress_cluster() call this.
   */
unsigned ext2_calc_free_ix(unsigned holemap_nbytes, u8 const *holemap,
			   unsigned c_nblk)
{
    unsigned i;

    assert(1 <= c_nblk);
    assert(c_nblk < EXT2_MAX_CLUSTER_BLOCKS);
    for (i = 0; (i < holemap_nbytes * 8) && (c_nblk > 0);) {
	assert(i < EXT2_MAX_CLUSTER_BLOCKS - 1);
	if ((holemap[i >> 3] & (1 << (i & 7))) == 0)
	    c_nblk--;
	i++;
    }
    i += c_nblk;
    assert(i < EXT2_MAX_CLUSTER_BLOCKS);
    return i;
}


/*  (): Prepare the blkaddr[] array for
   decompression by moving non-hole blocks to their proper positions
   (according to ubitmap) and zeroing any other blocks.

   Returns 0 on success, -errno on error.

   Note: We assume tht blkaddr[i] won't change under us forall
   clu_block0 =< i < clu_block0 + clu_nblocks.  Holding i_sem should
   guarantee this.

   Called by:
       ext2_decompress_cluster()
       ext2_file_write() */
int
ext2_unpack_blkaddrs(struct inode *inode,
		     struct buffer_head *bh[],
		     int mmcp,
		     unsigned holemap_nbytes,
		     u8 const *holemap,
		     unsigned c_nblk,
		     unsigned free_ix,
		     unsigned clu_block0, unsigned clu_nblocks)
{
    struct ext2_bkey key;
    u32 *blkaddr;
    unsigned si, di;

    assert(clu_nblocks <= EXT2_MAX_CLUSTER_BLOCKS);
    assert(1 <= c_nblk);
    assert(c_nblk <= free_ix);
    assert(free_ix < EXT2_MAX_CLUSTER_BLOCKS);
    if (!ext2_get_key(&key, inode, clu_block0))
	return -EIO;

    if (key.ptr[key.level] == NULL) {
	/* TODO: Call ext2_error(). */
	ext2_free_key(&key);
	return -EIO;
    }

    /* impl: Note tht we're relying on clusters not straddling
       address block boundaries. */
    blkaddr = &key.ptr[key.level][key.off[key.level]];
    memset(blkaddr + free_ix,
	   0, sizeof(*blkaddr) * (clu_nblocks - free_ix));
    si = c_nblk;
    for (di = free_ix; di > si;) {
	--di;
	if (((di >> 3) < holemap_nbytes)
	    && (holemap[di >> 3] & (1 << (di & 7)))) {
	    blkaddr[di] = 0;
	    bh[di]->b_blocknr = 0;
	    clear_bit(BH_Mapped, &bh[di]->b_state);
	} else {
	    if (si == 0) {
		break;
	    }
	    blkaddr[di] = blkaddr[--si];
	    assert(bh[di]->b_blocknr == 0);
	    assert(bh[si]->b_blocknr != 0);
	    assert(buffer_mapped(bh[si]));
#ifdef  EXT2_COMPR_REPORT_CPR
	    trace_e2c("unpack: di=%d sts=0x%x si=%d blk=%ld sts=0x%x\n",
		      di, (int) bh[di]->b_state, si, bh[si]->b_blocknr,
		      (int) bh[si]->b_state);
#endif
	    bh[di]->b_blocknr = bh[si]->b_blocknr;
	    set_bit(BH_Mapped, &bh[di]->b_state);
	    bh[si]->b_blocknr = 0;
	    clear_bit(BH_Mapped, &bh[si]->b_state);
	    set_bit(BH_Uptodate, &bh[di]->b_state);
	    if (mmcp) {
		restore_b_data_himem(bh[si]);
		restore_b_data_himem(bh[di]);
		memcpy(bh[di]->b_data, bh[si]->b_data,
		       inode->i_sb->s_blocksize);
	    }
	}
    }
    if (key.level > 0)
	mark_buffer_dirty(key.ibh[key.level]);
    return ext2_free_key(&key);
}


/*
 *    Decompress one cluster.  If already compressed, the cluster
 *      is decompressed in place, and the compress bitmap is updated.
 *
 *      Returns the size of decompressed data on success, a negative
 *      value in case of failure, or 0 if the cluster was not compressed.
 *
 *      The inode is supposed to be writable.
 *
 *      Called by :
 *
 *        ext2_decompress_inode()      [sem down]
 *        ext2_file_write()            [sem down]
 *        trunc_bitmap()               [sem down]
 */
int ext2_decompress_cluster(struct inode *inode, u32 cluster)
{
    struct buffer_head *bh[EXT2_MAX_CLUSTER_BLOCKS];
    struct buffer_head *bhc[EXT2_MAX_CLUSTER_BLOCKS];
    struct page *pg[EXT2_MAX_CLUSTER_PAGES], *epg[EXT2_MAX_CLUSTER_PAGES];
    int result, nbh;
    unsigned npg, c_nblk;
    struct ext2_cluster_head *head;
    int i = 0;
    unsigned free_ix, clu_block0, clu_nblocks;
    int d_npg = -1;		/* number of decompressed page  */
    unsigned long allpagesuptodate = 1;
    struct buffer_head *bh_writeout[EXT2_MAX_CLUSTER_BLOCKS];
    int bhn_writeout;
#ifdef CONFIG_HIGHMEM
    int kmapped = 0;
#endif

    for (i = 0; i < EXT2_MAX_CLUSTER_BLOCKS; i++) {
	bh_writeout[i] = NULL;
	bhn_writeout = 0;
    }

    assert(atomic_read(&inode->i_mutex.count) <= 0);	/* i.e. mutex_lock */

    for (i = 0; i < EXT2_MAX_CLUSTER_PAGES; i++)
	epg[i] = NULL;

    /*
       Get blocks from cluster.
       Assign to variables head, ubitmap, clu_block0, clu_nblocks.
       Shuffle blkaddr[] array and write zero to holes.
       Allocate new blocks.
       Get the working area.
       Decompress.
       Copy to bh[]->b_data (marking buffers uptodate and dirty).
       Release working area.
       Release bh[].
     */

    nbh = 0;
    npg = ext2_cluster_npages(inode, cluster);
    result = ext2_get_cluster_pages(inode, cluster, pg, NULL, 0);
    if (result <= 0) {
	for (i = 0; i < npg; i++)
	    epg[i] = NULL;
	goto out_err;
    }

    for (i = 0; i < npg; i++) {
	if ((pg[i]->index <= ((inode->i_size - 1) >> PAGE_CACHE_SHIFT)) &&
	    !PageUptodate(pg[i])) {
	    allpagesuptodate = 0;
	}
    }
    if (allpagesuptodate) {
	//printk("DecompressPages: Ino:%lu\n", inode->i_ino);
	result = ext2_decompress_pages(inode, cluster, pg);
	if (result != 0) {
	    for (i = 0; i < npg; i++)
		epg[i] = NULL;
	    if (result > 0)
		goto cleanup;
	    else
		goto out_err;
	}
	/*mw: if we continue here then in ext2_decompress_pages
	 * not all pages were up-to-date 
	 */
    }
    //printk("DecompressCluster: Ino:%lu\n", inode->i_ino);
    result = ext2_get_cluster_extra_pages(inode, cluster, pg, epg);
    if (result <= 0) {
	goto out_err;
    }
#ifdef CONFIG_HIGHMEM
    ext2_kmap_cluster_pages(NULL, pg, epg);
    kmapped = 1;
#endif

    result = ext2_get_cluster_blocks(inode, cluster, bh, pg, epg, 0);
    if (result <= 0) {
	goto out_err;
    }
    nbh = c_nblk = result;


#ifdef EXT2_COMPR_REPORT
    {
	int j;
	printk
	    (" > > > ext2_decompress_cluster %d: inode=%ld, size=%d nbh=%d\n",
	     cluster, inode->i_ino, (int) inode->i_size, nbh);
#ifdef EXT2_COMPR_REPORT_VERBOSE
	for (j = 0; j < nbh; j++) {
	    if (bh[j]) {
		printk("0buffer_head[%d]: blocknr=%lu,  addr=%p \n", j,
		       (unsigned long) bh[j]->b_blocknr, bh[j]);
		if (bh[j]->b_page)
		    printk("0:[page->index=%ld]\n", bh[j]->b_page->index);
		else
		    printk("[No page]\n");
	    } else
		printk("buffer_head[%d] is NULL\n", j);
	}
	while ((j < EXT2_MAX_CLUSTER_BLOCKS) && (bh[j] != NULL) && bh[j]->b_blocknr) {	/*Add by Yabo Ding */
	    printk
		("buffer_head[%d] is free but not NULL: blocknr=%lu, addr=%p\n",
		 j, (unsigned long) bh[j]->b_blocknr, bh[j]);
	    j++;
	}
#endif
    }
#endif
    for (i = 0; i < nbh; i++)
	assert(bh[i]->b_blocknr != 0);

    restore_b_data_himem(bh[0]);

    head = (struct ext2_cluster_head *) bh[0]->b_data;
    if (head->magic != cpu_to_le16(EXT2_COMPRESS_MAGIC_04X)) {
	ext2_msg(inode->i_sb,
		     "ext2_decompress_cluster: bad magic number",
		     "cluster %d: inode = %lu, magic = %#04x",
		     cluster, inode->i_ino, le16_to_cpu(head->magic));
	EXT2_I(inode)->i_flags |= EXT2_ECOMPR_FL;
	result = -EXT2_ECOMPR;
	goto out_err;
    }
    if (le32_to_cpu(head->ulen) -
	(c_nblk << inode->i_sb->s_blocksize_bits) <= 0) {
	ext2_error(inode->i_sb, "ext2_decompress_cluster",
		   "ulen too small for c_nblk.  ulen=%u, c_nblk=%u, bs=%lu",
		   le32_to_cpu(head->ulen), c_nblk,
		   inode->i_sb->s_blocksize);
	EXT2_I(inode)->i_flags |= EXT2_ECOMPR_FL;
	result = -EXT2_ECOMPR;
	goto out_err;
    }
    free_ix =
	ext2_calc_free_ix(head->holemap_nbytes, (u8 const *) (&head[1]),
			  c_nblk);
    clu_block0 = ext2_cluster_block0(inode, cluster);
    clu_nblocks = ext2_cluster_nblocks(inode, cluster);
    ext2_unpack_blkaddrs(inode, bh, 1,
			 head->holemap_nbytes, (u8 const *) (&head[1]),
			 c_nblk, free_ix, clu_block0, clu_nblocks);

    /* Allocate the extra blocks needed. */
    {
	int data_left = le32_to_cpu(head->ulen);

	data_left -= c_nblk << inode->i_sb->s_blocksize_bits;
	assert(data_left > 0);
	for (i = free_ix; i < clu_nblocks; i++)
	    if (((i >> 3) >= head->holemap_nbytes)
		|| !(head->holemap[i >> 3] & (1 << (i & 7)))) {
		result = ext2_get_block(inode,
					clu_block0 + i,
					bh[i], 1 /* create */ );
		if (bh[i]->b_blocknr == 0)
		    goto out_err;
		d_npg =
		    (i >>
		     (PAGE_CACHE_SHIFT - inode->i_sb->s_blocksize_bits)) +
		    1;
		nbh++;
		data_left -= inode->i_sb->s_blocksize;
		if (data_left <= 0)
		    break;
	    }
    }

    /* jmr 1998-10-28 Hope this is the last time I'm moving this code.
     * Module loading must be done _before_ we lock wa, just think what
     * can happen if we reallocate wa when somebody else uses it...
     */
    {
	unsigned meth;
#ifdef CONFIG_KMOD
	unsigned alg;
#endif

	meth = head->method;	/* only a byte, so no swabbing needed. */
	if (meth >= EXT2_N_METHODS) {
	    ext2_msg(inode->i_sb,
			 "Ass.: illegal method id",
			 "inode = %lu, id = %u", inode->i_ino, meth);
	    result = -EXT2_ECOMPR;
	    goto out_err;
	}
#ifdef CONFIG_KMOD
	alg = ext2_method_table[meth].alg;
	if (!ext2_algorithm_table[alg].avail) {
	    char str[32];

	    sprintf(str, "ext2-compr-%s", ext2_algorithm_table[alg].name);
	    request_module(str);
	}
#endif
    }

    result = -EINTR;

    /*
     *  Then, decompress and copy back data.
     */
    {
	int ic;

	for (ic = 0, i = 0; i < clu_nblocks; i++) {
	    if (bh[i]->b_blocknr != 0) {
		bhc[ic] = bh[i];
		ic++;
		if (ic == c_nblk) {
		    break;
		}
	    }
	}
    }


#ifdef EXT2_COMPR_REPORT_WA
    printk(KERN_DEBUG "pid %d locks wa\n", current->pid);
#endif
    if (get_cpu_var(ext2_rd_wa) == NULL)
    {
	 ext2_alloc_rd_wa();
    }
    assert(__get_cpu_var(ext2_rd_wa) != NULL);

    result = ext2_decompress_blocks(inode, bhc, c_nblk,
				    le32_to_cpu(head->ulen), cluster);
    if (result != (int) le32_to_cpu(head->ulen)) {
	if (result >= 0) {
	    /* I think this is impossible, as
	       ext2_decompress_blocks() checks against
	       head->ulen. */
	    printk(KERN_WARNING "Unexpected return value %d "
		   "from ext2_decompress_blocks()\n", result);
	    result = -EXT2_ECOMPR;
	}

#ifdef EXT2_COMPR_REPORT_WA
	printk(KERN_DEBUG "pid %d unlocks wa\n", current->pid);
#endif
	put_cpu_var(ext2_rd_wa);
	goto out_err;
    }

#ifdef EXT2_COMPR_REPORT
    printk(KERN_DEBUG "ext2: %04x:%lu: cluster %d+%d [%d] "
	   "decompressed into %d bytes\n",
	   inode->i_rdev,
	   inode->i_ino, clu_block0, clu_nblocks, c_nblk, result);
#endif

    /* Copy back decompressed data. */
    {
	int count = result;
	unsigned char const *src;
	int c, p;
	int cbh;
	int n;			/* block index in page  */
	struct buffer_head *bp;
	unsigned addr0, b_start, b_end;

	assert(count > 0);
	if (d_npg == -1) {
	    d_npg = ((count - 1) >> PAGE_CACHE_SHIFT) + 1;
	}
#ifdef  EXT2_COMPR_REPORT_CPR
	trace_e2c
	    ("ext2_decompress_cluster: cnt=%d free_ix=%d d_npg=%d nbh=%d\n",
	     count, free_ix, d_npg, nbh);
#endif
	result = -EXT2_ECOMPR;
	src =  __get_cpu_var(ext2_rd_wa)->u;
	cbh = 0;
	for (c = 0; c < clu_nblocks; c++) {

	    if (bh[c]->b_blocknr == 0) {
#ifdef  EXT2_COMPR_REPORT_CPR
		trace_e2c("\t clear buf %d sts=0x%x\n", c,
			  (int) bh[c]->b_state);
#endif
		restore_b_data_himem(bh[c]);
		memset(bh[c]->b_data, 0, inode->i_sb->s_blocksize);
		continue;
	    }
	    if (cbh >= (nbh - 1)) {
		break;
	    }
	    if (count < inode->i_sb->s_blocksize) {
		put_cpu_var(ext2_rd_wa);
		goto out_err;
	    }
	    cbh++;
	    count -= inode->i_sb->s_blocksize;
	    p = c >> (PAGE_CACHE_SHIFT - inode->i_sb->s_blocksize_bits);
	    if (!PageUptodate(pg[p])) {
		addr0 = (clu_block0 << inode->i_sb->s_blocksize_bits);
		b_start = addr0 + (c << inode->i_sb->s_blocksize_bits);
		b_end = b_start + inode->i_sb->s_blocksize;
#ifdef  EXT2_COMPR_REPORT_CPR
		trace_e2c("\t[%d] sts=0x%x e=%d s=%d sz=%d pg:%lu(%#x)\n",
			  c, (int) bh[c]->b_state, b_end, b_start,
			  (int) inode->i_size, pg[p]->index,
			  (unsigned int) pg[p]);
#endif
		if (b_end <= inode->i_size) {
		    /* Block is before end of file, copy data */
		    restore_b_data_himem(bh[c]);
		    memcpy(bh[c]->b_data, src, inode->i_sb->s_blocksize);

		} else if (b_start < inode->i_size) {
		    /* Block contains end of file, copy to end */
		    restore_b_data_himem(bh[c]);
		    memcpy(bh[c]->b_data, src, inode->i_size - b_start);

		}
		set_buffer_uptodate(bh[c]);
		set_buffer_dirty(bh[c]);
		bh_writeout[bhn_writeout] = bh[c];	//mw
		bhn_writeout++;	//mw
	    } else {
		//mw: DEBUG. buffer is uptodate now. compress will not reread! an get the compressed data!!!
		// clear flag in extra page!!!
		// clear_bit(BH_Uptodate, &bh[c]->b_state);

		n = c & ((PAGE_CACHE_SIZE - 1) >> inode->i_sb->
			 s_blocksize_bits);
		bp = page_buffers(pg[p]);
		for (i = 0; i < n; i++) {
		    bp = bp->b_this_page;
		}
		result = ext2_get_block(inode, clu_block0 + c, bp, 0);

		//mw: needed to do a writeback of the non-epg-buffers
		//no idea how it was done before
		set_buffer_uptodate(bp);
		set_buffer_dirty(bp);
		bh_writeout[bhn_writeout] = bp;	//mw
		bhn_writeout++;	//mw

		if (bp->b_blocknr == 0) {
			put_cpu_var(ext2_rd_wa);
			goto out_err;
		}
		assert(bp->b_blocknr == bh[c]->b_blocknr);
	    }
	    src += inode->i_sb->s_blocksize;
	}
	if (count > inode->i_sb->s_blocksize) {
		put_cpu_var(ext2_rd_wa);
		goto out_err;
	}
	p = c >> (PAGE_CACHE_SHIFT - inode->i_sb->s_blocksize_bits);
	if (!PageUptodate(pg[p])) {
	    addr0 = (clu_block0 << inode->i_sb->s_blocksize_bits);
	    b_start = addr0 + (c << inode->i_sb->s_blocksize_bits);
#ifdef  EXT2_COMPR_REPORT_CPR
	    trace_e2c("\t[%d] sts=0x%x c=%d s=%d sz=%d pg:%lu(%#x)\n", c,
		      (int) bh[c]->b_state, count, b_start,
		      (int) inode->i_size, pg[p]->index,
		      (unsigned int) pg[p]);
#endif
	    if (b_start >= inode->i_size) {
		restore_b_data_himem(bh[c]);
		memset(bh[c]->b_data, 0, inode->i_sb->s_blocksize);

	    } else {
		if ((inode->i_size - b_start) < count) {
		    restore_b_data_himem(bh[c]);
		    memcpy(bh[c]->b_data, src, inode->i_size - b_start);
		    memset(bh[c]->b_data + (inode->i_size - b_start), 0,
			   count - (inode->i_size - b_start));
		} else {
		    restore_b_data_himem(bh[c]);
		    memcpy(bh[c]->b_data, src, count);
		}
	    }
	    set_buffer_uptodate(bh[c]);
	    set_buffer_dirty(bh[c]);
	    bh_writeout[bhn_writeout] = bh[c];	//mw
	    bhn_writeout++;	//mw
	} else {
	    assert(epg[p] != NULL);	//mw
	    n = c & ((PAGE_CACHE_SIZE - 1) >> inode->i_sb->
		     s_blocksize_bits);
	    bp = page_buffers(pg[p]);
	    for (i = 0; i < n; i++) {
		bp = bp->b_this_page;
	    }
	    result = ext2_get_block(inode, clu_block0 + c, bp, 0);

	    //mw: needed to do a writeback of the non-epg-buffers
	    //no idea how it was done before
	    set_buffer_uptodate(bp);
	    set_buffer_dirty(bp);
	    bh_writeout[bhn_writeout] = bp;	//mw
	    bhn_writeout++;	//mw
	    if (bp->b_blocknr == 0) {
		put_cpu_var(ext2_rd_wa);
		goto out_err;
	    }
	    assert(bp->b_blocknr == bh[c]->b_blocknr);
	}
	result = (nbh - 1) * inode->i_sb->s_blocksize + count;
    }

    for (i = 0; i < EXT2_MAX_CLUSTER_PAGES; i++) {
	if (pg[i] == NULL)
	    break;
	if (i < d_npg)
	    SetPageUptodate(pg[i]);
    }

#ifdef EXT2_COMPR_REPORT_WA
    printk(KERN_DEBUG "pid %d unlocks wa\n", current->pid);
#endif
    put_cpu_var(ext2_rd_wa);

    inode->i_ctime = CURRENT_TIME;
    mark_inode_dirty_sync(inode);
    /* If needed, EXT2_DIRTY_FL is raised by the caller. */

#if 0
    /* TODO: SYNC */
    if (IS_SYNC(inode)) {
	generic_osync_inode(inode, inode->i_mapping,
			    OSYNC_METADATA | OSYNC_DATA);
    }
#endif
    assert(result >= 0);

    //Sync out changes:
    assert(bhn_writeout <= EXT2_MAX_CLUSTER_BLOCKS);
    assert(bhn_writeout >= 0);

    //mw: debug
    for (i = 0; i < bhn_writeout; i++) {
	if ((!buffer_mapped(bh_writeout[i]))
	    || (bh_writeout[i]->b_bdev == NULL)) {
	    u32 block = ext2_cluster_block0(inode, cluster);
	    ext2_get_block(inode, block + i, bh_writeout[i], 1);
	    //printk("ext2_get_block Block:%lu, Mapped:%i, Page:%lu, bdev: %#x\n", bh_writeout[i]->b_blocknr, (bh_writeout[i]->b_state & BH_Mapped), (bh_writeout[i]->b_page ? bh_writeout[i]->b_page->index : 0), bh_writeout[i]->b_bdev );
	}
	assert(buffer_mapped(bh_writeout[i]));
	assert(bh_writeout[i]->b_bdev != NULL);
	assert(bh_writeout[i]->b_bdev == inode->i_sb->s_bdev);
	/*if (bh_writeout[i]->b_bdev == NULL)
	   bh_writeout[i]->b_bdev = inode->i_sb->s_bdev; //fix bdev-bug */
    }

    ll_rw_block(WRITE, bhn_writeout, bh_writeout);
    //mw: seems we have to wait here, otherwise: crash!

    CHECK_NOT_ATOMIC
    for (i = 0; i < bhn_writeout; i++) {
	if (bh_writeout[i])
	    wait_on_buffer(bh_writeout[i]);
    }
    goto cleanup;

  out_err:
    printk("Error in Decompressing cluster: Err=%i\n", result);

  cleanup:

#ifdef CONFIG_HIGHMEM
    if (kmapped)
	ext2_kunmap_cluster_pages(NULL, pg, epg);
#endif

    for (i = 0; i < EXT2_MAX_CLUSTER_PAGES; i++) {
	if (pg[i] == NULL)
	    break;
	unlock_page(pg[i]);
	page_cache_release(pg[i]);
    }

    for (i = 0; i < EXT2_MAX_CLUSTER_PAGES; i++) {
	if (epg[i] != NULL) {
	    ClearPageDirty(epg[i]);
	    ClearPageUptodate(epg[i]);
	    try_to_free_buffers(epg[i]);
	    unlock_page(epg[i]);
	    assert(page_count(epg[i]) == 1);
	    page_cache_release(epg[i]);
	}
    }

    /*
     *  Release buffers, don't forget to unlock the locked ones.
     *  pjm 1998-01-14: TO_DO: Locked ones?
     */
    assert(nbh >= 0);
    assert(nbh <= EXT2_MAX_CLUSTER_BLOCKS);
#ifdef EXT2_COMPR_REPORT
    trace_e2c(" < < < ext2_decompress_cluster %d: inode=%ld, res=%i\n",
	      cluster, inode->i_ino, result);
#endif
    return result;
}


/*
 * Function to decompress the pages of a cluster.
 *
 *	Allocate buffers to pages what are not mapped on the device.
 *
 *      Returns the size of decompressed data on success, a negative
 *      value in case of failure, or 0 if some pages are not uptodate.
 *
 *      The inode is supposed to be writable.
 *	All the pages must be UPTODATE, 
 */
int ext2_decompress_pages(struct inode *inode, u32 cluster,
			  struct page *pg[])
{
    struct ext2_cluster_head *head;
    struct buffer_head *bh0;
    struct buffer_head *bh[EXT2_MAX_CLUSTER_BLOCKS];
    unsigned nbh, c_nblk;
    unsigned free_ix, clu_block0, clu_nblocks;
    int i, pagesPerCluster, data_left, size = 0;
    long status = 0;
    char *dp;
    struct buffer_head *bh_writeout[EXT2_MAX_CLUSTER_BLOCKS];
    int bhn_writeout;
#ifdef CONFIG_HIGHMEM
    int kmapped = 0;

    ext2_kmap_cluster_pages(NULL, pg, NULL);
    kmapped = 1;
#endif

    for (i = 0; i < EXT2_MAX_CLUSTER_BLOCKS; i++) {
	bh_writeout[i] = NULL;
	bhn_writeout = 0;
    }

    /* First, get cluster_head (For this, we need to re-read the first block of
       the cluster, without overwriting the data of the page the buffer point to... */
    /* This suppose that cluster are aligned with PAGE_SIZE... To be improved */

    /* Changed by Yabo Ding<bobfree_cn@yahoo.com.cn>,<yding@wyse.com>
       The old code cannot reread data from disk to a changed buffers data pointer in 2.6.x.
       So, I copy memory data(decompressed) to a temporary buffer;
       Then reread data(compressed) from disk, and copy to head;
       Then copy back the memory data from temporary buffer.
       It seems clumsy, but it works well.
     */

    bh0 = page_buffers(pg[0]);
    restore_b_data_himem(bh0);

    head = (struct ext2_cluster_head *) kmalloc(bh0->b_size, GFP_KERNEL);
    if (head == NULL) {
	ext2_msg(inode->i_sb, "no more memory", "inode = %lu",
		     inode->i_ino);
	status = -EIO;
	goto out_x;
    }
    dp = kmalloc(bh0->b_size, GFP_KERNEL);
    if (dp == NULL) {
	ext2_msg(inode->i_sb, "no more memory", "inode = %lu",
		     inode->i_ino);
	kfree(head);
	status = -EIO;
	goto out_x;
    }
    memcpy(dp, bh0->b_data, bh0->b_size);
    clear_bit(BH_Uptodate, &bh0->b_state);
    if (!buffer_mapped(bh0)) {
	status =
	    ext2_get_block(inode, ext2_cluster_block0(inode, cluster), bh0,
			   0);
	if (bh0->b_blocknr == 0) {
	    trace_e2c
		("ext2_decompress_pages: ext2_get_block error %ld (cluster = %u)\n",
		 status, cluster);
	    kfree(head);
	    memcpy(bh0->b_data, dp, bh0->b_size);
	    kfree(dp);
	    status = -EIO;
	    goto out;
	}
    }
    ll_rw_block(READ, 1, &bh0);

    CHECK_NOT_ATOMIC
    wait_on_buffer(bh0);
    //printk("RE-Read: Buffer: blocknr:%lu(%#x) \n", bh0->b_blocknr, bh0);
    if (!buffer_uptodate(bh0)) {	/* Read error ??? */
	trace_e2c("ext2_decompress_pages: IO error (cluster = %u)\n",
		  cluster);
	kfree(head);
	memcpy(bh0->b_data, dp, bh0->b_size);
	kfree(dp);
	status = -EIO;
	goto out;
    }
    /* This suppose that cluster are aligned with PAGE_SIZE... To be improved 
       bh0->b_data = page_address(pg[0]);                                    */
    memcpy((char *) head, bh0->b_data, bh0->b_size);
    memcpy(bh0->b_data, dp, bh0->b_size);
    kfree(dp);

    if (head->magic != cpu_to_le16(EXT2_COMPRESS_MAGIC_04X)) {
	ext2_msg(inode->i_sb,
		     "ext2_decompress_pages: bad magic number",
		     "inode = %lu, magic = %#04x", inode->i_ino,
		     le16_to_cpu(head->magic));
	kfree(head);
	status = -EIO;
	goto out;
    }
#ifdef  EXT2_COMPR_REPORT
    trace_e2c("ext2_decompress_pages: clt=%d i=%ld head=0x%x\n", cluster,
	      inode->i_ino, (unsigned) head);
#endif

    /* Now, try to do the same as in ext2_decompress_cluster for moving/allocating blocks */
    nbh = 0;
    pagesPerCluster = ext2_cluster_npages(inode, cluster);
    for (i = 0; i < pagesPerCluster && pg[i]; i++) {
	assert(PageLocked(pg[i]));
	//if (!(PageUptodate(pg[i]))) {                 
	//mw: do it like ext2_decompress_cluster to handle end of a file correctly
	if (!(PageUptodate(pg[i]))
	    && (pg[i]->index <= ((inode->i_size - 1) >> PAGE_CACHE_SHIFT))) {
	    kfree(head);
	    printk("should never happen: not all pages uptodate!\n");	//mw
	    status = 0;
	    goto out_x;
	}
    }

    for (i = 0; i < pagesPerCluster && pg[i]; i++) {
	struct buffer_head *bhead, *bhx;
	int idx = 0;

	/* assert(PageUptodate(pg[i])); with ftruncate() can be false */
	if (!page_has_buffers(pg[i])) {
	    ClearPageUptodate(pg[i]);	/*mw */
	    ClearPageDirty(pg[i]);	/*mw */
	    assert(0);
	    create_empty_buffers_e2c(pg[i], inode->i_sb->s_blocksize, 0,
				     inode);
	    if (unlikely(!page_has_buffers(pg[i])))
		printk("Error: NOMEM!\n");
	}
	bhead = page_buffers(pg[i]);
	for (bhx = bhead; bhx != bhead || !idx; bhx = bhx->b_this_page) {
	    idx++;
	    bh[nbh] = bhx;
	    nbh++;
	}
    }

    while ((nbh != 0) && (bh[nbh - 1]->b_blocknr == 0))
	--nbh;

    c_nblk = nbh;

    free_ix =
	ext2_calc_free_ix(head->holemap_nbytes, (u8 const *) (&head[1]),
			  c_nblk);
    clu_block0 = ext2_cluster_block0(inode, cluster);
    clu_nblocks = ext2_cluster_nblocks(inode, cluster);
    ext2_unpack_blkaddrs(inode, bh, 0, head->holemap_nbytes,
			 (u8 const *) (&head[1]), c_nblk, free_ix,
			 clu_block0, clu_nblocks);

    /* Allocate the extra blocks needed. */
    data_left = size = le32_to_cpu(head->ulen);

    data_left -= c_nblk << inode->i_sb->s_blocksize_bits;
    assert(data_left > 0);
    for (i = 0; i < free_ix; i++) {
	if (bh[i]->b_blocknr != 0) {
#ifdef  EXT2_COMPR_REPORT_CPR
	    trace_e2c("\t [%d] blk=%ld sts=0x%x\n", i, bh[i]->b_blocknr,
		      (int) bh[i]->b_state);
#endif
	    set_buffer_dirty(bh[i]);
	    bh_writeout[bhn_writeout] = bh[i];	//mw
	    bhn_writeout++;	//mw
	}
    }

    for (i = free_ix; i < clu_nblocks; i++) {
	if (((i >> 3) >= head->holemap_nbytes)
	    || !(head->holemap[i >> 3] & (1 << (i & 7)))) {
	    status =
		ext2_get_block(inode, clu_block0 + i, bh[i],
			       1 /* create */ );
	    if (status || bh[i]->b_blocknr == 0) {
		status = -EIO;
		goto out;
	    }
#ifdef  EXT2_COMPR_REPORT_CPR
	    trace_e2c("\t [%d] blk=%ld sts=0x%x\n", i, bh[i]->b_blocknr,
		      (int) bh[i]->b_state);
#endif
	    set_bit(BH_Uptodate, &bh[i]->b_state);
	    set_buffer_dirty(bh[i]);
	    bh_writeout[bhn_writeout] = bh[i];	//mw
	    bhn_writeout++;	//mw
	    nbh++;
	    data_left -= inode->i_sb->s_blocksize;
	    if (data_left <= 0)
		break;
	}
    }

  out:
    kfree(head);

  out_x:

    for (i = 0; i < bhn_writeout; i++) {

	if ((!buffer_mapped(bh_writeout[i]))
	    || (bh_writeout[i]->b_bdev == NULL)) {
	    u32 block = ext2_cluster_block0(inode, cluster);
	    ext2_get_block(inode, block + i, bh_writeout[i], 1);
	    //printk("ext2_get_block Block:%lu, Mapped:%i, Page:%lu, bdev: %#x\n", bh_writeout[i]->b_blocknr, (bh_writeout[i]->b_state & BH_Mapped), (bh_writeout[i]->b_page ? bh_writeout[i]->b_page->index : 0), bh_writeout[i]->b_bdev );
	}
	assert(buffer_mapped(bh_writeout[i]));
	assert(bh_writeout[i]->b_bdev != NULL);
	assert(bh_writeout[i]->b_bdev == inode->i_sb->s_bdev);
	/*if (bh_writeout[i]->b_bdev == NULL)
	   bh_writeout[i]->b_bdev = inode->i_sb->s_bdev; //fix bdev-bug */
    }
    //Sync out changes:
    ll_rw_block(WRITE, bhn_writeout, bh_writeout);
    //mw: seems we have to wait here, otherwise: crash!

    CHECK_NOT_ATOMIC
    for (i = 0; i < bhn_writeout; i++) {
	if (bh_writeout[i])
	    wait_on_buffer(bh_writeout[i]);
    }


#ifdef CONFIG_HIGHMEM
    if (kmapped)
	ext2_kunmap_cluster_pages(NULL, pg, NULL);
#endif

    return (status ? status : size);
}


/* Decompress every cluster that is still compressed.
   We stop and return -ENOSPC if we run out of space on device.

   The caller needs to check for EXT2_COMPRBLK_FL before calling.

   Returns 0 on success, -errno on failure.

   Called by ext2_ioctl(). */
int ext2_decompress_inode(struct inode *inode)
{
    u32 cluster;
    u32 n_clusters;
    int err = 0;
    struct ext2_inode_info *ei = EXT2_I(inode);

    assert(ei->i_flags & EXT2_COMPRBLK_FL);

    /* Quotas aren't otherwise kept if file is opened O_RDONLY. */
    dquot_initialize(inode);
    
    //mutex_lock(&inode->i_mutex); /* MW 5-16-07 */
    assert(atomic_read(&inode->i_mutex.count) <= 0);	/* i.e. mutex_lock */
    err = 0;
    /* This test can succeed because down() (and I think DQUOT_INIT) can block. */
    if (!(ei->i_flags & EXT2_COMPRBLK_FL))
	goto out;

    n_clusters = ext2_n_clusters(inode);
    for (cluster = 0; cluster < n_clusters; cluster++) {
	err = ext2_cluster_is_compressed_fn(inode, cluster);
	if (err > 0) {
	    err = ext2_decompress_cluster(inode, cluster);
	    /* If we later get an error, we'll need to recompress. */
	    ei->i_flags |= EXT2_DIRTY_FL;
	    ei->i_compr_flags |= EXT2_CLEANUP_FL;
	}
	if (err < 0)
	    goto error;
    }
    assert(err >= 0);
    err = 0;
    ei->i_flags &= ~(EXT2_COMPRBLK_FL | EXT2_DIRTY_FL);
    ei->i_compr_flags &= ~EXT2_CLEANUP_FL;
  error:
    inode->i_ctime = CURRENT_TIME;
    mark_inode_dirty_sync(inode);
  out:
//      mutex_unlock(&inode->i_mutex); /* MW 5-16-07 */
    return err;
}


/*
   TODO: SECRM_FL

   TODO: Avant de liberer les blocs, regarder si le compteur
   est a 1, et marquer le noeud si ce n'est pas le cas
   (pour preparer la recompression immediate).        

   pjm fixme translation.
   "Before freeing the blocks, check if the counter is 1, 
   and mark the inode if not (in order to prepare for
   immediate recompression)." */

/* This is called by ext2_compress_cluster to free the blocks now
   available due to compression.  We free ,nb` blocks beginning with
   block ,block`.  We set the address of each freed block to
   EXT2_COMPRESSED_BLKADDR, thus marking the cluster as compressed.
   N.B. It is up to the caller to adjust i_blocks. */

/* TODO: ext2_truncate() is much more careful than this routine.
   (E.g. it checks for bh->b_count > 1, and checks for things changing
   underneath it.  It also calls bforget instead of brelse if it's
   going to free it.)  Why?  Maybe we should copy it. */

/* effic: Reduce the number of calls to ext2_free_block() the way
   ext2_trunc_direct() does. */

/* fixme: I think tht we do indeed need to check if buffers are held by
   somebody else before freeing them. */
static int ext2_free_cluster_blocks(struct inode *inode, u32 block,
				    unsigned nb)
{
    u32 tmp;
    struct ext2_bkey key;
    int err;

/*
 * whitpa 04 Oct 2004: although it may be true that using e2compr in
 * conjunction with quotas is a Bad Idea, having quotas enabled for other
 * filesystems doesn't necessarily mean that the quota feature will actually be
 * used in this one, so many people find the following assertion very annoying.
 * I have therefore disabled it.
 */
/*	assert (!inode->i_sb->dq_op || (inode->i_flags & S_QUOTA)); */
    if (!nb)
	return 0;
    if (nb > EXT2_MAX_CLU_NBLOCKS) {
	assert((int) nb >= 0);
	assert(nb <= EXT2_MAX_CLU_NBLOCKS);
	return -EDOM;
    }
    assert(((block + nb) & 3) == 0);
    if (!ext2_get_key(&key, inode, block))
	return -EIO;

    while (nb-- > 0) {
	tmp = ext2_get_key_blkaddr(&key);
	err = ext2_set_key_blkaddr(&key, EXT2_COMPRESSED_BLKADDR);
	if (err)
	    goto out;
	if (tmp != 0) {
	    assert(tmp != EXT2_COMPRESSED_BLKADDR);
#ifdef EXT2_COMPR_REPORT_ALLOC
	    printk(KERN_DEBUG "ext2: free %d = (%d) %d:%d:%d:%d : %d\n",
		   key.block,
		   key.level,
		   key.off[0], key.off[1], key.off[2], key.off[3], tmp);
#endif
	    ext2_free_blocks(inode, tmp, 1);
	}
	if (!ext2_next_key(&key, 1))
	    break;
    }
    err = 0;
  out:
    ext2_free_key(&key);
    return err;
}

#ifdef EXT2_COMPR_DEBUG
static unsigned count_bits(unsigned char *p, unsigned nb)
{
    u32 x = le32_to_cpu(*(u32 *) p);
    unsigned n = 0;

    assert(nb <= 4);
    if (nb != 4)
	x &= (1 << (nb * 8)) - 1;
    while (x) {
	x &= (x - 1);
	n++;
    }
    return n;
}
#endif

/*
 * __remove_compr_assoc_queue is used in invalidate_inode_buffers
 * replacement code for ext2_compress_cluster(). TLL 02/21/07
 * Yeah, it is duplicate code, but using it does not require
 * patching fs/buffer.c/__remove_assoc_queue to export it.
 * The buffer's backing address_space's private_lock must be held.
 */
/*static inline void __remove_compr_assoc_queue(struct buffer_head *bh)
{
	list_del_init(&bh->b_assoc_buffers);
}*/

/* Compress one cluster.  If the cluster uses fewer blocks once
   compressed, it is stored in place of the original data.  Unused
   blocks are freed, and the cluster is marked as compressed.

   Returns a negative value on error,
   0 if the cluster does not compress well,
   positive if it is compressed (whether it was already compressed
   or whether we compressed it).
 
   Assume inode is writable.
 
   Called by :
 
         ext2_cleanup_compressed_inode () [i_sem] 

   If ever we acquire new callers, make sure that quotas are
   initialised, and COMPRBLK is handled correctly (i.e. such
   that ioctl() can't change the cluster size on us), and that caller
   tests for ext2_wa==NULL.
*/

int ext2_compress_cluster(struct inode *inode, u32 cluster)
{
    struct buffer_head *bh[EXT2_MAX_CLUSTER_BLOCKS + 1];
    struct page *pg[EXT2_MAX_CLUSTER_PAGES];
    int s_nblk;			/* Equals clu_nblocks less any trailing hole blocks. */
    unsigned u_nblk = (~(unsigned) 0), c_nblk;	/* Number of blocks occupied by
						   un/compressed data. */
    int result, n, x;
    int ulen, maxlen = 0, clen = 0;
    unsigned char *dst;
    u8 *src;
    unsigned meth, alg;
    int nbh = 0, npg, i;
    unsigned char holemap_nbytes = 0;
    unsigned last_hole_pos;
    struct ext2_cluster_head *head;
    unsigned r_nblk;
    struct ext2_inode_info *ei = EXT2_I(inode);
    unsigned long saved_isize;
    //int dotrunc = 1; //mw

#ifdef CONFIG_HIGHMEM
    int kmapped = 0;
#endif

    /* impl: Otherwise, ioctl() could change the cluster size
       beneath us. */
    /* TLL say not compressed and return -1 6-15-07 */
    if (!(ei->i_flags & EXT2_COMPRBLK_FL))
	return -1;

    //mw
    saved_isize = inode->i_size;

    assert(atomic_read(&inode->i_mutex.count) <= 0);	/* i.e. mutex_lock */
    assert(!mapping_mapped(inode->i_mapping));

    npg = ext2_cluster_npages(inode, cluster);

    result = ext2_get_cluster_pages(inode, cluster, pg, NULL, 1);
    if (result <= 0)
	goto done;

#ifdef CONFIG_HIGHMEM
    ext2_kmap_cluster_pages(NULL, pg, NULL);
    kmapped = 1;
#endif

    /* effic: We ought to use the page cache.  Using the page
       cache always costs extra CPU time, but saves I/O if the
       page is present.  We still need to detect holes, which
       unfortunately may still cause I/O.  Testing for all-zero
       could save us that I/O. */

    nbh = ext2_get_cluster_blocks(inode, cluster, bh, pg, NULL, 1);

    s_nblk = nbh;

#ifdef EXT2_COMPR_REPORT
    {
	int i;
	trace_e2c(" > > > ext2_compress_cluster %d: inode=%ld, size=%d\n",
		  cluster, inode->i_ino, (int) inode->i_size);
#ifdef EXT2_COMPR_REPORT_CPR
	for (i = 0; i < s_nblk; i++) {
	    if (bh[i]) {
		printk(KERN_DEBUG
		       "bbuffer_head[%d]: blocknr=%lu, addr=0x%p ", i,
		       (unsigned long) bh[i]->b_blocknr, bh[i]);
		if (bh[i]->b_page)
		    printk(KERN_DEBUG "bgn:[page->index=%ld]\n",
			   bh[i]->b_page->index);
		else
		    printk(KERN_DEBUG "[No page]\n");
	    } else
		printk("bbuffer_head[%d] is NULL\n", i);
	}
#endif
    }
#endif
    /*
     *  Did somebody else compress the cluster while we were waiting ?
     *  This should never arise ...
     */
    result = ext2_cluster_is_compressed_fn(inode, cluster);
    if (result != 0) {
	if (result > 0) {
	    ext2_msg(inode->i_sb,
			 "ext2_compress_cluster",
			 "compressing compressed cluster");
	}
	goto done;
    }

    /* I moved it here in case we need to load a module that
     * needs more heap that is currently allocated.
     * In such case "init_module" for that algorithm forces
     * re-allocation of ext2_wa. It should be safe here b/c the
     * first reference to ext2_wa comes just after and we have 
     * locked ext2_wa before.
     *
     * I know that we may not need the compression at all
     * (compressing 0 or 1 block) but it's better to sacrifice
     * a bit than do make a total mess of this code.
     *
     * FIXME: Totally separate working areas for reading and writing.
     *      Jan R.
     */

    meth = ei->i_compr_method;
    assert(meth < EXT2_N_METHODS);
    alg = ext2_method_table[meth].alg;
#ifdef CONFIG_KMOD
    if (!ext2_algorithm_table[alg].avail) {
	char str[32];

	sprintf(str, "ext2-compr-%s", ext2_algorithm_table[alg].name);
	request_module(str);
    }
#endif

    result = -EINTR;

    /*
     *  Try to get the working area.
     */
#ifdef EXT2_COMPR_REPORT_WA
    printk(KERN_DEBUG "pid %d enters critical region\n", current->pid);
#endif
    if (get_cpu_var(ext2_wr_wa) == NULL)
    {
	ext2_alloc_wr_wa();
    }
    assert(__get_cpu_var(ext2_wr_wa) != NULL);	


    /*
     * Now, we try to compress the cluster.  If the cluster does
     *    not compress well, we just give up.  Otherwise, we reuse
     *    the old blocks to store the compressed data (except that
     *    compressed data is contiguous in the file even if the
     *    uncompressed data had holes).
     */

    /*
     *  Compute the block bitmap, how many bytes of data we have
     *    in the cluster, and the maximum interesting length after
     *    compression.  The bitmap will be used to reallocate blocks
     *    when decompressing the cluster, so that we don't create blocks
     *    that were previously missing.  We also pack the buffers
     *    together.
     */

     head = (struct ext2_cluster_head *) __get_cpu_var(ext2_wr_wa)->c;
#if EXT2_MAX_CLUSTER_BLOCKS > 32
# error "We need to zero more bits than this."
#endif
    *(u32 *) (&head[1]) = 0;
    last_hole_pos = (unsigned) (-1);
    assert(head->holemap[0] == 0);
    assert(head->holemap[1] == 0);
    assert(head->holemap[2] == 0);
    assert(head->holemap[3] == 0);
    assert(*(u32 *) head->holemap == 0);
    assert(count_bits(head->holemap, 4) == 0);

    /* TODO: Check that i_size can't change beneath us.
       do_truncate() is safe because it uses i_sem around changing
       i_size.  For the moment, I do a runtime check. */

    saved_isize = inode->i_size;

#ifdef EXT2_COMPR_REPORT_VERBOSE
    printk
	("00 ext2_compress_cluster[%u]: i_size=%u, s_blocksize_bits=%u, s_nblk=%u\n",
	 __LINE__, (unsigned) inode->i_size, inode->i_sb->s_blocksize_bits,
	 s_nblk);
#endif
//      assert (ROUNDUP_RSHIFT(inode->i_size, inode->i_sb->s_blocksize_bits)
//              >= s_nblk);
    /* This initial guess at ulen doesn't take holes into account
       unless they're at end of cluster.  We ,compensate for other
       holes` during the loop below. */
    ulen = MIN(s_nblk << inode->i_sb->s_blocksize_bits,
	       inode->i_size - ext2_cluster_offset(inode, cluster));
    r_nblk = (((ulen - 1) >> inode->i_sb->s_blocksize_bits) + 1);
    if (r_nblk <= 1) {
	/* MW: required to remove Z flag, otherwise compress 
	 * is tried on each access */
	result = 0;
	goto no_compress;
    }
    /* Verify if more than 1 block to compress in the cluster       */
    nbh = 0;
    for (x = 0; x < s_nblk; x++) {
	if ((bh[x] != NULL) && (bh[x]->b_blocknr != 0)) {
	    nbh++;
	} else {
	    last_hole_pos = x;
	    head->holemap[x >> 3] |= 1 << (x & 7);
	    ulen -= inode->i_sb->s_blocksize;
	    /* impl: We know that it's a whole block because
	       ext2_get_cluster_blocks trims s_nblk for trailing
	       NULL blocks, and partial blocks only come at
	       the end, so there can't be partial NULL blocks. */
	}
    }
    /* We don't try to compress cluster that only have one block
       or no block at all.  (When fragments are implemented, this code
       should be changed.) */
    if (nbh <= 1) {
	/* MW: required to remove Z flag, otherwise compress 
	 * is tried on each access */
	goto no_compress;
    }

    u_nblk = nbh;
    /* Copy the data in the compression area        */
    dst =  __get_cpu_var(ext2_wr_wa)->u;
    for (x = 0; x < s_nblk; x++) {
	if ((bh[x] != NULL) && (bh[x]->b_blocknr != 0)) {
	    restore_b_data_himem(bh[x]);
	    memcpy(dst, bh[x]->b_data, bh[x]->b_size);
	    dst += bh[x]->b_size;
	}
    }

    assert(count_bits(head->holemap, 4) == s_nblk - u_nblk);

#if EXT2_GRAIN_SIZE != EXT2_MIN_BLOCK_SIZE
# error "this code ought to be changed"
#endif

    /* ,maxlen` is the maximum length that the compressed data can
       be while still taking up fewer blocks on disk. */
    holemap_nbytes = (last_hole_pos >> 3) + 1;
    /* impl: Remember that ,last_hole_pos` starts off as being -1,
       so the high 3 bits of ,last_hole_pos >> 3` can be wrong.
       This doesn't matter if holemap_nbytes discards the high
       bits. */

    assert(sizeof(holemap_nbytes) < sizeof(unsigned));
    assert((last_hole_pos == (unsigned) -1)
	   == (holemap_nbytes == 0));
    maxlen =
	((((r_nblk <
	    u_nblk) ? r_nblk : u_nblk) - 1) * inode->i_sb->s_blocksize -
	 sizeof(struct ext2_cluster_head)
	 - holemap_nbytes);
    clen = 0;
    /* Handling of EXT2_AUTO_METH at the moment is just that we
       use the kernel default algorithm.  I hope that in future
       this can be extended to the kernel deciding when to
       compress and what algorithm to use, based on available disk
       space, CPU time, algorithms currently used by the fs,
       etc. */
    if ((meth == EXT2_AUTO_METH)
	|| !ext2_algorithm_table[alg].avail) {
	meth = EXT2_DEFAULT_COMPR_METHOD;
	alg = ext2_method_table[meth].alg;
	assert(ext2_algorithm_table[alg].avail);
    }
    if (alg == EXT2_NONE_ALG)
	goto no_compress;

    clen = ext2_algorithm_table[alg].compress(__get_cpu_var(ext2_wr_wa)->u,
	 __get_cpu_var(ext2_wr_wa)->c + sizeof(struct ext2_cluster_head) + holemap_nbytes,
	 __get_cpu_var(ext2_wr_wa)->heap, ulen, maxlen, ext2_method_table[meth].xarg);

#ifdef EXT2_COMPR_REPORT_ALGORITHMS
    printk(KERN_DEBUG "03 ext2: %lu: cluster %d+%d [%d] compressed "
	   "into %d bytes (ulen = %d, maxlen = %d)\n",
	   inode->i_ino,
	   ext2_cluster_offset(inode, cluster),
	   ext2_cluster_nblocks(inode, cluster),
	   u_nblk, clen, ulen, maxlen);
#endif

    if ((clen == 0) || (clen > maxlen)) {
      no_compress:

	/* this chunk didn't compress. */
	assert(inode->i_size == saved_isize);
#ifdef EXT2_COMPR_REPORT_WA
	printk(KERN_DEBUG
	       "pid %d leaves critical region, nbh=%d, u_nblk=%d, "
	       "inode->i_size=%lu, saved_isize=%lu, clen=%d, ulen=%d, maxlen=%d\n",
	       current->pid, nbh, u_nblk,
	       (long unsigned) inode->i_size, saved_isize, clen, ulen,
	       maxlen);
#endif

	result = 0;
	put_cpu_var(ext2_wr_wa);
	goto done;
    }


#if EXT2_MAX_CLUSTER_BLOCKS > 32
# error "We need to zero more bits than this."
#endif
    assert(-1 <= (int) last_hole_pos);
    assert((int) last_hole_pos < 32);
    assert((le32_to_cpu(*(u32 *) head->holemap)
	    & (~0u << (1 + last_hole_pos))
	    & (~(~0u << (8 * holemap_nbytes))))
	   == 0);
    /* Don't change "~0u << (1 + last_hole_pos)" to "~1u << last_hole_pos" 
       as I almost did, as last_hole_pos can be -1 and cannot be 32. */
    assert(count_bits(head->holemap, holemap_nbytes) == s_nblk - u_nblk);

    /* Compress the blocks at the beginning of the cluster  */
    for (x = 0, nbh = 0; x < s_nblk; x++) {
	if ((bh[x] != NULL) && (bh[x]->b_blocknr != 0)) {
	    if (nbh != x) {
		restore_b_data_himem(bh[x]);
		bh[nbh]->b_blocknr = bh[x]->b_blocknr;
		set_bit(BH_Mapped, &bh[nbh]->b_state);
		bh[x]->b_blocknr = 0;
		assert(buffer_mapped(bh[x]));
		clear_bit(BH_Mapped, &bh[x]->b_state);
	    }
	    nbh++;
	}
    }
    assert(nbh == u_nblk);
    assert(count_bits(head->holemap, holemap_nbytes) == s_nblk - u_nblk);

    /*
     * Compression was successful, so add the header and copy to blocks.
     */

    /* Header. */
    {
	head->magic = cpu_to_le16(EXT2_COMPRESS_MAGIC_04X);
	head->method = meth;
	head->holemap_nbytes = holemap_nbytes;
	head->ulen = cpu_to_le32(ulen);
	head->clen = cpu_to_le32(clen);

	barrier(); //mw: "barrier" tells compiler not to re-order resulting asm statments, somehow.
	head->checksum =
	    cpu_to_le32(ext2_adler32
			(le32_to_cpu(*(u32 *) __get_cpu_var(ext2_wr_wa)->c),
			 __get_cpu_var(ext2_wr_wa)->c + 8,
			 (sizeof(struct ext2_cluster_head) - 8 +
			  head->holemap_nbytes + clen)));
    }

    assert((le32_to_cpu(*(u32 *) head->holemap)
	    & (~0 << (1 + last_hole_pos))
	    & ((1 << (8 * holemap_nbytes)) - 1)) == 0);
    result = clen += sizeof(struct ext2_cluster_head) + holemap_nbytes;
    c_nblk = ROUNDUP_RSHIFT(clen, inode->i_sb->s_blocksize_bits);

    /* Release unneeded buffer heads.  (Freeing is done later,
       after unlocking ext2_wr_wa.) */
    assert(nbh == u_nblk);
    nbh = c_nblk;

#ifdef  EXT2_COMPR_REPORT
    trace_e2c("ext2_compress_cluster: head->clen=%d, clen=%d\n", head->clen, clen);
#endif
    src = __get_cpu_var(ext2_wr_wa)->c;

    for (n = 0; (int) clen > 0; n++) {
	restore_b_data_himem(bh[n]);
	if (clen >= inode->i_sb->s_blocksize) {
	    memcpy(bh[n]->b_data, src, inode->i_sb->s_blocksize);
	} else {
	    memcpy(bh[n]->b_data, src, clen);
	}

	/* TO_DO: OSYNC.  means: write opertions are blocking until the
	 * the pages are written from page cache to disk */

	set_buffer_uptodate(bh[n]);
	set_buffer_dirty(bh[n]);
	src += inode->i_sb->s_blocksize;
	clen -= inode->i_sb->s_blocksize;
    }

    i = 0;
    assert(n == c_nblk);
    assert((le32_to_cpu(*(u32 *) head->holemap)
	    & (~0 << (1 + last_hole_pos))
	    & ((1 << (8 * holemap_nbytes)) - 1)) == 0);

    /* Runtime check that no-one can change i_size while i_sem is down.
       (See where saved_isize is set, above.) */
    assert(inode->i_size == saved_isize);
    assert(!mapping_mapped(inode->i_mapping));

    /* Free the remaining blocks, and shuffle used blocks to start
       of cluster in blkaddr array. */
    {
	u32 free_ix, curr;
	int err;

	/* Calculate free_ix.  There should be ,c_nblk`
	   non-hole blocks among the first ,free_ix`
	   blocks. */
	{
	    assert((le32_to_cpu(*(u32 *) head->holemap)
		    & (~0 << (1 + last_hole_pos))
		    & ((1 << (8 * holemap_nbytes)) - 1)) == 0);
	    assert(n == c_nblk);
	    for (free_ix = 0;
		 ((int) free_ix <= (int) last_hole_pos) && (n > 0);
		 free_ix++)
		if (!(head->holemap[free_ix >> 3]
		      & (1 << (free_ix & 7))))
		    n--;
	    free_ix += n;

	    if ((free_ix < c_nblk)
		|| (free_ix + u_nblk > s_nblk + c_nblk)
		|| (free_ix >= ext2_cluster_nblocks(inode, cluster))
		|| ((holemap_nbytes == 0) && (c_nblk != free_ix))) {
		assert(free_ix >= c_nblk);
		/*assert (free_ix - c_nblk <= s_nblk - u_nblk); */
		assert(free_ix + u_nblk <= s_nblk + c_nblk);
		assert(free_ix < ext2_cluster_nblocks(inode, cluster));
		assert((holemap_nbytes != 0) || (c_nblk == free_ix));
		assert(1 <= c_nblk);
		assert(c_nblk < u_nblk);
		assert(u_nblk <= s_nblk);
		assert(s_nblk <= ext2_cluster_nblocks(inode, cluster));
		assert(ext2_cluster_nblocks(inode, cluster) <=
		       EXT2_MAX_CLU_NBLOCKS);
		ext2_error(inode->i_sb, "ext2_compress_cluster",
			   "re assertions: c=%d, u=%d, f=%d, s=%d, n=%d, "
			   "lhp=%d, hm=%x, hnb=%d, " "ino=%lu, clu=%u",
			   (int) c_nblk, (int) u_nblk, (int) free_ix,
			   (int) s_nblk, (int) ext2_cluster_nblocks(inode,
								    cluster),
			   (int) last_hole_pos,
			   (unsigned) le32_to_cpu(*(u32 *) head->holemap),
			   (int) holemap_nbytes, inode->i_ino, cluster);
	    }
	}
        
	/*mw: put here: set all __get_cpu related pointers to NULL
	      as they become invalid with put_cpu */
    	head = NULL;		/* prevent any more stupid bugs */
    	src = NULL;
        dst = NULL;
    	put_cpu_var(ext2_wr_wa);

#ifdef EXT2_COMPR_DEBUG
	/* TODO: remove this TEST */
        /* mw: ext2_free_cluster_blocks can sleep: check we are not atomic */
	schedule();
#endif

	/* Free unneeded blocks, and mark cluster as
	   compressed. */
	err = ext2_free_cluster_blocks
	    (inode,
	     ext2_cluster_block0(inode, cluster) + free_ix,
	     ext2_cluster_nblocks(inode, cluster) - free_ix);
	/* pjm 1998-06-15: This should help reduce fragmentation.
	   Actually, we could set block to clu_block0 + clu_nbytes,
	   and goal to the last allocated blkaddr in the compressed
	   cluster.
	   It would be nice if we would transfer the freed blocks
	   to preallocation, while we're at it. */
//              write_lock(&ei->i_meta_lock);
	/* mw: i_next_alloc_goal and i_next_alloc_block were removed in 2.6.24.x
	 *     so we dont need to set them to 0 (they are anyway, somehow).
	 */
	//ei->i_next_alloc_goal = ei->i_next_alloc_block = 0;
//              write_unlock(&ei->i_meta_lock);
	if (err < 0) {
		goto done;
	}
	/* Note that ext2_free_cluster_blocks() marks the
	   cluster as compressed. */

	/* Shuffle used blocks to beginning of block-number array. */
	{
	    struct ext2_bkey key;
	    unsigned i;

	    if (!ext2_get_key(&key,
			      inode,
			      ext2_cluster_block0(inode, cluster))) {
		ei->i_flags |= EXT2_ECOMPR_FL;
		result = -EIO;
		free_ix = 0;
	    }
	    for (i = 0; i < free_ix; i++) {
		curr = ext2_get_key_blkaddr(&key);

		if ((c_nblk == free_ix)
		    && (curr != bh[i]->b_blocknr)) {
		    /* "Can't happen", yet has
		       happened a couple of times. */
		    ext2_error(inode->i_sb, "ext2_compress_cluster",
			       "c_nblk=free_ix=%d, "
			       "curr=%u, b_blocknr=%lu, "
			       "lhp=%d , hm=<noinfo>, "
			       "ino=%lu, blk=%u",
			       c_nblk, curr,
			       (unsigned long) bh[i]->b_blocknr,
			       (int) last_hole_pos,
			       /*mw: became invalid due put_cpu:
				(unsigned) le32_to_cpu(*(u32 *) head->
						      holemap),*/
			       inode->i_ino,
			       (unsigned) 
			       ext2_cluster_block0(inode, cluster) + i);
		}
		err = ext2_set_key_blkaddr(&key,
					   (i < c_nblk
					    ? bh[i]->b_blocknr
					    : EXT2_COMPRESSED_BLKADDR));
		if (err)
		    break;
		if (!ext2_next_key(&key, 1)) {
		    ei->i_flags |= EXT2_ECOMPR_FL;	/* sorry... */
		    result = -EIO;
		    break;
		}
	    }
	    ext2_free_key(&key);
	}
    }

    /*
     *        Unlock the working area.
     */

#ifdef EXT2_COMPR_REPORT_WA
    printk(KERN_DEBUG "pid %d leaves critical region\n", current->pid);
#endif

    assert(c_nblk < u_nblk);
    ext2_mark_algorithm_use(inode, alg);

    /* TLL update b_assoc_map per 2.6.20 6-07-07 */
    for (i = 0; i < c_nblk; i++)
	if (bh[i] != NULL) {
	    bh[i]->b_assoc_map = inode->i_mapping;
	    bh[i]->b_page->mapping = inode->i_mapping;	//Andreas 5-24-07 : necessary? WRONG?
	}
    //mw: we must force the writeback, otherwise ext2_readpage will get confused
    //    yaboo ding had similiar code above. but I think it makes more sense after
    //    the block shuffeling.
    //    Note: generic_oysnc_inode() made trouble with USB-Sticks and caused a lot
    //    of IO, stalled system ... therefore ll_rw_block() replace it. Anyway we already operate 
    //        with this low-level function. 

    /*mw: new "hole" fix. hole == bdev bug! */
    for (i = 0; i < c_nblk; i++) {

	/* this was a hole (uncompressed)
	 * at the beginning of the cluster.
	 * so NO block was yet associated with it.
	 * But now we need it, because a compressed
	 * cluster always starts at the cluster.*/
	if (!buffer_mapped(bh[i]) || bh[i]->b_bdev == NULL) {
	    u32 block = ext2_cluster_block0(inode, cluster);
	    ext2_get_block(inode, block + i, bh[i], 1);
	    //printk("ext2_get_block Block:%lu, Mapped:%i, Page:%lu, bdev: %#x\n", bh[i]->b_blocknr, (bh[i]->b_state & BH_Mapped), (bh[i]->b_page ? bh[i]->b_page->index : 0), bh[i]->b_bdev );
	}
	assert(buffer_mapped(bh[i]));
	assert(bh[i]->b_bdev != NULL);
	assert(bh[i]->b_bdev == inode->i_sb->s_bdev);
    }

    ll_rw_block(WRITE, c_nblk, bh);

    CHECK_NOT_ATOMIC
    //mw: seems we have to wait here, otherwise: crash!
    for (i = 0; i < c_nblk; i++) {
	if (bh[i])
	    wait_on_buffer(bh[i]);
	//printk("written compressed block: Block:%lu, Mapped:%i, Page:%lu, bdev: %#x\n", bh[i]->b_blocknr, (bh[i]->b_state & BH_Mapped), (bh[i]->b_page ? bh[i]->b_page->index : 0), bh[i]->b_bdev );
    }


#ifdef CONFIG_HIGHMEM
    if (kmapped)
	ext2_kunmap_cluster_pages(NULL, pg, NULL);
#endif

    inode->i_ctime = CURRENT_TIME;	//mw: these two come always together. So I also put it here.
    mark_inode_dirty_sync(inode);

    //ext2_update_inode(inode, inode_needs_sync(inode)); //mw: might be able to fix pipe_write vs. readpage. mutex-rec-locking

    /* COMPRBLK is already high, so no need to raise it. */
    {
	for (i = c_nblk; (i < EXT2_MAX_CLUSTER_BLOCKS) && (bh[i] != NULL);
	     i++) {
	    clear_buffer_dirty(bh[i]);
	    bh[i]->b_blocknr = 0;
	    clear_bit(BH_Mapped, &bh[i]->b_state);
	    clear_bit(BH_Uptodate, &bh[i]->b_state);
	}
	for (i = 0; i < EXT2_MAX_CLUSTER_PAGES; i++) {
	    if (pg[i] == NULL) {
		break;
	    }
	    assert(PageLocked(pg[i]));
	    ClearPageUptodate(pg[i]);
	    unlock_page(pg[i]);
	    page_cache_release(pg[i]);
	}

	/* invalidate_inode_buffers replacement code: TLL 02/21/07
	 * e2compr on post 2.6.10 kernels do not have an uptodate
	 * mapping->assoc_mapping (other Vm(?) changes require it be
	 * made explicit, 2.4 kernels have it implicit). Therefore, when
	 * umount is called, a GPF ensues from a NULL ops pointer.
	 * e2c on a USB thumbdrive mounted as the root fs does not
	 * support repeated compress/uncompress cycles on a given file.
	 * Inlined the flush list code to explicityly force update to
	 * disk with a known valid bh list.
	 */

	/* mw: I consider this code as ... not so good! */
	/*	  
	  if (inode_has_buffers(inode)) {
		//struct address_space *mapping = &inode->i_data;
		// struct address_space *buffer_mapping = mapping->assoc_mapping;
		// requires: inode->i_data->mapping->assoc_mapping; to be set
		invalidate_inode_buffers(inode);	// TLL do it proper 5-25-07
		//if (dotrunc)
		 //ext2_truncate(inode);	// TLL file size hack 6-19-07 
	  }
	*/

    }
#ifdef  EXT2_COMPR_REPORT
    trace_e2c(" < < < ext2_compress_cluster %i: [done cpr] inode=%ld\n", cluster, inode->i_ino);
#endif
    return result;


  done:

#ifdef CONFIG_HIGHMEM
    if (kmapped)
	ext2_kunmap_cluster_pages(NULL, pg, NULL);
#endif

    {
	for (i = 0; i < EXT2_MAX_CLUSTER_PAGES; i++) {
	    if (pg[i] == NULL) {
		break;
	    }
	    unlock_page(pg[i]);
	    page_cache_release(pg[i]);
	}
	/* TLL cp to compr dir bug fix 03-25-07
	   Truncate uncompressed files to their uncompressed
	   length, i.e. force kernel to update inode and sb */
	//if(dotrunc)
	//26.08.2011: ext2_truncate(inode) does not exist anymore
	ext2_truncate_blocks(inode, inode->i_size);
	
    }
#ifdef EXT2_COMPR_REPORT_VERBOSE
    {
	int i;

	printk(KERN_DEBUG "ext2_compress_cluster[end]: buffers kept for cluster=%d\n", cluster);
	for (i = 0; i < nbh; i++) {
	    if (bh[i]) {
		printk(KERN_DEBUG "2buffer_head[%d]: blocknr=%lu, addr=0x%p ", i, (unsigned long) bh[i]->b_blocknr, bh[i]);
		if (bh[i]->b_page)
		    printk(KERN_DEBUG "2:[page->index=%ld]\n", bh[i]->b_page->index);
		else
		    printk(KERN_DEBUG "[No page]\n");
	    } else
		printk(KERN_DEBUG "buffer_head[%d] is NULL\n", i);
	}
    }
#endif

#ifdef  EXT2_COMPR_REPORT
    trace_e2c(" < < < ext2_compress_cluster %i: [done NO cpr] inode=%ld\n", cluster, inode->i_ino);
#endif
    return result;
}


/* Go through all the clusters and compress them if not already
   compressed.

   This is called by ext2_put_inode() and ext2_release_file().  Later,
   we may have ext2_ioctl() call it (when EXT2_COMPR_FL rises).  None
   of the callers does any locking, so we do it here.

   Neither of the current callers uses the return code, but we get ready
   for if we start using it.

   Returns 0 on "success" (whether or not we cleared EXT2_CLEANUP_FL
   or EXT2_DIRTY_FL bits), -errno on error. */
int ext2_cleanup_compressed_inode(struct inode *inode)
{
    u32 cluster;
    u32 n_clusters;
    int dirty = 0;
    int err = 0;
    u32 comprblk_mask;
    atomic_t start_i_count = inode->i_count;
    int retry = 0;
    int have_downed;
    struct ext2_inode_info *ei = EXT2_I(inode);
#ifdef EXT2_COMPR_REPORT
    char bdn[BDEVNAME_SIZE];
#endif

    /* impl: Actually, this assertion could fail if the kernel
       isn't locked.  I haven't looked, but I suppose that the
       kernel always is locked when this is called. */
    assert(ei->i_compr_flags & EXT2_CLEANUP_FL);

#ifdef EXT2_COMPR_REPORT_PUT
    printk(KERN_DEBUG "ext2_cleanup_compressed_inode() called for pid %d; "
	   "dev=%s, ino=%lu, i_state=0x%lx, i_count=%u\n",
	   current->pid, bdevname(inode->i_sb->s_bdev, bdn), inode->i_ino,
	   inode->i_state, atomic_read(&inode->i_count));
#endif

    /* Do these tests twice: once before down() and once after. */
    for (have_downed = 0;; have_downed++) {
	if ((ei->i_flags & (EXT2_COMPR_FL | EXT2_DIRTY_FL))
	    != (EXT2_COMPR_FL | EXT2_DIRTY_FL)) {
	    if (have_downed)
		goto out;
	    /* TLL 5-25-07 changed from a warning to trace */
	    /*trace_e2c("ext2_cleanup_compressed_inode: trying to un/compress an "
	       "uncompressable file.\n"
	       "i_flags=%#x. (dev=%s, ino=%lu, down=%d)\n",
	       ei->i_flags, bdevname(inode->i_sb->s_bdev, bdn), 
	       inode->i_ino, have_downed); */
	    return 0;
	}

	/* test if file is mapped by mmap */
	if (mapping_mapped(inode->i_mapping))
	{
	    //trace_e2c("ext2_cleanup_compressed_inode: (dev. %s): ino=%ld: file mapped, does not compress cluster\n", bdevname(inode->i_sb->s_bdev, bdn), inode->i_ino);
	    if (have_downed)
		goto out;
	    else
		return 0;
	}

	if (IS_RDONLY(inode)
	    || (ei->i_flags & EXT2_ECOMPR_FL)) {
	    ei->i_compr_flags &= ~EXT2_CLEANUP_FL;
	    if (have_downed)
		goto out;
	    else
		return 0;
	}

	//mw            
	if (ext2_get_dcount(inode) > 1) {
	    err = 0;
	    //printk("Compress: file busy (dcount: %i>1)\n", ext2_get_dcount(inode));
	    if (have_downed)
		goto out;
	    else
		return 0;
	}

	if (have_downed)
	    break;

	/* Quotas aren't otherwise kept if file is opened O_RDONLY. */
	dquot_initialize(inode);
	
	/* Check whether OSYNC of inode is acutally running */
	//if (ei->i_compr_flags & EXT2_OSYNC_INODE)
	//printk(KERN_DEBUG "OSYNC!\n");

	/* I think:
	 * checking these flags should prevent that one Process aquires the MUTEX again, 
	 * e.g. in a recursive call
	 * BUT: what happens acutally: two processes are working on this inode: pdflush and the userprogramm
	 * SO: the check might be correct if:  ei->i_compr_flags & EXT2_OSYNC_INOD AND the same process already posesses this lock!!!
	 */
	//if (!(ei->i_compr_flags & EXT2_OSYNC_INODE))
	//{
	mutex_lock(&inode->i_mutex);
#ifdef EXT2_COMPR_REPORT_MUTEX
	printk(KERN_DEBUG "CLEANUP_LOCK of PID %u @ inode:%lu\n", current->pid, inode->i_ino);
#endif
	//}
    }
    n_clusters = ext2_n_clusters(inode);

#ifdef EXT2_COMPR_REPORT_PUT
    printk(KERN_DEBUG "ext2: inode:%lu: put compressed, clusters = %d, flags = %x, pid = %u\n",
	   inode->i_ino, n_clusters, ei->i_flags, current->pid);
#endif

    assert(atomic_read(&inode->i_mutex.count) <= 0);	/* i.e. mutex_lock */

    /* Try to compress the clusters.  We clear EXT2_DIRTY_FL only
       if we looked at every cluster and if there was no error.  */

    /* impl: We raise EXT2_COMPRBLK_FL now so that ext2_ioctl()
       doesn't try to change the cluster size beneath us.  If need
       be, we restore the bit to its original setting before
       returning.  Note that no-one else can _change_
       EXT2_COMPRBLK_FL while we work because i_sem is down. */
    /* impl: Note what's happening here with comprblk_mask.  The
       current state of COMPRBLK_FL (before we start) is that
       (comprblk == 1) || (no compressed clusters).  At the end of
       the procedure, comprblk == one if (at least one compressed
       cluster, or an error occurred preventing us from finding
       out). */
    comprblk_mask = ~EXT2_COMPRBLK_FL | ei->i_flags;
    ei->i_flags |= EXT2_COMPRBLK_FL;

    for (cluster = 0; cluster < n_clusters; cluster++) {
	if (atomic_read(&inode->i_count) > atomic_read(&start_i_count)) {
	    /* This is a poor way of doing this (and doubly
	       poor now that the only users of i_count are
	       the dentries), but the idea is not to
	       compress things tht are likely to be
	       decompressed soon.  I guess a better way of
	       doing this would be just to make sure tht
	       the stuff is in the page cache. */
	    retry = 1;
	    break;
	}
	err = ext2_cluster_is_compressed_fn(inode, cluster);
	if (err == 0) {
	    //mw: ext2_compress_cluster might clean EXT2_COMPRBLK_FL, therefore raise it for every new cluster
	    ei->i_flags |= EXT2_COMPRBLK_FL;

	    err = ext2_compress_cluster(inode, cluster);
	    if (err < 0)
		dirty = 1;
	    else if (err > 0)
		comprblk_mask = ~0ul;
	} else if (err < 0)
	    break;
	else {
	    err = 0;
	    assert(comprblk_mask == ~0ul);	/* i.e. that EXT2_COMPRBLK_FL was high. */
	}
    }

    if ((cluster >= n_clusters) && !dirty)
	ei->i_flags &= ~EXT2_DIRTY_FL;
    if (!retry) {
	ei->i_compr_flags &= ~EXT2_CLEANUP_FL;
	ei->i_flags &= comprblk_mask;
    }

    /* We clear EXT2_CLEANUP_FL because, otherwise, we'll get
       called again almost immediately. */

    /*
     *  The CLEANUP flag *MUST* be cleared, otherwise the iput routine
     *  calls ext2_put_inode() again (because i_dirt is set) and there
     *  is a loop.  The control scheme (CLEANUP + DIRTY flags) could 
     *  probably be improved.  On the other hand, i_dirt MUST be set
     *  because we may have sleeped, and we must force the iput routine
     *  to look again at the i_count ...
     */
    /* TODO: Have a look at this cleanup scheme.  The above
       comment sounds wrong. */

    inode->i_ctime = CURRENT_TIME;
    mark_inode_dirty_sync(inode);
  out:

#ifdef EXT2_COMPR_REPORT_MUTEX
    printk(KERN_DEBUG "CLEANUP_UNLOCK of PID %u @ inode:%lu\n", current->pid, inode->i_ino);
#endif

//      if (!(ei->i_compr_flags & EXT2_OSYNC_INODE)) {  /* MW 5-16-07 */
    mutex_unlock(&inode->i_mutex);
//      }       /* MW 5-16-07 */                  
    return err;			/* TODO: Check that ,err` is appropriate. */
}


int ext2_recognize_compressed(struct inode *inode, unsigned cluster)
{
    /* ext2_recognize_compressed(): Check tht the cluster is valid
       in every way, and then do the EXT2_COMPRESSED_BLKADDR
       thing. */
    /* nyi, fixme.  All of the userspace stuff (EXT2_NOCOMPR_FL
       etc.) needs work, so I might as well leave this.  See
       ioctl.c for a description of what it's supposed to do. */
    return -ENOSYS;
}


/* Look for compressed clusters.  If none, then clear EXT2_COMPRBLK_FL.

   Called by:
       ext2_truncate().
       */
void ext2_update_comprblk(struct inode *inode)
{
    unsigned block, last_block;
    struct ext2_bkey key;
    struct ext2_inode_info *ei = EXT2_I(inode);

    assert(ei->i_flags & EXT2_COMPRBLK_FL);
    if (inode->i_size == 0) {
	ei->i_flags &= ~EXT2_COMPRBLK_FL;
	trace_e2c("ext2_update_comprblk 1: inode: %lu removed EXT2_COMPRBLK_FL!\n", inode->i_ino);
	return;
    }
    last_block = ROUNDUP_RSHIFT(inode->i_size,
				inode->i_sb->s_blocksize_bits) - 1;
    block = ext2_first_cluster_nblocks(inode) - 1;

    assert(atomic_read(&inode->i_mutex.count) <= 0);

    if (!ext2_get_key(&key, inode, block))
	return;
    for (;;) {
	if (ext2_get_key_blkaddr(&key) == EXT2_COMPRESSED_BLKADDR)
	    goto out;
	if (block >= last_block)
	    goto clear;
	if (!ext2_next_key(&key, ei->i_clu_nblocks))
	    goto out;
	block += ei->i_clu_nblocks;
    }
  clear:
    trace_e2c("ext2_update_comprblk 2: inode: %lu removed EXT2_COMPRBLK_FL!\n", inode->i_ino);
    ei->i_flags &= ~EXT2_COMPRBLK_FL;
  out:
    ext2_free_key(&key);
    assert(atomic_read(&inode->i_mutex.count) <= 0);

}


/*
 * allocate working areas
 */

DEFINE_PER_CPU(struct ext2_wa_S *, ext2_rd_wa) = NULL;
DEFINE_PER_CPU(struct ext2_wa_S *, ext2_wr_wa) = NULL;

/* SMP, setup wa's. caller must hold wa already via get_cpu_var */
void ext2_alloc_rd_wa(){
	if ((__get_cpu_var(ext2_rd_wa) == NULL) ) {
		size_t rsize =  2 * EXT2_MAX_CLUSTER_BYTES; //mw: just guessing

		__get_cpu_var(ext2_rd_wa) = vmalloc (rsize);
		if (__get_cpu_var(ext2_rd_wa) == NULL)
			printk ("EXT2-fs: can't allocate working area; compression turned off.\n");
		else {
			printk ("ext2-compression: allocated read buffer for CPU%i at %p-%p (%zu bytes)\n",
				get_cpu(), __get_cpu_var(ext2_rd_wa), (char *)__get_cpu_var(ext2_rd_wa) + rsize, rsize);
#  ifdef EXT2_COMPR_REPORT_WA
			printk (KERN_INFO "EXT2-fs: rd_wa=%p--%p (%d)\n",
				ext2_rd_wa, (char *)ext2_rd_wa + rsize, rsize);
#  endif
			put_cpu();
		}
	}
}

void ext2_alloc_wr_wa(){

	if ((__get_cpu_var(ext2_wr_wa) == NULL) ) {
		size_t wsize = 2 * EXT2_MAX_CLUSTER_BYTES; //mw: just guessing
		__get_cpu_var(ext2_wr_wa) = vmalloc (wsize);

		if (__get_cpu_var(ext2_wr_wa) == NULL)
			printk ("EXT2-fs: can't allocate working area; "
				"compression turned off.\n");
		else {
			printk ("ext2-compression: allocated write buffer for CPU%i at %p-%p (%zu bytes)\n",
				get_cpu(), __get_cpu_var(ext2_wr_wa), (char *)__get_cpu_var(ext2_wr_wa) + wsize, wsize);
#ifdef EXT2_COMPR_REPORT_WA
			printk (KERN_INFO "EXT2-fs: wr_wa=%p--%p (%d)\n",
				ext2_wr_wa, (char *)ext2_wr_wa + wsize, wsize);
#endif
			put_cpu();	
		}
	}
}


