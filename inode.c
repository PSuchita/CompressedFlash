/*
 *  linux/fs/ext2/inode.c
 *
 * Copyright (C) 1992, 1993, 1994, 1995
 * Remy Card (card@masi.ibp.fr)
 * Laboratoire MASI - Institut Blaise Pascal
 * Universite Pierre et Marie Curie (Paris VI)
 *
 *  from
 *
 *  linux/fs/minix/inode.c
 *
 *  Copyright (C) 1991, 1992  Linus Torvalds
 *
 *  Goal-directed block allocation by Stephen Tweedie
 * 	(sct@dcs.ed.ac.uk), 1993, 1998
 *  Big-endian to little-endian byte-swapping/bitmaps by
 *        David S. Miller (davem@caip.rutgers.edu), 1995
 *  64-bit file support on 64-bit platforms by Jakub Jelinek
 * 	(jj@sunsite.ms.mff.cuni.cz)
 *
 *  Assorted race fixes, rewrite of ext2_get_block() by Al Viro, 2000
 */

#include <linux/time.h>
#include <linux/highuid.h>
#include <linux/pagemap.h>
#include <linux/quotaops.h>
#include <linux/writeback.h>
#include <linux/buffer_head.h>
#include <linux/mpage.h>
#include <linux/fiemap.h>
#include <linux/namei.h>
#ifdef CONFIG_EXT2_COMPRESS
#include <linux/kmod.h>
#include <linux/ext2_fs_c.h>
#include <linux/spinlock.h>
#include <linux/pagevec.h>
#else
#include "ext2.h"
#endif
#include "acl.h"
#include "xip.h"

#ifdef CONFIG_EXT2_COMPRESS
/* mw: this function counts all references 
 * to this inode. this is necessary to
 * refuse un/compression if the file has
 * more than one refernce, I guess. */ 
int ext2_get_dcount(struct inode *inode)
{        
	struct dentry *dentry;
    	struct list_head *head, *next, *tmp;
    	int count;
    	
    	head = &inode->i_dentry;
    	next = inode->i_dentry.next;
    	count = 0;
    	while (next != head) {
    		dentry = list_entry(next, struct dentry, d_alias);
    		tmp = next;
    		next = tmp->next;
		spin_lock(&dentry->d_lock);
    		count += dentry->d_count;
		spin_unlock(&dentry->d_lock);
    		//mw: similar to fs/dcache.c
    	}
        
        return count;
}
#endif

static int __ext2_write_inode(struct inode *inode, int do_sync);

/*
 * Test whether an inode is a fast symlink.
 */
static inline int ext2_inode_is_fast_symlink(struct inode *inode)
{
	int ea_blocks = EXT2_I(inode)->i_file_acl ?
		(inode->i_sb->s_blocksize >> 9) : 0;

	return (S_ISLNK(inode->i_mode) &&
		inode->i_blocks - ea_blocks == 0);
}

#ifndef CONFIG_EXT2_COMPRESS
static void ext2_truncate_blocks(struct inode *inode, loff_t offset);
#endif

static void ext2_write_failed(struct address_space *mapping, loff_t to)
{
	struct inode *inode = mapping->host;

	if (to > inode->i_size) {
		truncate_pagecache(inode, to, inode->i_size);
		ext2_truncate_blocks(inode, inode->i_size);
	}
}

/*
 * Called at the last iput() if i_nlink is zero.
 */
void ext2_evict_inode(struct inode * inode)
{
	struct ext2_block_alloc_info *rsv;
	int want_delete = 0;

	if (!inode->i_nlink && !is_bad_inode(inode)) {
		want_delete = 1;
		dquot_initialize(inode);
	} else {
		dquot_drop(inode);
	}

	truncate_inode_pages(&inode->i_data, 0);

	if (want_delete) {
		/* set dtime */
		EXT2_I(inode)->i_dtime	= get_seconds();
		mark_inode_dirty(inode);
		__ext2_write_inode(inode, inode_needs_sync(inode));
		/* truncate to 0 */
		inode->i_size = 0;
		if (inode->i_blocks)
			ext2_truncate_blocks(inode, 0);
	}

	invalidate_inode_buffers(inode);
	end_writeback(inode);

	ext2_discard_reservation(inode);
	rsv = EXT2_I(inode)->i_block_alloc_info;
	EXT2_I(inode)->i_block_alloc_info = NULL;
	if (unlikely(rsv))
		kfree(rsv);

	if (want_delete)
		ext2_free_inode(inode);
}

typedef struct {
	__le32	*p;
	__le32	key;
	struct buffer_head *bh;
} Indirect;

static inline void add_chain(Indirect *p, struct buffer_head *bh, __le32 *v)
{
	p->key = *(p->p = v);
	p->bh = bh;
}

static inline int verify_chain(Indirect *from, Indirect *to)
{
	while (from <= to && from->key == *from->p)
		from++;
	return (from > to);
}

/**
 *	ext2_block_to_path - parse the block number into array of offsets
 *	@inode: inode in question (we are only interested in its superblock)
 *	@i_block: block number to be parsed
 *	@offsets: array to store the offsets in
 *      @boundary: set this non-zero if the referred-to block is likely to be
 *             followed (on disk) by an indirect block.
 *	To store the locations of file's data ext2 uses a data structure common
 *	for UNIX filesystems - tree of pointers anchored in the inode, with
 *	data blocks at leaves and indirect blocks in intermediate nodes.
 *	This function translates the block number into path in that tree -
 *	return value is the path length and @offsets[n] is the offset of
 *	pointer to (n+1)th node in the nth one. If @block is out of range
 *	(negative or too large) warning is printed and zero returned.
 *
 *	Note: function doesn't find node addresses, so no IO is needed. All
 *	we need to know is the capacity of indirect blocks (taken from the
 *	inode->i_sb).
 */

/*
 * Portability note: the last comparison (check that we fit into triple
 * indirect block) is spelled differently, because otherwise on an
 * architecture with 32-bit longs and 8Kb pages we might get into trouble
 * if our filesystem had 8Kb blocks. We might use long long, but that would
 * kill us on x86. Oh, well, at least the sign propagation does not matter -
 * i_block would have to be negative in the very beginning, so we would not
 * get there at all.
 */

static int ext2_block_to_path(struct inode *inode,
			long i_block, int offsets[4], int *boundary)
{
	int ptrs = EXT2_ADDR_PER_BLOCK(inode->i_sb);
	int ptrs_bits = EXT2_ADDR_PER_BLOCK_BITS(inode->i_sb);
	const long direct_blocks = EXT2_NDIR_BLOCKS,
		indirect_blocks = ptrs,
		double_blocks = (1 << (ptrs_bits * 2));
	int n = 0;
	int final = 0;

	if (i_block < 0) {
		ext2_msg(inode->i_sb, KERN_WARNING,
			"warning: %s: block < 0", __func__);
	} else if (i_block < direct_blocks) {
		offsets[n++] = i_block;
		final = direct_blocks;
	} else if ( (i_block -= direct_blocks) < indirect_blocks) {
		offsets[n++] = EXT2_IND_BLOCK;
		offsets[n++] = i_block;
		final = ptrs;
	} else if ((i_block -= indirect_blocks) < double_blocks) {
		offsets[n++] = EXT2_DIND_BLOCK;
		offsets[n++] = i_block >> ptrs_bits;
		offsets[n++] = i_block & (ptrs - 1);
		final = ptrs;
	} else if (((i_block -= double_blocks) >> (ptrs_bits * 2)) < ptrs) {
		offsets[n++] = EXT2_TIND_BLOCK;
		offsets[n++] = i_block >> (ptrs_bits * 2);
		offsets[n++] = (i_block >> ptrs_bits) & (ptrs - 1);
		offsets[n++] = i_block & (ptrs - 1);
		final = ptrs;
	} else {
		ext2_msg(inode->i_sb, KERN_WARNING,
			"warning: %s: block is too big", __func__);
	}
	if (boundary)
		*boundary = final - 1 - (i_block & (ptrs - 1));

	return n;
}

/**
 *	ext2_get_branch - read the chain of indirect blocks leading to data
 *	@inode: inode in question
 *	@depth: depth of the chain (1 - direct pointer, etc.)
 *	@offsets: offsets of pointers in inode/indirect blocks
 *	@chain: place to store the result
 *	@err: here we store the error value
 *
 *	Function fills the array of triples <key, p, bh> and returns %NULL
 *	if everything went OK or the pointer to the last filled triple
 *	(incomplete one) otherwise. Upon the return chain[i].key contains
 *	the number of (i+1)-th block in the chain (as it is stored in memory,
 *	i.e. little-endian 32-bit), chain[i].p contains the address of that
 *	number (it points into struct inode for i==0 and into the bh->b_data
 *	for i>0) and chain[i].bh points to the buffer_head of i-th indirect
 *	block for i>0 and NULL for i==0. In other words, it holds the block
 *	numbers of the chain, addresses they were taken from (and where we can
 *	verify that chain did not change) and buffer_heads hosting these
 *	numbers.
 *
 *	Function stops when it stumbles upon zero pointer (absent block)
 *		(pointer to last triple returned, *@err == 0)
 *	or when it gets an IO error reading an indirect block
 *		(ditto, *@err == -EIO)
 *	or when it notices that chain had been changed while it was reading
 *		(ditto, *@err == -EAGAIN)
 *	or when it reads all @depth-1 indirect blocks successfully and finds
 *	the whole chain, all way to the data (returns %NULL, *err == 0).
 */
static Indirect *ext2_get_branch(struct inode *inode,
				 int depth,
				 int *offsets,
				 Indirect chain[4],
				 int *err)
{
	struct super_block *sb = inode->i_sb;
	Indirect *p = chain;
	struct buffer_head *bh;

	*err = 0;
	/* i_data is not going away, no lock needed */
	add_chain (chain, NULL, EXT2_I(inode)->i_data + *offsets);
#ifdef CONFIG_EXT2_COMPRESS
	if (HOLE_BLKADDR(p->key))
#else
	if (!p->key)
#endif
		goto no_block;
	while (--depth) {
		bh = sb_bread(sb, le32_to_cpu(p->key));
		if (!bh)
			goto failure;
		read_lock(&EXT2_I(inode)->i_meta_lock);
		if (!verify_chain(chain, p))
			goto changed;
		add_chain(++p, bh, (__le32*)bh->b_data + *++offsets);
		read_unlock(&EXT2_I(inode)->i_meta_lock);
#ifdef CONFIG_EXT2_COMPRESS
		if (HOLE_BLKADDR(p->key))
#else
		if (!p->key)
#endif
			goto no_block;
	}
	return NULL;

changed:
	read_unlock(&EXT2_I(inode)->i_meta_lock);
	brelse(bh);
	*err = -EAGAIN;
	goto no_block;
failure:
	*err = -EIO;
no_block:
	return p;
}

/**
 *	ext2_find_near - find a place for allocation with sufficient locality
 *	@inode: owner
 *	@ind: descriptor of indirect block.
 *
 *	This function returns the preferred place for block allocation.
 *	It is used when heuristic for sequential allocation fails.
 *	Rules are:
 *	  + if there is a block to the left of our position - allocate near it.
 *	  + if pointer will live in indirect block - allocate near that block.
 *	  + if pointer will live in inode - allocate in the same cylinder group.
 *
 * In the latter case we colour the starting block by the callers PID to
 * prevent it from clashing with concurrent allocations for a different inode
 * in the same block group.   The PID is used here so that functionally related
 * files will be close-by on-disk.
 *
 *	Caller must make sure that @ind is valid and will stay that way.
 */

static ext2_fsblk_t ext2_find_near(struct inode *inode, Indirect *ind)
{
	struct ext2_inode_info *ei = EXT2_I(inode);
	__le32 *start = ind->bh ? (__le32 *) ind->bh->b_data : ei->i_data;
	__le32 *p;
	ext2_fsblk_t bg_start;
	ext2_fsblk_t colour;

	/* Try to find previous block */
	for (p = ind->p - 1; p >= start; p--)
#ifdef CONFIG_EXT2_COMPRESS
		if (!HOLE_BLKADDR(*p))
#else
		if (*p)
#endif
			return le32_to_cpu(*p);

	/* No such thing, so let's try location of indirect block */
	if (ind->bh)
		return ind->bh->b_blocknr;

	/*
	 * It is going to be referred from inode itself? OK, just put it into
	 * the same cylinder group then.
	 */
	bg_start = ext2_group_first_block_no(inode->i_sb, ei->i_block_group);
	colour = (current->pid % 16) *
			(EXT2_BLOCKS_PER_GROUP(inode->i_sb) / 16);
	return bg_start + colour;
}

/**
 *	ext2_find_goal - find a preferred place for allocation.
 *	@inode: owner
 *	@block:  block we want
 *	@partial: pointer to the last triple within a chain
 *
 *	Returns preferred place for a block (the goal).
 */

static inline ext2_fsblk_t ext2_find_goal(struct inode *inode, long block,
					  Indirect *partial)
{
	struct ext2_block_alloc_info *block_i;

	block_i = EXT2_I(inode)->i_block_alloc_info;

	/*
	 * try the heuristic for sequential allocation,
	 * failing that at least try to get decent locality.
	 */
	if (block_i && (block == block_i->last_alloc_logical_block + 1)
		&& (block_i->last_alloc_physical_block != 0)) {
		return block_i->last_alloc_physical_block + 1;
	}

	return ext2_find_near(inode, partial);
}

/**
 *	ext2_blks_to_allocate: Look up the block map and count the number
 *	of direct blocks need to be allocated for the given branch.
 *
 * 	@branch: chain of indirect blocks
 *	@k: number of blocks need for indirect blocks
 *	@blks: number of data blocks to be mapped.
 *	@blocks_to_boundary:  the offset in the indirect block
 *
 *	return the total number of blocks to be allocate, including the
 *	direct and indirect blocks.
 */
static int
ext2_blks_to_allocate(Indirect * branch, int k, unsigned long blks,
		int blocks_to_boundary)
{
	unsigned long count = 0;

	/*
	 * Simple case, [t,d]Indirect block(s) has not allocated yet
	 * then it's clear blocks on that path have not allocated
	 */
	if (k > 0) {
		/* right now don't hanel cross boundary allocation */
		if (blks < blocks_to_boundary + 1)
			count += blks;
		else
			count += blocks_to_boundary + 1;
		return count;
	}

	count++;
	while (count < blks && count <= blocks_to_boundary
		&& le32_to_cpu(*(branch[0].p + count)) == 0) {
		count++;
	}
	return count;
}

/**
 *	ext2_alloc_blocks: multiple allocate blocks needed for a branch
 *	@indirect_blks: the number of blocks need to allocate for indirect
 *			blocks
 *
 *	@new_blocks: on return it will store the new block numbers for
 *	the indirect blocks(if needed) and the first direct block,
 *	@blks:	on return it will store the total number of allocated
 *		direct blocks
 */
static int ext2_alloc_blocks(struct inode *inode,
			ext2_fsblk_t goal, int indirect_blks, int blks,
			ext2_fsblk_t new_blocks[4], int *err)
{
	int target, i;
	unsigned long count = 0;
	int index = 0;
	ext2_fsblk_t current_block = 0;
	int ret = 0;

	/*
	 * Here we try to allocate the requested multiple blocks at once,
	 * on a best-effort basis.
	 * To build a branch, we should allocate blocks for
	 * the indirect blocks(if not allocated yet), and at least
	 * the first direct block of this branch.  That's the
	 * minimum number of blocks need to allocate(required)
	 */
	target = blks + indirect_blks;

	while (1) {
		count = target;
		/* allocating blocks for indirect blocks and direct blocks */
		current_block = ext2_new_blocks(inode,goal,&count,err);
		if (*err)
			goto failed_out;

		target -= count;
		/* allocate blocks for indirect blocks */
		while (index < indirect_blks && count) {
			new_blocks[index++] = current_block++;
			count--;
		}

		if (count > 0)
			break;
	}

	/* save the new block number for the first direct block */
	new_blocks[index] = current_block;

	/* total number of blocks allocated for direct blocks */
	ret = count;
	*err = 0;
	return ret;
failed_out:
	for (i = 0; i <index; i++)
		ext2_free_blocks(inode, new_blocks[i], 1);
	if (index)
		mark_inode_dirty(inode);
	return ret;
}

/**
 *	ext2_alloc_branch - allocate and set up a chain of blocks.
 *	@inode: owner
 *	@num: depth of the chain (number of blocks to allocate)
 *	@offsets: offsets (in the blocks) to store the pointers to next.
 *	@branch: place to store the chain in.
 *
 *	This function allocates @num blocks, zeroes out all but the last one,
 *	links them into chain and (if we are synchronous) writes them to disk.
 *	In other words, it prepares a branch that can be spliced onto the
 *	inode. It stores the information about that chain in the branch[], in
 *	the same format as ext2_get_branch() would do. We are calling it after
 *	we had read the existing part of chain and partial points to the last
 *	triple of that (one with zero ->key). Upon the exit we have the same
 *	picture as after the successful ext2_get_block(), except that in one
 *	place chain is disconnected - *branch->p is still zero (we did not
 *	set the last link), but branch->key contains the number that should
 *	be placed into *branch->p to fill that gap.
 *
 *	If allocation fails we free all blocks we've allocated (and forget
 *	their buffer_heads) and return the error value the from failed
 *	ext2_alloc_block() (normally -ENOSPC). Otherwise we set the chain
 *	as described above and return 0.
 */

static int ext2_alloc_branch(struct inode *inode,
			int indirect_blks, int *blks, ext2_fsblk_t goal,
			int *offsets, Indirect *branch)
{
	int blocksize = inode->i_sb->s_blocksize;
	int i, n = 0;
	int err = 0;
	struct buffer_head *bh;
	int num;
	ext2_fsblk_t new_blocks[4];
	ext2_fsblk_t current_block;

	num = ext2_alloc_blocks(inode, goal, indirect_blks,
				*blks, new_blocks, &err);
	if (err)
		return err;

	branch[0].key = cpu_to_le32(new_blocks[0]);
	/*
	 * metadata blocks and data blocks are allocated.
	 */
	for (n = 1; n <= indirect_blks;  n++) {
		/*
		 * Get buffer_head for parent block, zero it out
		 * and set the pointer to new one, then send
		 * parent to disk.
		 */
		bh = sb_getblk(inode->i_sb, new_blocks[n-1]);
		branch[n].bh = bh;
#ifndef CONFIG_EXT2_COMPRESS
		lock_buffer(bh);
#else
		CHECK_NOT_ATOMIC
		if (!buffer_uptodate(bh))
 			wait_on_buffer(bh);
#endif	
		memset(bh->b_data, 0, blocksize);
		branch[n].p = (__le32 *) bh->b_data + offsets[n];
		branch[n].key = cpu_to_le32(new_blocks[n]);
		*branch[n].p = branch[n].key;
		if ( n == indirect_blks) {
			current_block = new_blocks[n];
			/*
			 * End of chain, update the last new metablock of
			 * the chain to point to the new allocated
			 * data blocks numbers
			 */
			for (i=1; i < num; i++)
				*(branch[n].p + i) = cpu_to_le32(++current_block);
		}
		set_buffer_uptodate(bh);
#ifndef CONFIG_EXT2_COMPRESS
		unlock_buffer(bh);
#endif
		mark_buffer_dirty_inode(bh, inode);
		/* We used to sync bh here if IS_SYNC(inode).
		 * But we now rely upon generic_write_sync()
		 * and b_inode_buffers.  But not for directories.
		 */
		if (S_ISDIR(inode->i_mode) && IS_DIRSYNC(inode))
			sync_dirty_buffer(bh);
	}
	*blks = num;
	return err;
}

/**
 * ext2_splice_branch - splice the allocated branch onto inode.
 * @inode: owner
 * @block: (logical) number of block we are adding
 * @where: location of missing link
 * @num:   number of indirect blocks we are adding
 * @blks:  number of direct blocks we are adding
 *
 * This function fills the missing link and does all housekeeping needed in
 * inode (->i_blocks, etc.). In case of success we end up with the full
 * chain to new block and return 0.
 */
static void ext2_splice_branch(struct inode *inode,
			long block, Indirect *where, int num, int blks)
{
	int i;
	struct ext2_block_alloc_info *block_i;
	ext2_fsblk_t current_block;

	block_i = EXT2_I(inode)->i_block_alloc_info;

	/* XXX LOCKING probably should have i_meta_lock ?*/
	/* That's it */

	*where->p = where->key;

	/*
	 * Update the host buffer_head or inode to point to more just allocated
	 * direct blocks blocks
	 */
	if (num == 0 && blks > 1) {
		current_block = le32_to_cpu(where->key) + 1;
		for (i = 1; i < blks; i++)
			*(where->p + i ) = cpu_to_le32(current_block++);
	}

	/*
	 * update the most recently allocated logical & physical block
	 * in i_block_alloc_info, to assist find the proper goal block for next
	 * allocation
	 */
	if (block_i) {
		block_i->last_alloc_logical_block = block + blks - 1;
		block_i->last_alloc_physical_block =
				le32_to_cpu(where[num].key) + blks - 1;
	}

	/* We are done with atomic stuff, now do the rest of housekeeping */

	/* had we spliced it onto indirect block? */
	if (where->bh)
		mark_buffer_dirty_inode(where->bh, inode);

	inode->i_ctime = CURRENT_TIME_SEC;
	mark_inode_dirty(inode);
}

/*
 * Allocation strategy is simple: if we have to allocate something, we will
 * have to go the whole way to leaf. So let's do it before attaching anything
 * to tree, set linkage between the newborn blocks, write them if sync is
 * required, recheck the path, free and repeat if check fails, otherwise
 * set the last missing link (that will protect us from any truncate-generated
 * removals - all blocks on the path are immune now) and possibly force the
 * write on the parent block.
 * That has a nice additional property: no special recovery from the failed
 * allocations is needed - we simply release blocks and do not touch anything
 * reachable from inode.
 *
 * `handle' can be NULL if create == 0.
 *
 * return > 0, # of blocks mapped or allocated.
 * return = 0, if plain lookup failed.
 * return < 0, error case.
 */
static int ext2_get_blocks(struct inode *inode,
			   sector_t iblock, unsigned long maxblocks,
			   struct buffer_head *bh_result,
			   int create)
{
	int err = -EIO;
	int offsets[4];
	Indirect chain[4];
	Indirect *partial;
	ext2_fsblk_t goal;
	int indirect_blks;
	int blocks_to_boundary = 0;
	int depth;
	struct ext2_inode_info *ei = EXT2_I(inode);
	int count = 0;
	ext2_fsblk_t first_block = 0;

	depth = ext2_block_to_path(inode,iblock,offsets,&blocks_to_boundary);

	if (depth == 0)
		return (err);

	partial = ext2_get_branch(inode, depth, offsets, chain, &err);
	/* Simplest case - block found, no allocation needed */
	if (!partial) {
		first_block = le32_to_cpu(chain[depth - 1].key);
		clear_buffer_new(bh_result); /* What's this do? */
		count++;
		/*map more blocks*/
		while (count < maxblocks && count <= blocks_to_boundary) {
			ext2_fsblk_t blk;

			if (!verify_chain(chain, chain + depth - 1)) {
				/*
				 * Indirect block might be removed by
				 * truncate while we were reading it.
				 * Handling of that case: forget what we've
				 * got now, go to reread.
				 */
				err = -EAGAIN;
				count = 0;
				break;
			}
			blk = le32_to_cpu(*(chain[depth-1].p + count));
			if (blk == first_block + count)
				count++;
			else
				break;
		}
		if (err != -EAGAIN)
			goto got_it;
	}

	/* Next simple case - plain lookup or failed read of indirect block */
	if (!create || err == -EIO)
		goto cleanup;

	mutex_lock(&ei->truncate_mutex);
	/*
	 * If the indirect block is missing while we are reading
	 * the chain(ext2_get_branch() returns -EAGAIN err), or
	 * if the chain has been changed after we grab the semaphore,
	 * (either because another process truncated this branch, or
	 * another get_block allocated this branch) re-grab the chain to see if
	 * the request block has been allocated or not.
	 *
	 * Since we already block the truncate/other get_block
	 * at this point, we will have the current copy of the chain when we
	 * splice the branch into the tree.
	 */
	if (err == -EAGAIN || !verify_chain(chain, partial)) {
		while (partial > chain) {
			brelse(partial->bh);
//			bforget(partial->bh); /*mw: e2c-pre-2.6.30.4 used bforget here*/
			partial--;
		}
		partial = ext2_get_branch(inode, depth, offsets, chain, &err);
		if (!partial) {
			count++;
			mutex_unlock(&ei->truncate_mutex);
			if (err)
				goto cleanup;
			clear_buffer_new(bh_result);
			goto got_it;
		}
	}

	/*
	 * Okay, we need to do block allocation.  Lazily initialize the block
	 * allocation info here if necessary
	*/
	if (S_ISREG(inode->i_mode) && (!ei->i_block_alloc_info))
		ext2_init_block_alloc_info(inode);

	goal = ext2_find_goal(inode, iblock, partial);

	/* the number of blocks need to allocate for [d,t]indirect blocks */
	indirect_blks = (chain + depth) - partial - 1;
	/*
	 * Next look up the indirect map to count the totoal number of
	 * direct blocks to allocate for this branch.
	 */
	count = ext2_blks_to_allocate(partial, indirect_blks,
					maxblocks, blocks_to_boundary);
	/*
	 * XXX ???? Block out ext2_truncate while we alter the tree
	 */
	err = ext2_alloc_branch(inode, indirect_blks, &count, goal,
				offsets + (partial - chain), partial);

	if (err) {
		mutex_unlock(&ei->truncate_mutex);
		goto cleanup;
	}

	if (ext2_use_xip(inode->i_sb)) {
		/*
		 * we need to clear the block
		 */
		err = ext2_clear_xip_target (inode,
			le32_to_cpu(chain[depth-1].key));
		if (err) {
			mutex_unlock(&ei->truncate_mutex);
			goto cleanup;
		}
	}

	ext2_splice_branch(inode, iblock, partial, indirect_blks, count);
	mutex_unlock(&ei->truncate_mutex);
	set_buffer_new(bh_result);
got_it:
	map_bh(bh_result, inode->i_sb, le32_to_cpu(chain[depth-1].key));
	if (count > blocks_to_boundary)
		set_buffer_boundary(bh_result);
	err = count;
	/* Clean up and exit */
	partial = chain + depth - 1;	/* the whole chain */
cleanup:
	while (partial > chain) {
		brelse(partial->bh);
		partial--;
	}
	return err;
}

int ext2_get_block(struct inode *inode, sector_t iblock, struct buffer_head *bh_result, int create)
{
	unsigned max_blocks = bh_result->b_size >> inode->i_blkbits;
	int ret = ext2_get_blocks(inode, iblock, max_blocks,
			      bh_result, create);
	if (ret > 0) {
		bh_result->b_size = (ret << inode->i_blkbits);
		ret = 0;
	}
	return ret;

}

int ext2_fiemap(struct inode *inode, struct fiemap_extent_info *fieinfo,
		u64 start, u64 len)
{
	return generic_block_fiemap(inode, fieinfo, start, len,
				    ext2_get_block);
}

#ifdef CONFIG_EXT2_COMPRESS
/*
 *    Readpage method that will take care of decompression.
 */
/* effic: I (pjm) think tht at present, reading a 32KB cluster 4KB at
   a time does `decompress 4KB' for the first 4KB, then `decompress
   8KB' for the second, and so on.  See if we can provide the page
   cache with all the pages in a cluster.  The problem is, we don't
   want to erase anything tht hasn't been written to disk, so we can't
   just call update_vm_cache().  The plan at present is to remember
   what the contents of ext2_rd_wa.u come from, and don't bother
   decompressing anything if the working area already contains the
   right data.  However, this is only a win where adjacent calls to
   ext2_decompress_blocks() request the same cluster.  We could force
   that by copying some code from generic_file_read() (but check for
   deadlocks before doing anything like that), but instead I'm taking
   the more passive approach of hoping for the best. */
static int ext2_readpage(struct file *file, struct page *page)
{
	struct inode *inode = page->mapping->host;
	struct page *pg[EXT2_MAX_CLUSTER_PAGES], *epg[EXT2_MAX_CLUSTER_PAGES];
	u32 cluster0, max_cluster;
	int i, blockOfCluster, blocksToDo, npg;
	const int inc = PAGE_SIZE >> inode->i_sb->s_blocksize_bits;
	struct ext2_inode_info *ei = EXT2_I(page->mapping->host);
#ifdef CONFIG_HIGHMEM	
	int kmapped = 0; //mw
#endif	
	
	int iClusterCnt;

	/* For directories, fall out through default routine */
	if (S_ISDIR(inode->i_mode))
	{	
		int rc;
		
		rc = block_read_full_page(page,ext2_get_block);
		assert(!rc);
		return rc;
	}

	/* The semaphore prevents us trying to compress and decompress
	   the cluster at the same time, or compress a cluster in the
	   middle of reading it (thinking it to be uncompressed).

	   You may not like the fact that we hold the semaphore across
	   readpage (given that it isn't held without e2compr compiled
	   in), but it does guarantee that we won't compress the
	   cluster during readpage.  (OTOH, it's unlikely, if not
	   impossible, for someone to ,compress a cluster and rewrite
	   the blocks` before the readpage completes.) */
	/* This procedure used to have `#ifndef EXT2_LOCK_BUFFERS'
	   around all the semaphore stuff, and unlocked each buffer
	   before brelsing them ifdef EXT2_LOCK_BUFFERS.  I (pjm,
	   1998-01-20) have removed that because (a) EXT2_LOCK_BUFFERS
	   isn't #defined anywhere, and doesn't appear outside of this
	   function, and (b) I haven't looked at what effect locking
	   the buffers has.  You may like to reintroduce the idea of
	   buffer locking to this function if you're more familiar
	   with buffer locking than I, and believe that the full i_sem
	   isn't necessary to protect from races (people seeing raw
	   compressed data) between readpage and ext2_file_write(),
	   ext2_compress_cluster() and ext2_truncate(). */
	unlock_page(page);
	mutex_lock(&inode->i_mutex);
	
	assert (atomic_read(&inode->i_mutex.count) <= 0); /* i.e. mutex_lock */
	
	//mw: added EXT2_COMPR_FL, because EXT2_COMPRBLK_FL mit change without mutex !!!
	if ( !(ei->i_flags & (EXT2_COMPRBLK_FL|EXT2_COMPR_FL))  
		||  (ei->i_flags & EXT2_NOCOMPR_FL) )
	{
		goto readpage_uncompressed;
	}

	{
		register u32 blockOfFile
		    = (page->index << PAGE_CACHE_SHIFT) >> inode->i_sb->s_blocksize_bits;

		blocksToDo = PAGE_SIZE >> inode->i_sb->s_blocksize_bits;
		cluster0 = ext2_block_to_cluster(inode, blockOfFile);
		max_cluster = ext2_block_to_cluster
		    (inode, blockOfFile + blocksToDo - 1);
		blockOfCluster
		    = blockOfFile - ext2_cluster_block0(inode, cluster0);
	}

	/* return -???, any idea which code. do_generic_file_read() cares, ext2_readpages() doesn't. 
	   maybe I should look at the "generic" readpage() and see what it returns in this case */

	/* Check if any part of the requested area contains part of a
	   compressed cluster.  If not, we can use default ext2_readpage().

	   (Note that we don't have to worry about a cluster becoming
	   compressed in the meantime, because we have the semaphore.)

	   A page can cover up to 9 clusters.  (The maximum can only
	   occur with 32KB pages, 4KB clusters, and a non-page-aligned
	   offset.  Thanks go to Kurt Fitzner for reporting that
	   page offsets needn't be aligned; see generic_file_mmap().)  */
	{
	int isCmp[(PAGE_SIZE >> 12) + 1];
	u8 *dst;
	unsigned clu_ix;

	assert (max_cluster - cluster0 < sizeof(isCmp)/sizeof(*isCmp));
	for (clu_ix = 0; cluster0 + clu_ix <= max_cluster; clu_ix++) {
		isCmp[clu_ix] = ext2_cluster_is_compressed_fn (inode, cluster0 + clu_ix);
		if (isCmp[clu_ix] < 0){
			printk("IO-ERROR: isCmp\n");
			goto io_error;
		}
	} 
	
	for (clu_ix = 0; cluster0 + clu_ix <= max_cluster; clu_ix++)
		if (isCmp[clu_ix] > 0)
			goto readpage_compressed;
	/* fall through */
	readpage_uncompressed:
	{
		int rc=0;
		lock_page(page);

		/* Did somebody else fill it already? */
		if (PageUptodate(page) ){  //mw: necessary for DEBUG! anyway checked in do_generic_mapping_read
			unlock_page(page);
		}
		else {
			//try_to_free_buffers(page);
			rc = block_read_full_page(page,ext2_get_block);
		}
		mutex_unlock(&inode->i_mutex);
		assert(!rc);			
		return rc;
	}

	readpage_compressed:

	/* Copied from block_read_full_page */
	/* if (!PageLocked(page)) */
	/* PAGE_BUG(page); */
	lock_page(page);
	if (PageUptodate(page)) { 
		unlock_page(page);
		mutex_unlock(&inode->i_mutex);
		return(0);
	}
	get_page(page);

	ClearPageUptodate(page);
	ClearPageError(page);

	dst = (u8 *) page_address(page);
	for (clu_ix = 0; cluster0 + clu_ix <= max_cluster; clu_ix++) {
		struct buffer_head *bh[EXT2_MAX_CLUSTER_BLOCKS];
		int nbh, blocksThisClu;
		
		for (i = 0; i < EXT2_MAX_CLUSTER_PAGES; i++) {
			pg[i] = NULL;
			epg[i] = NULL;
		}
			
		/* clear_bit(PG_locked, &page->flags); */
		npg = ext2_cluster_npages(inode, cluster0 + clu_ix);
		nbh = ext2_get_cluster_pages(inode, cluster0 + clu_ix, pg, page, 0);
		
		if (nbh <= 0) {
			for (i = 0; i < EXT2_MAX_CLUSTER_PAGES; i++)
			printk("no pages\n");
			goto out;
		}		
		iClusterCnt =  ext2_cluster_npages(inode, cluster0);
		
		nbh =  ext2_get_cluster_extra_pages(inode, cluster0 + clu_ix, pg, epg);
		if (nbh <= 0) 
		{
			for (i = 0; i < EXT2_MAX_CLUSTER_PAGES; i++)
				epg[i] = NULL;
			printk("no extra pages\n");
			goto out;
		}
		assert (iClusterCnt = ext2_cluster_npages(inode, cluster0));
		
#ifdef CONFIG_HIGHMEM
		ext2_kmap_cluster_pages(page, pg, epg);
		kmapped = 1;
#endif
		
		nbh = ext2_get_cluster_blocks(inode, cluster0 + clu_ix, bh, pg, epg, 0);
		if (nbh <= 0)
		{
			printk("no blocks\n");
			goto out;
		}	
			
		/* How many blocks (including holes) we need from this cluster. */
		{
			blocksThisClu = (ext2_cluster_nblocks(inode, cluster0 +
			    clu_ix) - blockOfCluster);
			if (blocksThisClu > blocksToDo)
				blocksThisClu = blocksToDo;
		}

		if (isCmp[clu_ix]) {
			u8 const *src;
			int n, nbytes_wanted;
			struct ext2_cluster_head *head;
			unsigned meth;
# ifdef CONFIG_KMOD
			unsigned alg;
# endif

			bh[0]->b_data = page_address(bh[0]->b_page);
            		head = (struct ext2_cluster_head *) bh[0]->b_data;

			/* jmr 1998-10-28 Hope this is the last time I'm moving this code.
			 * Module loading must be done _before_ we lock wa, just think what
			 * can happen if we reallocate wa when somebody else uses it...
			 */
			meth = head->method; /* only a byte, so no swabbing needed. */
			if (meth >= EXT2_N_METHODS) {
				printk("illegal method id\n");
				ext2_msg(inode->i_sb,
				    "illegal method id",
				    "inode = %lu, id = %u",
				    inode->i_ino, meth);
				goto out;
			}
# ifdef CONFIG_KMOD
			alg = ext2_method_table[meth].alg;
			if (!ext2_algorithm_table[alg].avail) {
				char str[32];

				sprintf(str, "ext2-compr-%s", ext2_algorithm_table[alg].name);
				request_module(str);
			}
# endif /* CONFIG_KMOD */

			/* Calculate nbytes_wanted. */
			{
				unsigned nblk_wanted, i;

				/* We want to decompress the whole cluster */
				//nblk_wanted = ext2_cluster_nblocks(inode, cluster0 + clu_ix);
				nblk_wanted = npg << (PAGE_CACHE_SHIFT - inode->i_sb->s_blocksize_bits); /*mw: FIXED */
				
				for (i = nblk_wanted; i != 0;)
					if (((--i >> 3) < head->holemap_nbytes)
					    && (head->holemap[i >> 3] & (1 << (i & 7))))
						--nblk_wanted;
				nbytes_wanted = (nblk_wanted
				    << inode->i_sb->s_blocksize_bits);
			}

			/* Decompress. */
			get_cpu_var(ext2_rd_wa);
			if (__get_cpu_var(ext2_rd_wa) == NULL)
			{
				ext2_alloc_rd_wa();
			}
			assert(__get_cpu_var(ext2_rd_wa) != NULL);

			n = ext2_decompress_blocks(inode, bh, nbh, nbytes_wanted, cluster0 + clu_ix);
			if (n < 0) {
			    assert(nbh >= 0);
			    printk("ext2_readpage: noblocks decompressed\n");
			    put_cpu_var(ext2_rd_wa);	
			    goto out;
			}

# ifdef EXT2_COMPR_REPORT_VERBOSE_INODE
			if (ei->i_flags & EXT2_COMPR_FL)
				printk(KERN_DEBUG "ext2: mmap %04x:%lu: blocksToDo=%d, blockOfCluster=%d, blocksThisClu=%d, clu_nblocks=%d\n",
				    inode->i_rdev,
				    inode->i_ino,
				    blocksToDo,
				    blockOfCluster,
				    blocksThisClu,
				    ext2_cluster_nblocks(inode, cluster0 + clu_ix));
# endif

			/* */
			{
			unsigned i;
			int ipg;
			
			i = ext2_cluster_nblocks(inode, cluster0 + clu_ix) - 1;
			//i = (npg << (PAGE_CACHE_SHIFT - inode->i_sb->s_blocksize_bits)) - 1; /*mw: FIXED!!!  (here: shift = 2Bit) */
			//if(i+1 != ext2_cluster_nblocks(inode, cluster0 + clu_ix))
			//printk("npg=%i, nbh=%i, npgf=%i, nbhf =%i, cluster:%i, dec_blk:%i, b_wanted:%i, size:%i\n ", ext2_cluster_npages(inode, cluster0 + clu_ix),  ext2_cluster_nblocks(inode, cluster0 + clu_ix), npgtest, i+1, cluster0 + clu_ix, n, nbytes_wanted, inode->i_size);
			blockOfCluster = 0;
			assert(n > 0);
			src = __get_cpu_var(ext2_rd_wa)->u + nbytes_wanted - inode->i_sb->s_blocksize;
#ifdef  EXT2_COMPR_REPORT			
			trace_e2c("ext2_readpage: copy data inc=%d blocksThisClu=%d, n=%d\n", inc, blocksThisClu, n);
#endif			
			for (ipg = npg - 1; ipg >= 0; ipg--) {
				if (pg[ipg] == NULL) {
				i -= inc;
				src -= PAGE_SIZE;
				continue;
			}
			if (((inode->i_size-1) >> PAGE_SHIFT) ==  pg[ipg]->index) {
				n = ((inode->i_size-1) & (PAGE_SIZE -1)) >> inode->i_sb->s_blocksize_bits;
				i -= ((blocksThisClu-1) - n);
				src -= ((blocksThisClu-1) - n) << inode->i_sb->s_blocksize_bits;
			} else {
				n = blocksThisClu - 1;
			}
			if (PageUptodate(pg[ipg])  ) {
				for (;n >= 0;n--, i--) {
					if (((i >> 3) >= head->holemap_nbytes)
					    || !(head->holemap[i >> 3] & (1 << (i & 7)))) {
						src -= inode->i_sb->s_blocksize;
					}
				}
			} else {
                
				dst = (u8 *) page_address(pg[ipg]) + (n << inode->i_sb->s_blocksize_bits);
				
				for (;
				    n >= 0;
				    n--, i--, dst -= inode->i_sb->s_blocksize) {
					assert(!buffer_dirty(bh[i]));
					clear_buffer_dirty(bh[i]);   //mw: had a refile_buffer in 2.4 
					if (((i >> 3) >= head->holemap_nbytes)
					    || !(head->holemap[i >> 3] & (1 << (i & 7)))) {
						assert(i >= 0);
						memcpy(dst, src, inode->i_sb->s_blocksize);
						src -= inode->i_sb->s_blocksize;
					} else {
						assert(i >= 0);
						memset (dst, 0, inode->i_sb->s_blocksize);
					}
					//clear_bit(BH_Uptodate, &bh[i]->b_state);
 				}
				SetPageUptodate(pg[ipg]);
			}
			} 
			}
			put_cpu_var(ext2_rd_wa);
		} else {
			/* Uncompressed cluster.  Just copy the data.  */
			int n;

# ifdef EXT2_COMPR_REPORT_VERBOSE_INODE
			if (ei->i_flags & EXT2_COMPR_FL)
				printk(KERN_DEBUG
				    "ext2: mmap %lu: blocksToDo = %d, "
				    "blockOfCluster = %d, clu_nblocks = %d\n",
				    inode->i_ino, blocksToDo, blockOfCluster,
				    ext2_cluster_nblocks(inode, cluster0 +
				    clu_ix));
# endif

			for (n = 0;
			    n < blocksThisClu;
			    n++, dst += inode->i_sb->s_blocksize) {
				if ((blockOfCluster + n < nbh)
				    && (bh[blockOfCluster + n] != NULL))
				{
					memcpy(dst,
					    bh[blockOfCluster + n]->b_data,
					    inode->i_sb->s_blocksize);
				}
				else
				{
					memset(dst, 0, inode->i_sb->s_blocksize);
				}
			}
			blockOfCluster = 0;
		} // end uncompressed Cluster

		blocksToDo -= blocksThisClu;

#ifdef CONFIG_HIGHMEM
		if (kmapped)
			ext2_kunmap_cluster_pages(page, pg, epg);
#endif

		for (i = 0; i < EXT2_MAX_CLUSTER_PAGES; i++) {
			if (epg[i] != NULL) {
				
				ClearPageDirty(epg[i]);			
				ClearPageUptodate(epg[i]);	
				try_to_free_buffers(epg[i]);
				unlock_page(epg[i]);
				assert(page_count(epg[i]) <= 1);
				page_cache_release(epg[i]);
			}
		}

		for (i = 0; i < EXT2_MAX_CLUSTER_PAGES; i++) {
			if (pg[i] == NULL)
				break;		
			if (pg[i] == page)
				continue;
			unlock_page(pg[i]);
			page_cache_release(pg[i]);
		}		
		//mw
		assert (isCmp[clu_ix] == ext2_cluster_is_compressed_fn (inode, cluster0 + clu_ix));
	} // end for-loop: Cluster	
	} 

	SetPageUptodate(page);		
	unlock_page(page);
	atomic_dec(&page->_count);
	mutex_unlock(&inode->i_mutex);
	return 0;
	
	out:

#ifdef CONFIG_HIGHMEM
	if (kmapped)
		ext2_kunmap_cluster_pages(page, pg, epg);
#endif

	for (i = 0; i < EXT2_MAX_CLUSTER_PAGES; i++) {
		if (epg[i] != NULL) {
			
			ClearPageDirty(epg[i]);			
			ClearPageUptodate(epg[i]);		
			try_to_free_buffers(epg[i]);
			unlock_page(epg[i]);
			assert(page_count(epg[i]) <= 1);
			page_cache_release(epg[i]);
		}
	}

	for (i = 0; i < EXT2_MAX_CLUSTER_PAGES; i++) {
		if (pg[i] == NULL)
			break;
		if (pg[i] == page)
			continue;
		unlock_page(pg[i]);
		page_cache_release(pg[i]);	
	}
 	mutex_unlock(&inode->i_mutex);
 	return 0;

	io_error:
#ifdef CONFIG_HIGHMEM
	if (kmapped)
		ext2_kunmap_cluster_pages(page, pg, epg);
#endif
	SetPageError(page);
	unlock_page(page);
	atomic_dec(&page->_count);
	mutex_unlock(&inode->i_mutex);
	printk("Readpage: IOERROR\n");
	return -EIO; /* it is tested in do_generic_file_read(), ...     */
}
#endif /* CONFIG_EXT2_COMPRESS */

static int ext2_writepage(struct page *page, struct writeback_control *wbc)
{
/* mw (24/06/2008):
 *         WRITEPAGE: this code was also in e2compr 2.4 and once removed by yaboo ding.
 *         ext2_writepage() is also called for dirty pages. Usually we write using file_write() which 
 *         wraps correctly to compressed files. BUT: a writeable memory map might
 *         produce dirty pages, which will be written back normally. this should/might fail.
 *         The following code should fix this bug, but this was not tested yet. 
 */
#ifdef CONFIG_EXT2_COMPRESS
#undef USE_WRITEPAGE
//#define USE_WRITEPAGE	
#ifdef USE_WRITEPAGE	

    struct ext2_inode_info *ei = EXT2_I(page->mapping->host);
    int retval;

    struct inode *inode = page->mapping->host;
    u32 cluster0, max_cluster;                                        
    int blocksToDo;                                              
                                                                                 
    unlock_page(page);  
    //mw: do we need this ???
    //if (!(ei->i_compr_flags & EXT2_OSYNC_INODE)) {          
    /* trace_e2c("ext2_writepage: inode"); */                                          
    mutex_lock(&inode->i_mutex);           
    /* trace_e2c(" down\n"); */                                   
    //}                                                                  
    if (!(ei->i_flags & EXT2_COMPRBLK_FL)                       
    	|| (ei->i_flags & EXT2_NOCOMPR_FL) )
    {                               
    	//mw: do we need this ???
    	//if (!(ei->i_compr_flags & EXT2_OSYNC_INODE)) {           
    		/* trace_e2c("ext2_writepage: inode up 1\n"); */                                 
    		mutex_unlock(&inode->i_mutex);                                       
    	//}
    	lock_page(page);                                           
    	return block_write_full_page(page, ext2_get_block, wbc);
    }                                                                
    /* */   
    {                                                                                  
    	register u32 blockOfFile                              
        	= (page->index << PAGE_CACHE_SHIFT) >> inode->i_sb->s_blocksize_bits;
                                                                        
    	blocksToDo = PAGE_SIZE >> inode->i_sb->s_blocksize_bits;      
        cluster0 = ext2_block_to_cluster(inode, blockOfFile);          
        max_cluster = ext2_block_to_cluster(inode, blockOfFile + blocksToDo - 1);      
    }                                                             
                                                                         
    /* Check if any part of the requested area contains part of a                        
       compressed cluster.  If not, we can use default ext2_writepage(). 
       
       (Note that we don't have to worry about a cluster becoming
       compressed in the meantime, because we have the semaphore.)            
                                                                              
       A page can cover up to 9 clusters.  (The maximum can only              
       occur with 32KB pages, 4KB clusters, and a non-page-aligned   
       offset.  Thanks go to Kurt Fitzner for reporting that               
       page offsets needn't be aligned; see generic_file_mmap().)  */                  
                                                                                            
    {                                                                   
    	int isCmp[(PAGE_SIZE >> 12) + 1];                                 
        unsigned clu_ix;                                          
                                                                     
        assert (max_cluster - cluster0 < sizeof(isCmp)/sizeof(*isCmp));
        for (clu_ix = 0; cluster0 + clu_ix <= max_cluster; clu_ix++) { 
        	isCmp[clu_ix] = ext2_cluster_is_compressed_fn (inode, cluster0 + clu_ix);     
            if (isCmp[clu_ix] < 0) {                                            
            	//mw: do we need this ???if (!(ei->i_compr_flags & EXT2_OSYNC_INODE)) {  
            		/* trace_e2c("ext2_writepage: inode up 2\n"); */
            		lock_page(page);
                	mutex_unlock(&inode->i_mutex);                                            
                //}                                                               
                return -EIO;                                                    
            }                                                     
         }
        
         for (clu_ix = 0; cluster0 + clu_ix <= max_cluster; clu_ix++)   
        	if (isCmp[clu_ix] > 0)                                                     
                    ext2_decompress_cluster(inode, cluster0 + clu_ix);
         
         //mw: do we need this ???
         //if (!(ei->i_compr_flags & EXT2_OSYNC_INODE)) {                            
        	 /* trace_e2c("ext2_writepage: inode up 3\n"); */
         mutex_unlock(&inode->i_mutex);                                                   
         //}
         lock_page(page);
         
         /* fall through */                                            
    }                                                                                 
#endif /* CONFIG_EXT2_COMPRESS */
#endif       
	return block_write_full_page(page, ext2_get_block, wbc);
}

#ifndef CONFIG_EXT2_COMPRESS
static int ext2_readpage(struct file *file, struct page *page)
{
	return mpage_readpage(page, ext2_get_block);
}
#endif

static int
ext2_readpages(struct file *file, struct address_space *mapping,
		struct list_head *pages, unsigned nr_pages)
{
#ifdef CONFIG_EXT2_COMPRESS
/*
 * For now, just read each page into cache and don't worry about emitting BIOs.
 * (whitpa 02 Aug 2004).
 */

	unsigned page_idx;
	struct pagevec lru_pvec;
	int iError;
	
	pagevec_init(&lru_pvec, 0);

	for (page_idx = 0; page_idx < nr_pages; page_idx++) {
		struct page *page = list_entry(pages->prev, struct page, lru);
        
		prefetchw(&page->flags);
		list_del(&page->lru);
		
		iError = add_to_page_cache(page, mapping, page->index,	GFP_KERNEL);
		if (!iError) {
			if (!PageUptodate(page))  
			{
				(void) ext2_readpage(file, page);
			}
			else
			{
				unlock_page(page);
			}
			if (!pagevec_add(&lru_pvec, page))
				__pagevec_lru_add_file(&lru_pvec);
		} else {
			page_cache_release(page);
		}
		
	}
	pagevec_lru_add_file(&lru_pvec);
	BUG_ON(!list_empty(pages));
	return 0;
#else
	return mpage_readpages(mapping, pages, nr_pages, ext2_get_block);
#endif
}

static int
ext2_write_begin(struct file *file, struct address_space *mapping,
		loff_t pos, unsigned len, unsigned flags,
		struct page **pagep, void **fsdata)
{
	int ret;

	ret = block_write_begin(mapping, pos, len, flags, pagep,
				ext2_get_block);
	if (ret < 0)
		ext2_write_failed(mapping, pos + len);
	return ret;
}

static int ext2_write_end(struct file *file, struct address_space *mapping,
			loff_t pos, unsigned len, unsigned copied,
			struct page *page, void *fsdata)
{
	int ret;

	ret = generic_write_end(file, mapping, pos, len, copied, page, fsdata);
	if (ret < len)
		ext2_write_failed(mapping, pos + len);
	return ret;
}

static int
ext2_nobh_write_begin(struct file *file, struct address_space *mapping,
		loff_t pos, unsigned len, unsigned flags,
		struct page **pagep, void **fsdata)
{
	int ret;

	ret = nobh_write_begin(mapping, pos, len, flags, pagep, fsdata,
			       ext2_get_block);
	if (ret < 0)
		ext2_write_failed(mapping, pos + len);
	return ret;
}

static int ext2_nobh_writepage(struct page *page,
			struct writeback_control *wbc)
{
	return nobh_writepage(page, ext2_get_block, wbc);
}

#ifdef CONFIG_EXT2_COMPRESS
static sector_t ext2_do_bmap(struct address_space *mapping, sector_t block)
#else
static sector_t ext2_bmap(struct address_space *mapping, sector_t block)
#endif
{
	return generic_block_bmap(mapping,block,ext2_get_block);
}

#ifdef CONFIG_EXT2_COMPRESS
/* Return 0 instead of EXT2_COMPRESSED_BLKADDR if EXT2_NOCOMPR_FL
 * high.  This is necessary for us to be able to use
 * generic_readpage() when EXT2_NOCOMPR_FL is high.
 */
static sector_t ext2_bmap(struct address_space *mapping, sector_t block)
{
	sector_t result;
	struct inode *inode = mapping->host;

	if ((EXT2_I(inode)->i_flags & (EXT2_COMPRBLK_FL | EXT2_NOCOMPR_FL))
	    == (EXT2_COMPRBLK_FL | 0)) {
		int err;

		err = ext2_cluster_is_compressed_fn
		    (inode, ext2_block_to_cluster(inode, block));
		if (err > 0)
			ext2_msg (inode->i_sb, "ext2_bmap",
			    "compressed cluster, inode %lu",
			    inode->i_ino);
		if (err != 0)
			return 0;
	}

	result = ext2_do_bmap(mapping, block);
	if (result != EXT2_COMPRESSED_BLKADDR)
		return result;

	if (!(EXT2_SB(inode->i_sb)->s_es->s_feature_incompat
	    & cpu_to_le32(EXT2_FEATURE_INCOMPAT_COMPRESSION)))
		ext2_error(inode->i_sb, "ext2_bmap",
		    "compressed_blkaddr (ino %lu, blk %lu) "
		    "on non-compressed fs",
		    inode->i_ino, (unsigned long) block);
	if (!S_ISREG(inode->i_mode))
		ext2_error(inode->i_sb, "ext2_bmap",
		    "compressed_blkaddr for non-regular file "
		    "(ino %lu, blk %lu)",
		    inode->i_ino, (unsigned long) block);
	return 0;
}
#endif /* CONFIG_EXT2_COMPRESS */

static ssize_t
ext2_direct_IO(int rw, struct kiocb *iocb, const struct iovec *iov,
			loff_t offset, unsigned long nr_segs)
{
	struct file *file = iocb->ki_filp;
	struct address_space *mapping = file->f_mapping;
	struct inode *inode = mapping->host;
	ssize_t ret;

	ret = blockdev_direct_IO(rw, iocb, inode, iov, offset, nr_segs,
				 ext2_get_block);
	if (ret < 0 && (rw & WRITE))
		ext2_write_failed(mapping, offset + iov_length(iov, nr_segs));
	return ret;
}

static int
ext2_writepages(struct address_space *mapping, struct writeback_control *wbc)
{
#ifdef CONFIG_EXT2_COMPRESS
#ifdef USE_WRITEPAGE
	struct ext2_inode_info *ei = EXT2_I(mapping->host);
	if ( (ei->i_flags & EXT2_COMPRBLK_FL) 
			 && !(ei->i_flags & EXT2_NOCOMPR_FL))
	{
		//NULL will invoke ext2_writepage for writeback, hopefully.
		return mpage_writepages(mapping, wbc, NULL);
	}
	else
#endif	
#endif
	return mpage_writepages(mapping, wbc, ext2_get_block);
}

const struct address_space_operations ext2_aops = {
	.readpage		= ext2_readpage,
	.readpages		= ext2_readpages,
	.writepage		= ext2_writepage,
	.write_begin		= ext2_write_begin,
	.write_end		= ext2_write_end,
	.bmap			= ext2_bmap,
	.direct_IO		= ext2_direct_IO,
	.writepages		= ext2_writepages,
	.migratepage		= buffer_migrate_page,
	.is_partially_uptodate	= block_is_partially_uptodate,
	.error_remove_page	= generic_error_remove_page,
};

const struct address_space_operations ext2_aops_xip = {
	.bmap			= ext2_bmap,
	.get_xip_mem		= ext2_get_xip_mem,
};

const struct address_space_operations ext2_nobh_aops = {
	.readpage		= ext2_readpage,
	.readpages		= ext2_readpages,
	.writepage		= ext2_nobh_writepage,
	.write_begin		= ext2_nobh_write_begin,
	.write_end		= nobh_write_end,
	.bmap			= ext2_bmap,
	.direct_IO		= ext2_direct_IO,
	.writepages		= ext2_writepages,
	.migratepage		= buffer_migrate_page,
	.error_remove_page	= generic_error_remove_page,
};

/*
 * Probably it should be a library function... search for first non-zero word
 * or memcmp with zero_page, whatever is better for particular architecture.
 * Linus?
 */
static inline int all_zeroes(__le32 *p, __le32 *q)
{
	while (p < q)
		if (*p++)
			return 0;
	return 1;
}

/**
 *	ext2_find_shared - find the indirect blocks for partial truncation.
 *	@inode:	  inode in question
 *	@depth:	  depth of the affected branch
 *	@offsets: offsets of pointers in that branch (see ext2_block_to_path)
 *	@chain:	  place to store the pointers to partial indirect blocks
 *	@top:	  place to the (detached) top of branch
 *
 *	This is a helper function used by ext2_truncate().
 *
 *	When we do truncate() we may have to clean the ends of several indirect
 *	blocks but leave the blocks themselves alive. Block is partially
 *	truncated if some data below the new i_size is referred from it (and
 *	it is on the path to the first completely truncated data block, indeed).
 *	We have to free the top of that path along with everything to the right
 *	of the path. Since no allocation past the truncation point is possible
 *	until ext2_truncate() finishes, we may safely do the latter, but top
 *	of branch may require special attention - pageout below the truncation
 *	point might try to populate it.
 *
 *	We atomically detach the top of branch from the tree, store the block
 *	number of its root in *@top, pointers to buffer_heads of partially
 *	truncated blocks - in @chain[].bh and pointers to their last elements
 *	that should not be removed - in @chain[].p. Return value is the pointer
 *	to last filled element of @chain.
 *
 *	The work left to caller to do the actual freeing of subtrees:
 *		a) free the subtree starting from *@top
 *		b) free the subtrees whose roots are stored in
 *			(@chain[i].p+1 .. end of @chain[i].bh->b_data)
 *		c) free the subtrees growing from the inode past the @chain[0].p
 *			(no partially truncated stuff there).
 */

static Indirect *ext2_find_shared(struct inode *inode,
				int depth,
				int offsets[4],
				Indirect chain[4],
				__le32 *top)
{
	Indirect *partial, *p;
	int k, err;

	*top = 0;
	for (k = depth; k > 1 && !offsets[k-1]; k--)
		;
	partial = ext2_get_branch(inode, k, offsets, chain, &err);
	if (!partial)
		partial = chain + k-1;
	/*
	 * If the branch acquired continuation since we've looked at it -
	 * fine, it should all survive and (new) top doesn't belong to us.
	 */
	write_lock(&EXT2_I(inode)->i_meta_lock);
	if (!partial->key && *partial->p) {
		write_unlock(&EXT2_I(inode)->i_meta_lock);
		goto no_top;
	}
	for (p=partial; p>chain && all_zeroes((__le32*)p->bh->b_data,p->p); p--)
		;
	/*
	 * OK, we've found the last block that must survive. The rest of our
	 * branch should be detached before unlocking. However, if that rest
	 * of branch is all ours and does not grow immediately from the inode
	 * it's easier to cheat and just decrement partial->p.
	 */
	if (p == chain + k - 1 && p > chain) {
		p->p--;
	} else {
		*top = *p->p;
		*p->p = 0;
	}
	write_unlock(&EXT2_I(inode)->i_meta_lock);

	while(partial > p)
	{
		brelse(partial->bh);
		partial--;
	}
no_top:
	return partial;
}

/**
 *	ext2_free_data - free a list of data blocks
 *	@inode:	inode we are dealing with
 *	@p:	array of block numbers
 *	@q:	points immediately past the end of array
 *
 *	We are freeing all blocks referred from that array (numbers are
 *	stored as little-endian 32-bit) and updating @inode->i_blocks
 *	appropriately.
 */
static inline void ext2_free_data(struct inode *inode, __le32 *p, __le32 *q)
{
	unsigned long block_to_free = 0, count = 0;
	unsigned long nr;

	for ( ; p < q ; p++) {
		nr = le32_to_cpu(*p);
#ifdef CONFIG_EXT2_COMPRESS
		if (nr == EXT2_COMPRESSED_BLKADDR) {
			*p = 0;
			continue;
		}
#endif
		if (nr) {
			*p = 0;
			/* accumulate blocks to free if they're contiguous */
			if (count == 0)
				goto free_this;
			else if (block_to_free == nr - count)
				count++;
			else {
				ext2_free_blocks (inode, block_to_free, count);
				mark_inode_dirty(inode);
			free_this:
				block_to_free = nr;
				count = 1;
			}
		}
	}
	if (count > 0) {
		ext2_free_blocks (inode, block_to_free, count);
		mark_inode_dirty(inode);
	}
}

/**
 *	ext2_free_branches - free an array of branches
 *	@inode:	inode we are dealing with
 *	@p:	array of block numbers
 *	@q:	pointer immediately past the end of array
 *	@depth:	depth of the branches to free
 *
 *	We are freeing all blocks referred from these branches (numbers are
 *	stored as little-endian 32-bit) and updating @inode->i_blocks
 *	appropriately.
 */
static void ext2_free_branches(struct inode *inode, __le32 *p, __le32 *q, int depth)
{
	struct buffer_head * bh;
	unsigned long nr;

	if (depth--) {
		int addr_per_block = EXT2_ADDR_PER_BLOCK(inode->i_sb);
		for ( ; p < q ; p++) {
			nr = le32_to_cpu(*p);
			if (!nr)
				continue;
#ifdef CONFIG_EXT2_COMPRESS
			if (nr == EXT2_COMPRESSED_BLKADDR) {
			  *p = 0;
			  continue;
			}
#endif
			*p = 0;
			bh = sb_bread(inode->i_sb, nr);
			/*
			 * A read failure? Report error and clear slot
			 * (should be rare).
			 */ 
			if (!bh) {
				ext2_error(inode->i_sb, "ext2_free_branches",
					"Read failure, inode=%ld, block=%ld",
					inode->i_ino, nr);
				continue;
			}
			ext2_free_branches(inode,
					   (__le32*)bh->b_data,
					   (__le32*)bh->b_data + addr_per_block,
					   depth);
			bforget(bh);
			ext2_free_blocks(inode, nr, 1);
			mark_inode_dirty(inode);
		}
	} else
		ext2_free_data(inode, p, q);
}

/* pjm 1998-01-14: As far as I can tell, "I don't do any locking" is
  no longer correct, as i_sem is downed for all write() and
  truncate() stuff except where it doesn't matter (e.g. new inode). */

#ifdef CONFIG_EXT2_COMPRESS
/* If the EXT2_ECOMPR_FL bit is high, then things can go rather badly.
  This can only happen if access permission was obtained before the
  flag was raised.  Also, it shouldn't be too much of a problem
  unless the end point of truncation is a compressed cluster with a
  compression error. */

 /* From what I (Antoine) understand, the complexity of the truncate
    code is due to the fact that we don't want to free blocks that
    are still referenced.  It does not ensure that concurrent read
    operation will terminate properly, i.e., the semantic of reading
    while somebody truncates is undefined (you can either get the old
    data if you got the blocks before, or get plenty of zeros
    otherwise). */

/* todo: Provide error trapping in readiness for when i_op->truncate
  allows a return code. */
static void fix_compression (struct inode * inode)
{
	struct ext2_inode_info *ei = EXT2_I(inode);		
	/*if (atomic_read(&inode->i_mutex.count) > 0)
	{
		printk("Assert Mutex failed for file: %s \n", inode_name(inode, 0));
		dump_stack();
	}*/
		
	assert (ei->i_flags & EXT2_COMPRBLK_FL);  /* one or more compressed clusters */
	assert ((atomic_read(&inode->i_mutex.count) < 1)
		|| ((inode->i_nlink == 0)
		    && (atomic_read(&inode->i_count) == 0)));
	/* pjm 1998-01-14: I think the below comment can safely be removed, as
	   it's impossible for someone to be compressing during truncate(), because
	   i_sem is down. */
	/*   Dans le cas ou les clusters peuvent etre compresses, cela pose
	     un probleme : il faudrait stopper aussi si le cluster est
	     comprime et ne contient pas plus de donnees que i_size ne
	     permet. Sinon, on peut passer son temps a decompresser un
	     cluster que quelqu'un d'autre compresse en meme
	     temps... (TODO).  Cela ne peut arriver que si on reverifie apres
	     coup si le cluster est non compresse (ce qu'on fait a l'heure
	     actuelle) => faire autrement.

	     pjm fixme tr

	     If the clusters can be compressed, we'd have a problem: we'd
	     also need to stop if the cluster is compressed and doesn't
	     contain more data than i_size permits.  Otherwise we can spend
	     time decompressing a cluster that someone else is compressing
	     at the same time.  (TODO.)  This can only happen if we reverify
	     "apres coup" ("after the event"? "after each time"?) "si" ("if"
	     or "that") the cluster is not compressed (as we are currently
	     doing) => do differently. */

	/* todo: Handle errors from ext2_cluster_is_compressed().
	   (Except ext2_truncate() currently silently ignores errors
	   anyway.) */

	if (!ext2_offset_is_clu_boundary(inode, inode->i_size)
	    && (! ( ei->i_flags & EXT2_NOCOMPR_FL))
	    && (ext2_cluster_is_compressed_fn
		  (inode, ext2_offset_to_cluster (inode, inode->i_size))
		> 0)) {
		trace_e2c("fix_compression: inode:%lu decompress_cluster!\n", inode->i_ino);
		ext2_decompress_cluster(inode, ext2_offset_to_cluster(inode, inode->i_size));
		/* todo: Check the return code of
		   ext2_decompress_cluster().  (Then again, I don't
		   know how to report an error anyway.
		   ext2_truncate() silently ignores errors.) */
	  
		/* Organise for the cluster to be recompressed later. */
		assert (ei->i_flags & EXT2_COMPR_FL);
		
		ei->i_flags |= EXT2_DIRTY_FL;
		ei->i_compr_flags |= EXT2_CLEANUP_FL;
		mark_inode_dirty(inode);
	} else
		/* If there are no more compressed clusters, then
		   remove the EXT2_COMPRBLK_FL.  Not essential from a
		   safety point of view, but friendlier.  We only do
		   this in the `else' because the cleanup function
		   will handle it in the `if' case. */
		ext2_update_comprblk(inode);
}
#endif


static void __ext2_truncate_blocks(struct inode *inode, loff_t offset)
{
	__le32 *i_data = EXT2_I(inode)->i_data;
	struct ext2_inode_info *ei = EXT2_I(inode);
	int addr_per_block = EXT2_ADDR_PER_BLOCK(inode->i_sb);
	int offsets[4];
	Indirect chain[4];
	Indirect *partial;
	__le32 nr = 0;
	int n;
	long iblock;
	unsigned blocksize;

#ifdef CONFIG_EXT2_COMPRESS
	/* If the new size is in the middle of a compressed cluster,
	   then we decompress it, and set things up to be recompressed
	   later.

	   todo: It isn't very nice to get ENOSPC on truncate.  We
	   can't completely remove the possibility (unless the
	   compression algorithms obey the rule `shorter input never
	   gives longer output') but we could greatly reduce the
	   possibility, e.g. by moving the fix_compression() function
	   to compress.c, and have it decompress and immediately
	   recompress the cluster, without allocating blocks for the
	   full decompressed data. */
	if (EXT2_I(inode)->i_flags & EXT2_COMPRBLK_FL) {
	  	trace_e2c("ext2_truncate: ino=%ld sz=%d\n", inode->i_ino, (int)inode->i_size);
	 	fix_compression(inode);
		truncate_inode_pages(inode->i_mapping, inode->i_size);	    
	}
#endif

	blocksize = inode->i_sb->s_blocksize;
	iblock = (offset + blocksize-1) >> EXT2_BLOCK_SIZE_BITS(inode->i_sb);

	n = ext2_block_to_path(inode, iblock, offsets, NULL);
	if (n == 0)
		return;

	/*
	 * From here we block out all ext2_get_block() callers who want to
	 * modify the block allocation tree.
	 */
	mutex_lock(&ei->truncate_mutex);

	if (n == 1) {
		ext2_free_data(inode, i_data+offsets[0],
					i_data + EXT2_NDIR_BLOCKS);
		goto do_indirects;
	}

	partial = ext2_find_shared(inode, n, offsets, chain, &nr);
	/* Kill the top of shared branch (already detached) */
	if (nr) {
		if (partial == chain)
			mark_inode_dirty(inode);
		else
			mark_buffer_dirty_inode(partial->bh, inode);
		ext2_free_branches(inode, &nr, &nr+1, (chain+n-1) - partial);
	}
	/* Clear the ends of indirect blocks on the shared branch */
	while (partial > chain) {
		ext2_free_branches(inode,
				   partial->p + 1,
				   (__le32*)partial->bh->b_data+addr_per_block,
				   (chain+n-1) - partial);
		mark_buffer_dirty_inode(partial->bh, inode);
		brelse (partial->bh);
		partial--;
	}
do_indirects:
	/* Kill the remaining (whole) subtrees */
	switch (offsets[0]) {
		default:
			nr = i_data[EXT2_IND_BLOCK];
			if (nr) {
				i_data[EXT2_IND_BLOCK] = 0;
				mark_inode_dirty(inode);
				ext2_free_branches(inode, &nr, &nr+1, 1);
			}
		case EXT2_IND_BLOCK:
			nr = i_data[EXT2_DIND_BLOCK];
			if (nr) {
				i_data[EXT2_DIND_BLOCK] = 0;
				mark_inode_dirty(inode);
				ext2_free_branches(inode, &nr, &nr+1, 2);
			}
		case EXT2_DIND_BLOCK:
			nr = i_data[EXT2_TIND_BLOCK];
			if (nr) {
				i_data[EXT2_TIND_BLOCK] = 0;
				mark_inode_dirty(inode);
				ext2_free_branches(inode, &nr, &nr+1, 3);
			}
		case EXT2_TIND_BLOCK:
			;
	}

	ext2_discard_reservation(inode);

	mutex_unlock(&ei->truncate_mutex);
}
#ifdef CONFIG_EXT2_COMPRESS
void ext2_truncate_blocks(struct inode *inode, loff_t offset)
#else
static void ext2_truncate_blocks(struct inode *inode, loff_t offset)
#endif
{
	/*
	 * XXX: it seems like a bug here that we don't allow
	 * IS_APPEND inode to have blocks-past-i_size trimmed off.
	 * review and fix this.
	 *
	 * Also would be nice to be able to handle IO errors and such,
	 * but that's probably too much to ask.
	 */
	if (!(S_ISREG(inode->i_mode) || S_ISDIR(inode->i_mode) ||
	    S_ISLNK(inode->i_mode)))
		return;
	if (ext2_inode_is_fast_symlink(inode))
		return;
	if (IS_APPEND(inode) || IS_IMMUTABLE(inode))
		return;
	__ext2_truncate_blocks(inode, offset);
}

static int ext2_setsize(struct inode *inode, loff_t newsize)
{
	int error;

	if (!(S_ISREG(inode->i_mode) || S_ISDIR(inode->i_mode) ||
	    S_ISLNK(inode->i_mode)))
		return -EINVAL;
	if (ext2_inode_is_fast_symlink(inode))
		return -EINVAL;
	if (IS_APPEND(inode) || IS_IMMUTABLE(inode))
		return -EPERM;

	inode_dio_wait(inode);

	if (mapping_is_xip(inode->i_mapping))
		error = xip_truncate_page(inode->i_mapping, newsize);
	else if (test_opt(inode->i_sb, NOBH))
		error = nobh_truncate_page(inode->i_mapping,
				newsize, ext2_get_block);
	else
		error = block_truncate_page(inode->i_mapping,
				newsize, ext2_get_block);
	if (error)
		return error;

	truncate_setsize(inode, newsize);
	__ext2_truncate_blocks(inode, newsize);

	inode->i_mtime = inode->i_ctime = CURRENT_TIME_SEC;
	if (inode_needs_sync(inode)) {
		sync_mapping_buffers(inode->i_mapping);
		sync_inode_metadata(inode, 1);
	} else {
		mark_inode_dirty(inode);
	}

	return 0;
}

static struct ext2_inode *ext2_get_inode(struct super_block *sb, ino_t ino,
					struct buffer_head **p)
{
	struct buffer_head * bh;
	unsigned long block_group;
	unsigned long block;
	unsigned long offset;
	struct ext2_group_desc * gdp;

	*p = NULL;
	if ((ino != EXT2_ROOT_INO && ino < EXT2_FIRST_INO(sb)) ||
	    ino > le32_to_cpu(EXT2_SB(sb)->s_es->s_inodes_count))
		goto Einval;

	block_group = (ino - 1) / EXT2_INODES_PER_GROUP(sb);
	gdp = ext2_get_group_desc(sb, block_group, NULL);
	if (!gdp)
		goto Egdp;
	/*
	 * Figure out the offset within the block group inode table
	 */
	offset = ((ino - 1) % EXT2_INODES_PER_GROUP(sb)) * EXT2_INODE_SIZE(sb);
	block = le32_to_cpu(gdp->bg_inode_table) +
		(offset >> EXT2_BLOCK_SIZE_BITS(sb));
	if (!(bh = sb_bread(sb, block)))
		goto Eio;

	*p = bh;
	offset &= (EXT2_BLOCK_SIZE(sb) - 1);
	return (struct ext2_inode *) (bh->b_data + offset);

Einval:
	ext2_error(sb, "ext2_get_inode", "bad inode number: %lu",
		   (unsigned long) ino);
	return ERR_PTR(-EINVAL);
Eio:
	ext2_error(sb, "ext2_get_inode",
		   "unable to read inode block - inode=%lu, block=%lu",
		   (unsigned long) ino, block);
Egdp:
	return ERR_PTR(-EIO);
}

void ext2_set_inode_flags(struct inode *inode)
{
	unsigned int flags = EXT2_I(inode)->i_flags;

	inode->i_flags &= ~(S_SYNC|S_APPEND|S_IMMUTABLE|S_NOATIME|S_DIRSYNC);
	if (flags & EXT2_SYNC_FL)
		inode->i_flags |= S_SYNC;
	if (flags & EXT2_APPEND_FL)
		inode->i_flags |= S_APPEND;
	if (flags & EXT2_IMMUTABLE_FL)
		inode->i_flags |= S_IMMUTABLE;
	if (flags & EXT2_NOATIME_FL)
		inode->i_flags |= S_NOATIME;
	if (flags & EXT2_DIRSYNC_FL)
		inode->i_flags |= S_DIRSYNC;
}

/* Propagate flags from i_flags to EXT2_I(inode)->i_flags */
void ext2_get_inode_flags(struct ext2_inode_info *ei)
{
	unsigned int flags = ei->vfs_inode.i_flags;

	ei->i_flags &= ~(EXT2_SYNC_FL|EXT2_APPEND_FL|
			EXT2_IMMUTABLE_FL|EXT2_NOATIME_FL|EXT2_DIRSYNC_FL);
	if (flags & S_SYNC)
		ei->i_flags |= EXT2_SYNC_FL;
	if (flags & S_APPEND)
		ei->i_flags |= EXT2_APPEND_FL;
	if (flags & S_IMMUTABLE)
		ei->i_flags |= EXT2_IMMUTABLE_FL;
	if (flags & S_NOATIME)
		ei->i_flags |= EXT2_NOATIME_FL;
	if (flags & S_DIRSYNC)
		ei->i_flags |= EXT2_DIRSYNC_FL;
}

struct inode *ext2_iget (struct super_block *sb, unsigned long ino)
{
	struct ext2_inode_info *ei;
	struct buffer_head * bh;
	struct ext2_inode *raw_inode;
	struct inode *inode;
	long ret = -EIO;
	int n;

	inode = iget_locked(sb, ino);
	if (!inode)
		return ERR_PTR(-ENOMEM);
	if (!(inode->i_state & I_NEW))
		return inode;

	ei = EXT2_I(inode);
	ei->i_block_alloc_info = NULL;

	raw_inode = ext2_get_inode(inode->i_sb, ino, &bh);
	if (IS_ERR(raw_inode)) {
		ret = PTR_ERR(raw_inode);
 		goto bad_inode;
	}

	inode->i_mode = le16_to_cpu(raw_inode->i_mode);
	inode->i_uid = (uid_t)le16_to_cpu(raw_inode->i_uid_low);
	inode->i_gid = (gid_t)le16_to_cpu(raw_inode->i_gid_low);
	if (!(test_opt (inode->i_sb, NO_UID32))) {
		inode->i_uid |= le16_to_cpu(raw_inode->i_uid_high) << 16;
		inode->i_gid |= le16_to_cpu(raw_inode->i_gid_high) << 16;
	}
	set_nlink(inode, le16_to_cpu(raw_inode->i_links_count));
	inode->i_size = le32_to_cpu(raw_inode->i_size);
	inode->i_atime.tv_sec = (signed)le32_to_cpu(raw_inode->i_atime);
	inode->i_ctime.tv_sec = (signed)le32_to_cpu(raw_inode->i_ctime);
	inode->i_mtime.tv_sec = (signed)le32_to_cpu(raw_inode->i_mtime);
	inode->i_atime.tv_nsec = inode->i_mtime.tv_nsec = inode->i_ctime.tv_nsec = 0;
	ei->i_dtime = le32_to_cpu(raw_inode->i_dtime);
	/* We now have enough fields to check if the inode was active or not.
	 * This is needed because nfsd might try to access dead inodes
	 * the test is that same one that e2fsck uses
	 * NeilBrown 1999oct15
	 */
	if (inode->i_nlink == 0 && (inode->i_mode == 0 || ei->i_dtime)) {
		/* this inode is deleted */
		brelse (bh);
		ret = -ESTALE;
		goto bad_inode;
	}
	inode->i_blocks = le32_to_cpu(raw_inode->i_blocks);
#ifdef CONFIG_EXT2_COMPRESS
	ei->i_flags = 0x807fffff & le32_to_cpu(raw_inode->i_flags);
	ei->i_compr_flags = 0;
	if (S_ISREG(inode->i_mode) || S_ISDIR(inode->i_mode)) {
		
		if (S_ISDIR(inode->i_mode)) 
		{
			//mw:
			//mutex_lock(&inode->i_mutex);
			if (S_ISDIR(inode->i_mode))
			{
				ei->i_flags &= ~(EXT2_COMPRBLK_FL | EXT2_DIRTY_FL);  //modify!!!
			}
			//mutex_unlock(&inode->i_mutex);
		}
		
		/* The above shouldn't be necessary unless someone's
		 * been playing with EXT2_IOC_SETFLAGS on a non-e2compr
		 * kernel, or the inode has been scribbled on.
		 */
		if (ei->i_flags & (EXT2_COMPR_FL | EXT2_COMPRBLK_FL)) {
			ei->i_compr_method
			    = (le32_to_cpu(raw_inode->i_flags) >> 26) & 0x1f;
			ei->i_log2_clu_nblocks
			    = (le32_to_cpu(raw_inode->i_flags) >> 23) & 0x7;
			if ((ei->i_log2_clu_nblocks < 2)
			    || (ei->i_log2_clu_nblocks > 5)) {
				if ((ei->i_log2_clu_nblocks == 0)
				    && !(ei->i_flags & EXT2_COMPRBLK_FL)) {
					/* The EXT2_COMPR_FL flag was
					 * raised under a kernel
					 * without e2compr support.
					 */
					if (S_ISREG(inode->i_mode))
						ei->i_flags |= EXT2_DIRTY_FL;
					/* Todo: once we're sure the kernel can
					 * handle [log2_]clu_nblocks==0, get rid
					 * of the next statement.
					 */
					ei->i_log2_clu_nblocks
					    = EXT2_DEFAULT_LOG2_CLU_NBLOCKS;
				} else {
					ei->i_flags |= EXT2_ECOMPR_FL;
					ext2_error(inode->i_sb,
					    "ext2_read_inode",
					    "inode %lu is corrupted: "
					    "log2_clu_nblocks=%u",
					    inode->i_ino,
					    ei->i_log2_clu_nblocks);
				}
			}
		} else {
			ei->i_compr_method = EXT2_DEFAULT_COMPR_METHOD;
			ei->i_log2_clu_nblocks
			    = EXT2_DEFAULT_LOG2_CLU_NBLOCKS;
		}
	if (ei->i_log2_clu_nblocks >
	    (EXT2_LOG2_MAX_CLUSTER_BYTES - inode->i_sb->s_blocksize_bits))
		ei->i_log2_clu_nblocks = (EXT2_LOG2_MAX_CLUSTER_BYTES
		    - inode->i_sb->s_blocksize_bits);
	ei->i_clu_nblocks = 1 << ei->i_log2_clu_nblocks;
	if (ei->i_flags & EXT2_DIRTY_FL)
		ei->i_compr_flags = EXT2_CLEANUP_FL;
	}
#else /* !CONFIG_EXT2_COMPRESS */
	ei->i_flags = le32_to_cpu(raw_inode->i_flags);
#endif
	ei->i_faddr = le32_to_cpu(raw_inode->i_faddr);
	ei->i_frag_no = raw_inode->i_frag;
	ei->i_frag_size = raw_inode->i_fsize;
	ei->i_file_acl = le32_to_cpu(raw_inode->i_file_acl);
	ei->i_dir_acl = 0;
	if (S_ISREG(inode->i_mode))
		inode->i_size |= ((__u64)le32_to_cpu(raw_inode->i_size_high)) << 32;
	else
		ei->i_dir_acl = le32_to_cpu(raw_inode->i_dir_acl);
	ei->i_dtime = 0;
	inode->i_generation = le32_to_cpu(raw_inode->i_generation);
	ei->i_state = 0;
	ei->i_block_group = (ino - 1) / EXT2_INODES_PER_GROUP(inode->i_sb);
	ei->i_dir_start_lookup = 0;

	/*
	 * NOTE! The in-memory inode i_data array is in little-endian order
	 * even on big-endian machines: we do NOT byteswap the block numbers!
	 */
	for (n = 0; n < EXT2_N_BLOCKS; n++)
		ei->i_data[n] = raw_inode->i_block[n];

	if (S_ISREG(inode->i_mode)) {
		inode->i_op = &ext2_file_inode_operations;
		if (ext2_use_xip(inode->i_sb)) {
			inode->i_mapping->a_ops = &ext2_aops_xip;
			inode->i_fop = &ext2_xip_file_operations;
		} else if (test_opt(inode->i_sb, NOBH)) {
			inode->i_mapping->a_ops = &ext2_nobh_aops;
			inode->i_fop = &ext2_file_operations;
		} else {
			inode->i_mapping->a_ops = &ext2_aops;
			inode->i_fop = &ext2_file_operations;
		}
	} else if (S_ISDIR(inode->i_mode)) {
		inode->i_op = &ext2_dir_inode_operations;
		inode->i_fop = &ext2_dir_operations;
		if (test_opt(inode->i_sb, NOBH))
			inode->i_mapping->a_ops = &ext2_nobh_aops;
		else
			inode->i_mapping->a_ops = &ext2_aops;
	} else if (S_ISLNK(inode->i_mode)) {
		if (ext2_inode_is_fast_symlink(inode)) {
			inode->i_op = &ext2_fast_symlink_inode_operations;
			nd_terminate_link(ei->i_data, inode->i_size,
				sizeof(ei->i_data) - 1);
		} else {
			inode->i_op = &ext2_symlink_inode_operations;
			if (test_opt(inode->i_sb, NOBH))
				inode->i_mapping->a_ops = &ext2_nobh_aops;
			else
				inode->i_mapping->a_ops = &ext2_aops;
		}
	} else {
		inode->i_op = &ext2_special_inode_operations;
		if (raw_inode->i_block[0])
			init_special_inode(inode, inode->i_mode,
			   old_decode_dev(le32_to_cpu(raw_inode->i_block[0])));
		else 
			init_special_inode(inode, inode->i_mode,
			   new_decode_dev(le32_to_cpu(raw_inode->i_block[1])));
	}
	brelse (bh);
	ext2_set_inode_flags(inode);
	unlock_new_inode(inode);
	return inode;
	
bad_inode:
	iget_failed(inode);
	return ERR_PTR(ret);
}

static int __ext2_write_inode(struct inode *inode, int do_sync)
{
	struct ext2_inode_info *ei = EXT2_I(inode);
	struct super_block *sb = inode->i_sb;
	ino_t ino = inode->i_ino;
	uid_t uid = inode->i_uid;
	gid_t gid = inode->i_gid;
	struct buffer_head * bh;
	struct ext2_inode * raw_inode = ext2_get_inode(sb, ino, &bh);
	int n;
	int err = 0;

	if (IS_ERR(raw_inode))
 		return -EIO;

	/* For fields not not tracking in the in-memory inode,
	 * initialise them to zero for new inodes. */
	if (ei->i_state & EXT2_STATE_NEW)
		memset(raw_inode, 0, EXT2_SB(sb)->s_inode_size);

	ext2_get_inode_flags(ei);
	raw_inode->i_mode = cpu_to_le16(inode->i_mode);
	if (!(test_opt(sb, NO_UID32))) {
		raw_inode->i_uid_low = cpu_to_le16(low_16_bits(uid));
		raw_inode->i_gid_low = cpu_to_le16(low_16_bits(gid));
/*
 * Fix up interoperability with old kernels. Otherwise, old inodes get
 * re-used with the upper 16 bits of the uid/gid intact
 */
		if (!ei->i_dtime) {
			raw_inode->i_uid_high = cpu_to_le16(high_16_bits(uid));
			raw_inode->i_gid_high = cpu_to_le16(high_16_bits(gid));
		} else {
			raw_inode->i_uid_high = 0;
			raw_inode->i_gid_high = 0;
		}
	} else {
		raw_inode->i_uid_low = cpu_to_le16(fs_high2lowuid(uid));
		raw_inode->i_gid_low = cpu_to_le16(fs_high2lowgid(gid));
		raw_inode->i_uid_high = 0;
		raw_inode->i_gid_high = 0;
	}
	raw_inode->i_links_count = cpu_to_le16(inode->i_nlink);
	raw_inode->i_size = cpu_to_le32(inode->i_size);
	raw_inode->i_atime = cpu_to_le32(inode->i_atime.tv_sec);
	raw_inode->i_ctime = cpu_to_le32(inode->i_ctime.tv_sec);
	raw_inode->i_mtime = cpu_to_le32(inode->i_mtime.tv_sec);

	raw_inode->i_blocks = cpu_to_le32(inode->i_blocks);
	raw_inode->i_dtime = cpu_to_le32(ei->i_dtime);
#ifdef CONFIG_EXT2_COMPRESS
	if ((S_ISREG(inode->i_mode) || S_ISDIR(inode->i_mode))
	    && (ei->i_flags & (EXT2_COMPR_FL | EXT2_COMPRBLK_FL))) {
		if ((ei->i_log2_clu_nblocks < 2)
		    || (ei->i_log2_clu_nblocks > 5)) {
			ei->i_flags |= EXT2_ECOMPR_FL;
			ext2_error (inode->i_sb, "ext2_write_inode",
			    "inode %lu is corrupted: log2_clu_nblocks=%u",
			    inode->i_ino, ei->i_log2_clu_nblocks);
		}
		assert (ei->i_clu_nblocks == (1 << ei->i_log2_clu_nblocks));
		assert (ei->i_compr_method < 0x20);
		raw_inode->i_flags = cpu_to_le32
		    ((ei->i_flags & 0x807fffff)
		    | (ei->i_compr_method << 26)
		    | (ei->i_log2_clu_nblocks << 23));
	} else
	{
		//mw: i_mutex was introduced and disabled again: deadlock with lilo
		//	mutex_lock(&inode->i_mutex); //mw
		raw_inode->i_flags = cpu_to_le32	//modify !!!
		   (ei->i_flags
		   & 0x807fffff /* no compr meth/size */
		   & ~(EXT2_COMPR_FL | EXT2_COMPRBLK_FL | EXT2_IMMUTABLE_FL | EXT2_ECOMPR_FL | EXT2_NOCOMPR_FL));
		//	mutex_unlock(&inode->i_mutex); //mw
	}
#else
	raw_inode->i_flags = cpu_to_le32(ei->i_flags);
#endif
	raw_inode->i_faddr = cpu_to_le32(ei->i_faddr);
	raw_inode->i_frag = ei->i_frag_no;
	raw_inode->i_fsize = ei->i_frag_size;
	raw_inode->i_file_acl = cpu_to_le32(ei->i_file_acl);
	if (!S_ISREG(inode->i_mode))
		raw_inode->i_dir_acl = cpu_to_le32(ei->i_dir_acl);
	else {
		raw_inode->i_size_high = cpu_to_le32(inode->i_size >> 32);
		if (inode->i_size > 0x7fffffffULL) {
			if (!EXT2_HAS_RO_COMPAT_FEATURE(sb,
					EXT2_FEATURE_RO_COMPAT_LARGE_FILE) ||
			    EXT2_SB(sb)->s_es->s_rev_level ==
					cpu_to_le32(EXT2_GOOD_OLD_REV)) {
			       /* If this is the first large file
				* created, add a flag to the superblock.
				*/
				spin_lock(&EXT2_SB(sb)->s_lock);
				ext2_update_dynamic_rev(sb);
				EXT2_SET_RO_COMPAT_FEATURE(sb,
					EXT2_FEATURE_RO_COMPAT_LARGE_FILE);
				spin_unlock(&EXT2_SB(sb)->s_lock);
				ext2_write_super(sb);
			}
		}
	}
	
	raw_inode->i_generation = cpu_to_le32(inode->i_generation);
	if (S_ISCHR(inode->i_mode) || S_ISBLK(inode->i_mode)) {
		if (old_valid_dev(inode->i_rdev)) {
			raw_inode->i_block[0] =
				cpu_to_le32(old_encode_dev(inode->i_rdev));
			raw_inode->i_block[1] = 0;
		} else {
			raw_inode->i_block[0] = 0;
			raw_inode->i_block[1] =
				cpu_to_le32(new_encode_dev(inode->i_rdev));
			raw_inode->i_block[2] = 0;
		}
	} else for (n = 0; n < EXT2_N_BLOCKS; n++)
		raw_inode->i_block[n] = ei->i_data[n];
	mark_buffer_dirty(bh);
	if (do_sync) {
		sync_dirty_buffer(bh);
		if (buffer_req(bh) && !buffer_uptodate(bh)) {
			printk ("IO error syncing ext2 inode [%s:%08lx]\n",
				sb->s_id, (unsigned long) ino);
			err = -EIO;
		}
	}
	ei->i_state &= ~EXT2_STATE_NEW;
	brelse (bh);
	return err;
}

int ext2_write_inode(struct inode *inode, struct writeback_control *wbc)
{
	return __ext2_write_inode(inode, wbc->sync_mode == WB_SYNC_ALL);
}

int ext2_setattr(struct dentry *dentry, struct iattr *iattr)
{
	struct inode *inode = dentry->d_inode;
	int error;

	error = inode_change_ok(inode, iattr);
	if (error)
		return error;

	if (is_quota_modification(inode, iattr))
		dquot_initialize(inode);
	if ((iattr->ia_valid & ATTR_UID && iattr->ia_uid != inode->i_uid) ||
	    (iattr->ia_valid & ATTR_GID && iattr->ia_gid != inode->i_gid)) {
		error = dquot_transfer(inode, iattr);
		if (error)
			return error;
	}
	if (iattr->ia_valid & ATTR_SIZE && iattr->ia_size != inode->i_size) {
		error = ext2_setsize(inode, iattr->ia_size);
		if (error)
			return error;
	}
	setattr_copy(inode, iattr);
	if (iattr->ia_valid & ATTR_MODE)
		error = ext2_acl_chmod(inode);
	mark_inode_dirty(inode);

	return error;
}
