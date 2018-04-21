#
# Makefile for the linux ext2-filesystem routines.
#

ifeq ($(CONFIG_EXT2_COMPRESS),y)

COMPRESS_STUFF := adler32.o compress.o e2zlib.o\
		  $($(obj-y):%/=%/ext2-compr-%.o)
endif

obj-$(CONFIG_EXT2_FS) += ext2.o

ext2-y := balloc.o dir.o file.o ialloc.o inode.o \
 	  ioctl.o namei.o super.o symlink.o $(COMPRESS_STUFF)
 

ext2-$(CONFIG_EXT2_FS_XATTR)	 += xattr.o xattr_user.o xattr_trusted.o
ext2-$(CONFIG_EXT2_FS_POSIX_ACL) += acl.o
ext2-$(CONFIG_EXT2_FS_SECURITY)	 += xattr_security.o
ext2-$(CONFIG_EXT2_FS_XIP)	 += xip.o
