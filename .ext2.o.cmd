cmd_fs/ext2/ext2.o := ld -m elf_x86_64   -r -o fs/ext2/ext2.o fs/ext2/balloc.o fs/ext2/dir.o fs/ext2/file.o fs/ext2/ialloc.o fs/ext2/inode.o fs/ext2/ioctl.o fs/ext2/namei.o fs/ext2/super.o fs/ext2/symlink.o fs/ext2/adler32.o fs/ext2/compress.o fs/ext2/e2zlib.o 