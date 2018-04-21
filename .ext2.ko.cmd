cmd_fs/ext2/ext2.ko := ld -r -m elf_x86_64 -T /users/Suchita/linux-3.4.1/scripts/module-common.lds --build-id  -o fs/ext2/ext2.ko fs/ext2/ext2.o fs/ext2/ext2.mod.o
