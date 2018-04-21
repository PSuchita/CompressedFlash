cmd_fs/ext2/adler32.o := gcc -Wp,-MD,fs/ext2/.adler32.o.d  -nostdinc -isystem /usr/lib/gcc/x86_64-linux-gnu/4.8/include -I/users/Suchita/linux-3.4.1/arch/x86/include -Iarch/x86/include/generated -Iinclude  -include /users/Suchita/linux-3.4.1/include/linux/kconfig.h -D__KERNEL__ -Wall -Wundef -Wstrict-prototypes -Wno-trigraphs -fno-strict-aliasing -fno-common -Werror-implicit-function-declaration -Wno-format-security -fno-delete-null-pointer-checks -O2 -m64 -mtune=generic -mno-red-zone -mcmodel=kernel -funit-at-a-time -maccumulate-outgoing-args -fstack-protector -DCONFIG_X86_X32_ABI -DCONFIG_AS_CFI=1 -DCONFIG_AS_CFI_SIGNAL_FRAME=1 -DCONFIG_AS_CFI_SECTIONS=1 -DCONFIG_AS_FXSAVEQ=1 -pipe -Wno-sign-compare -fno-asynchronous-unwind-tables -mno-sse -mno-mmx -mno-sse2 -mno-3dnow -mno-avx -Wframe-larger-than=1024 -Wno-unused-but-set-variable -fno-omit-frame-pointer -fno-optimize-sibling-calls -g -pg -Wdeclaration-after-statement -Wno-pointer-sign -fno-strict-overflow -fconserve-stack -DCC_HAVE_ASM_GOTO  -DMODULE  -D"KBUILD_STR(s)=\#s" -D"KBUILD_BASENAME=KBUILD_STR(adler32)"  -D"KBUILD_MODNAME=KBUILD_STR(ext2)" -c -o fs/ext2/.tmp_adler32.o fs/ext2/adler32.c

source_fs/ext2/adler32.o := fs/ext2/adler32.c

deps_fs/ext2/adler32.o := \

fs/ext2/adler32.o: $(deps_fs/ext2/adler32.o)

$(deps_fs/ext2/adler32.o):
