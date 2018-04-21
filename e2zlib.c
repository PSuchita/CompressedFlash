
#include <linux/string.h>
#include <linux/types.h>
#include <linux/fs.h>
#include <linux/ext2_fs_c.h>
#include <linux/module.h>
#include <linux/crypto.h>
#include <linux/zlib.h>
#include <linux/vmalloc.h>

static DEFINE_PER_CPU(struct crypto_comp *, tfm) = NULL;

size_t ext2_iZLIB(int action)
{
	/*mw: we init tfm when we need it...*/
	return 0;
}


size_t ext2_wZLIB(__u8 * ibuf, __u8 * obuf, void *heap,
		  size_t ilen, size_t olen, int level)
{
    int ret, dlen;
  
    if (!try_module_get(THIS_MODULE))
	return 0;
    
    /*check if we already have a tfm*/    
    get_cpu_var(tfm);
    if (__get_cpu_var(tfm) == NULL){
	 __get_cpu_var(tfm) = crypto_alloc_comp("deflate", 0, CRYPTO_ALG_ASYNC);
     }
    assert(__get_cpu_var(tfm) != NULL);			
    
    dlen = olen;
    ret = crypto_comp_compress(__get_cpu_var(tfm) , ibuf, ilen, obuf, &dlen);

    put_cpu_var(tfm);

    if (ret) {
	//printk(KERN_DEBUG "ext2_wZLIB: crypto_comp_compress failed: %d, ilen: %d, olen: %d\n", ret, ilen, olen);
	return 0;
    }
    return dlen;
}


size_t ext2_rZLIB(__u8 * ibuf, __u8 * obuf, void *heap,
		  size_t ilen, size_t olen, int ignored)
{
    int ret, dlen;
  
    if (!try_module_get(THIS_MODULE))
	return 0;
    
    /*check if we already have a tfm*/    
    get_cpu_var(tfm);
    if (__get_cpu_var(tfm) == NULL){
	 __get_cpu_var(tfm) = crypto_alloc_comp("deflate", 0, CRYPTO_ALG_ASYNC);
     }
    assert(__get_cpu_var(tfm) != NULL);				

    dlen = olen;
    ret = crypto_comp_decompress(__get_cpu_var(tfm), ibuf, ilen, obuf, &dlen);

    put_cpu_var(tfm);

    if (ret) {
	//printk(KERN_DEBUG "ext2_wZLIB: crypto_comp_decompress failed: %d, ilen: %d, olen: %d\n", ret, ilen, olen);
	return 0;
    }
	
    return dlen;
}
