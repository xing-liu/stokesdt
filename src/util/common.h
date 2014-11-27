/**
 * @file   common.h
 * @brief  Platform-specific constants and functions.
 */

#ifndef COMMON_H_
#define COMMON_H_


#include <malloc.h>


namespace stokesdt {

namespace detail {
    
/// the SIMD width
#if defined (__MIC__)
const int kSimdWidth = 64;
#elif defined (__AVX__)
const int kSimdWidth = 32;
#elif defined (__SSE__)
const int kSimdWidth = 16;
#else
const int kSimdWidth = 64;
#endif


/// the alignment length
const int kAlignLen = 64;


/// Allocates aligned memory
inline void *AlignMalloc(size_t size)
{
#ifdef __INTEL_COMPILER
    void * addr = _mm_malloc(size, kAlignLen);
#else
    void * addr = memalign(kAlignLen, size);
#endif
    return addr;
}


/// Destroys aligned memory
inline void AlignFree(void *addr)
{
    if (addr != NULL) {
    #ifdef __INTEL_COMPILER
        _mm_free(addr);
    #else
        free(addr);
    #endif
    }
}


/** \brief  Returns the smallest integer that is a multiple
 *  of <code>size</code> and larger than <code>len</code>.
 */
inline size_t PadLen(size_t len, size_t size)
{
    size_t len0 = detail::kAlignLen/size;
    return ((len + len0 - 1)/len0 * len0);
}


/// the maximum line length
const int kMaxLine = 1024;

} // namespace detail

} // namespace stokesdt


/// Disable copy constructor and assignement operator.
#define DISALLOW_COPY_AND_ASSIGN(TypeName)   \
    TypeName(const TypeName&);               \
    void operator=(const TypeName&)


#endif // COMMON_H_
