#ifndef GOPPA_SIMD_H
#define GOPPA_SIMD_H

#include <stddef.h>
#include <stdint.h>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

static inline void xor_u16(unsigned short *dst, const unsigned short *src, size_t len)
{
  size_t i = 0;

#if defined(__AVX2__)
  for (; i + 16 <= len; i += 16)
  {
    __m256i a = _mm256_loadu_si256((const __m256i *)(const void *)(dst + i));
    __m256i b = _mm256_loadu_si256((const __m256i *)(const void *)(src + i));
    _mm256_storeu_si256((__m256i *)(void *)(dst + i), _mm256_xor_si256(a, b));
  }
#endif

  for (; i < len; i++)
    dst[i] ^= src[i];
}

static inline void xor_u16_mod2(unsigned short *dst, const unsigned short *src, size_t len)
{
  size_t i = 0;

#if defined(__AVX2__)
  const __m256i lsb_mask = _mm256_set1_epi16(1);
  for (; i + 16 <= len; i += 16)
  {
    __m256i a = _mm256_loadu_si256((const __m256i *)(const void *)(dst + i));
    __m256i b = _mm256_loadu_si256((const __m256i *)(const void *)(src + i));
    b = _mm256_and_si256(b, lsb_mask);
    _mm256_storeu_si256((__m256i *)(void *)(dst + i), _mm256_xor_si256(a, b));
  }
#endif

  for (; i < len; i++)
    dst[i] ^= src[i] % 2;
}

static inline void xor_u8(uint8_t *dst, const uint8_t *src, size_t len)
{
  size_t i = 0;

#if defined(__AVX2__)
  for (; i + 32 <= len; i += 32)
  {
    __m256i a = _mm256_loadu_si256((const __m256i *)(const void *)(dst + i));
    __m256i b = _mm256_loadu_si256((const __m256i *)(const void *)(src + i));
    _mm256_storeu_si256((__m256i *)(void *)(dst + i), _mm256_xor_si256(a, b));
  }
#endif

  for (; i < len; i++)
    dst[i] ^= src[i];
}

#endif
