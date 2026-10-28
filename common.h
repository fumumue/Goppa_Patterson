#ifndef GOPPA_COMMON_H
#define GOPPA_COMMON_H

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#include "global.h"
#include "struct.h"
#include "simd.h"

#define SIZE_OF_ARRAY(array) (sizeof(array) / sizeof((array)[0]))

unsigned long xor128(RNG *state);

/* common.h */
static inline int mlt(int x, int y)
{
  int s = x + y - 2;
  s -= (s >= M - 1) ? M - 1 : 0;
  return (x && y) ? s + 1 : 0;
}

static inline int mltn(int n, int x)
{
  if (n <= 0) return 1;   /* 元のコードは n<=0 のとき 1 を返す */
  if (x == 0) return 0;   /* 0 の正のべき乗は 0 */
  return (int)(((unsigned long long)n * (unsigned)(x - 1)) % (M - 1)) + 1;
}


int mltn(int n, int x);
void gen_gf(int deg, int order, int opt_sage);
void random_shuffle(unsigned short *array, size_t size);
void print_trace(void);

vec o2v(OP f);
OP v2o(vec a);
OP conv(OP f);
void Pgen(unsigned short permutation[N]);
vec pubkeygen(unsigned short permutation[N],
              MTX *inverse_scrambler,
              RNG *rng,
              unsigned short syndrome_matrix[N][K]);
vec zin(const MTX *inverse_scrambler,
        const unsigned short syndrome_matrix[N][K],
        unsigned short zz[]);
vec paloma(OP w,OP f);
vec patterson(vec w, vec f);
int ero2(vec v, const unsigned short permutation[N]);
void printpol(vec a);
vec vadd(vec a, vec b);
vec vmul(vec a, vec b);
vec fmul(vec a, vec b);
OP oadd(OP f, OP g);
OP omul(OP f, OP g);

int mkg(unsigned short matrix[N][K], vec *result);
int is_reg(const MTX *cc, MTX *R, RNG *rng);
int random_full_rank_binary_matrix(MTX *matrix, RNG *rng);
int random_full_rank_binary_matrix_inv(MTX *S, MTX *Sinv, RNG *rng);
static int check_inverse(const MTX *S, const MTX *Sinv, RNG *rng);

int mkerr(unsigned short *z1, int num, RNG *rng);
int mkS(const MTX *cc, MTX *R);
int mkRS(const MTX *cc, MTX *R);
int binv(const MTX *cc, MTX *L, int Y);
int gauss(const MTX *input, MTX *inverse);
int matinv(const MTX *input, int n, MTX *inverse);
/* Zero-initialize MTX objects before first use; init replaces existing storage. */
int matrix_init(MTX *matrix, size_t rows, size_t cols);
void matrix_free(MTX *matrix);
int matrix_copy(MTX *destination, const MTX *source);
int mulmat(const MTX *A, const MTX *B, int flg, MTX *result);
void mmul(const MTX *A, const MTX *B);
int pubkey_pack(const MTX *scrambler,
                unsigned short field_matrix[N][K],
                const unsigned short permutation[N],
                unsigned short syndrome_matrix[N][K]);
#endif
