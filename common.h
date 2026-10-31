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
int mlt(int x, int y);
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

int is_reg(const MTX *cc, MTX *R, RNG *rng);
int random_full_rank_binary_matrix(MTX *matrix, RNG *rng);
int mkerr(unsigned short *z1, int num, RNG *rng);
int mkS(const MTX *cc, MTX *R);
int mkRS(const MTX *cc, MTX *R);
int binv(const MTX *cc, MTX *L, int Y);
MTX gauss(const MTX *input);
MTX matinv(const MTX *input, int n);
MTX mulmat(const MTX *A, const MTX *B, int flg);
void mmul(const MTX *A, const MTX *B);

#endif
