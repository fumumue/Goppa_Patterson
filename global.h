#ifndef GOPPA_GLOBAL_H
#define GOPPA_GLOBAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//符号のパラメータの指定。通常[N,K,T]として、
//Nは符号の長さ、Kが符号の次元、Tは訂正エラー数
//を表す。ここではDは符号長にしている。
#define N 8192 // 符号長
#define M N //有限体の元の数
#define K (256) //符号の次元
#define E (13)    //拡大体のビット数
#define GF_TABLE_SIZE 32768

#define NB ((size_t)K * E)            // 3328 行(シンドローム長)
#define RW ((size_t)N / 64)           // 128 語: 長さ N の行
#define TW ((size_t)(N - NB) / 64)    // 76 語: T の行 (4864 ビット)
_Static_assert(((size_t)K * E) % 64 == 0, "K*E must be a multiple of 64");
_Static_assert(N % 64 == 0, "N must be a multiple of 64");


#if (K != 4) && (K != 8) && (K != 16) && (K != 32) && \
    (K != 64) && (K != 96) && (K != 128) && (K != 256) && (K != 512)
#error "Unsupported K: no field reduction polynomial is defined for this value."
#endif

#if (K % 2) != 0
#error "Unsupported K: use an even value for T = K / 2."
#endif

#if N > GF_TABLE_SIZE
#error "Unsupported N: the configured finite-field tables are too small."
#endif

#if ((K) * (E)) > N
#error "Unsupported configuration: K * E must not exceed N."
#endif

#define DEG K*E //(K * E)
#define T (K/2) //エラーの数

extern unsigned short gf[GF_TABLE_SIZE];
extern unsigned short fg[GF_TABLE_SIZE];

#endif
