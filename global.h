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
#define K (128) //符号の次元
#define E (13)    //拡大体のビット数

/*
 * This implementation is intentionally fixed-parameter.
 * The original code hard-codes special-case reductions for a few K values
 * and assumes K <= 32 with a K*E-sized check matrix.
 */
#if K > 256
#error "Unsupported K: this build is only safe for K <= 32."
#endif

#if (K % 2) != 0
#error "Unsupported K: use an even value for T = K / 2."
#endif

#if ((K) * (E)) > N
#error "Unsupported configuration: K * E must not exceed N."
#endif

#define DEG K*E //(K * E)
#define T (K/2) //エラーの数

#define GF_TABLE_SIZE 32768

extern unsigned short gf[GF_TABLE_SIZE];
extern unsigned short fg[GF_TABLE_SIZE];

#endif
