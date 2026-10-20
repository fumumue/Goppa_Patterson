#ifndef GOPPA_GLOBAL_H
#define GOPPA_GLOBAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//符号のパラーメータの指定。通常[N,K,T]として、
//Nは符号の長さ、Kが符号の次元、Tは訂正エラー数
//を表す。ここではDは符号長にしている。
#define N 256 // 符号長
#define M N //有限体の元の数
#define K (32) //符号の次元
#define E (8)    //拡大体のビット数
#define DEG K*E //(K * E)
#define T (K/2) //エラーの数

#define GF_TABLE_SIZE 32768

extern unsigned short gf[GF_TABLE_SIZE];
extern unsigned short fg[GF_TABLE_SIZE];
extern unsigned short P[N];
extern unsigned short inv_P[N];
extern unsigned short uu;

#endif
