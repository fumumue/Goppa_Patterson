
/* -*- mode: C; coding:utf-8 -*- */
#include <xmmintrin.h>
#include <immintrin.h>

//monomial
typedef struct
{
  unsigned short n; //単項式の次数
  unsigned short a; //単項式の係数
} oterm;

//polynomial
typedef struct
{
  oterm t[DEG]; //単項式の配列として多項式を表現する
} OP;

typedef union 
{
  unsigned short x[DEG]; //配列の添字を次数に、配列の値を係数に持つ多項式の表現
  unsigned long long d[DEG/4];
  __m256 a[DEG/16];
} vec;

typedef struct {
  OP q;
  OP r;
} rem;

typedef struct
{
  unsigned short v[N];
  int f;
} MT;

//extra gcd
typedef struct
{
  OP u; //inverse of polynomial?
  OP v; //error locater
 OP d; //gcd
} EX;

typedef union
{ //test(SIMN)
  unsigned long long int u[K/4];
  unsigned short s[K];
} SU;

typedef union
{
  __uint128_t z[K/16];
  unsigned long long int u[K/8];
  unsigned int t[K/4];
  unsigned short x[K/2];
  unsigned char d[K];
} arrayul;

typedef union
{
  unsigned long long int u[N];
  unsigned char d[N * 8];
} arrayull;

typedef union a4
{
  unsigned char ar[4];
  unsigned int n;
} array;

typedef struct a8
{
  unsigned char ar[8];
} array8;

typedef union
{
  unsigned int h[16];
  unsigned char c[64];
} array16;

typedef union aN
{
  unsigned int d[64];
  unsigned long long int u[32];
  unsigned char ar[256];
  //
} arrayn;

typedef struct pub
{
  unsigned char a[K];
  unsigned char b[K];
} set;


#define I8T char
#define U8C(v) (v##U)

#define U8V(v) ((unsigned char)(v)&U8C(0xFF))
#define ROTL8(v, n) \
  (U8V((v) << (n)) | ((v) >> (8 - (n))))

#define R(x, n) (((x) << (n)) | ((x) >> (32 - (n))))

/*
unsigned int rotate_left(unsigned int x, int n)
{
  assert(0 < n && n < 32);
  return (x << n) | (x >> (32 - n));
}


typedef union
{
  unsigned long long int u[8];
  unsigned char d[64];
} arrayul;

typedef struct a4
{
  unsigned char ar[4];
} array;

typedef struct a8
{
  unsigned char ar[8];
} array8;

typedef struct
{
  unsigned int h[16];
} array16;

typedef struct aN
{
  unsigned char ar[N];
} arrayn;

typedef struct pub
{
  unsigned char a[N];
  unsigned char b[N];
} set;
*/

typedef struct {
vec v[N];
} ext;

typedef union 
{
  unsigned short x[N][N];
  vec a[N];
} maya;

typedef struct {
  int row; //行
  int col; //列
} ja;
 
typedef union {
  unsigned short x[16];
  unsigned long long d[4];
  __m256 a;
} Uh;

typedef union {
  unsigned short x[N][N];
  unsigned long long d[N][N/4];
  __m256 c[N][N/16];
  vec a[N];
  ja b;
} MTX;

/*
typedef struct {
  unsigned short x[N][K];
  unsigned char z[N][K*E];
  unsigned char w[K*E][K*E];
  int i;
  int rank;
} MAT;
*/