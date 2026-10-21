// Niederreiter cryptosystem with Patterson decoding.

#include "common.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#define REG 64
extern int mlt(int x, int y);
extern int mltn(int n, int a);
extern void makeS();

// 有限体の元の逆数
unsigned short
oinv(unsigned short a)
{
  if (a == 0)
    return -1;

  return (N - fg[a]) % (N-1) + 1;

  printf("no return \n");

  exit(1);
}

// 有限体の元の逆数
unsigned short
oinv2(unsigned short a)
{
  int i;

  if (a == 0)
    return -1;

  for (i = 0; i < N; i++)
  {
    if (gf[mlt(fg[a], i)] == 1)
      return (unsigned short)i;
  }

  printf("no return \n");
  exit (1);
}

// aに何をかけたらbになるか
unsigned short
equ(unsigned short a, unsigned short b)
{
  int i;

  for (i = 0; i < N; i++)
  {
    if (gf[mlt(fg[a], fg[i])] == b)
      break;
  }
  return i;
}

// 停止コマンド
void wait(void)
{
  int a;                                     // 読み込む変数はローカルに取るべし
  printf(" (enter number and hit return) "); // 何か表示させたほうが良いだろう
  fflush(stdout);                            // just in case
  if (scanf("%d", &a) != 1)
    fprintf(stderr, "Failed to read wait input.\n");
}

// 多項式の次数(default)
int deg(vec a)
{
  int i, n = 0, flg = 0;

  // #pragma omp parallel for
  for (i = 0; i < DEG; i++)
  {
    if (a.x[i] > 0)
    {
      n = i;
      flg = 1;
    }
  }
  if (flg == 0)
    return 0;

  return n;
}


vec vadd(vec a, vec b)
{
  vec c = {0};
  xor_u16(c.x, a.x, DEG);
  xor_u16(c.x, b.x, DEG);

  return c;
}

vec vmul(vec a, vec b)
{
  int i, j, k, l;
  vec c = {0};

  k = deg(a);
  l = deg(b);

  for (i = 0; i < k; i++)
  {
    for (j = 0; j < l; j++)
    {
      if (a.x[i] > 0)
      {
        c.x[i + j] ^= gf[mlt(fg[a.x[i]], fg[b.x[j]])];
      }
    }
  }

  return c;
}

vec fmul(vec a, vec b)
{
  int i, j;
  vec c = {0};

  for (i = 0; i < DEG; i++)
  {
    if (a.x[i] == 0)
      continue;

    for (j = 0; j < DEG - i; j++)
    {
      if (b.x[j] > 0)
        c.x[i + j] ^= gf[mlt(fg[a.x[i]], fg[b.x[j]])];
    }
  }

  return c;
}

static int poly_degree(vec f)
{
  for (int i = DEG - 1; i >= 0; i--)
    if (f.x[i] != 0)
      return i;
  return -1;
}

static unsigned short leading_coefficient(vec f)
{
  int degree = poly_degree(f);
  return degree >= 0 ? f.x[degree] : 0;
}

static unsigned short lowest_coefficient(vec f)
{
  for (int i = 0; i < DEG; i++)
    if (f.x[i] != 0)
      return f.x[i];
  return 0;
}

vec poly_add(vec f, vec g)
{
  return vadd(f, g);
}

vec poly_mul(vec f, vec g)
{
  return fmul(f, g);
}

static vec coeff(vec f, unsigned short divisor)
{
  unsigned short inverse;

  if (divisor == 0)
  {
    fprintf(stderr, "cannot normalize a polynomial by zero\n");
    exit(EXIT_FAILURE);
  }
  inverse = oinv(divisor);
  for (int i = 0; i < DEG; i++)
    if (f.x[i] != 0)
      f.x[i] = gf[mlt(fg[f.x[i]], fg[inverse])];
  return f;
}

// 多項式を表示する(default)
void printpol(vec a)
{
  int i, n;

  n = deg(a);

  // printf ("baka\n");
  assert(n >= 0);

  for (i = n; i > -1; i--)
  {
    if (a.x[i] > 0)
    {
      printf("%u", a.x[i]);
      // if (i > 0)
      printf("x^%d", i);
      // if (i > 0)
      printf("+");
    }
  }
  //  printf("\n");

  return;
}

vec poly_div(vec f, vec g)
{
  vec quotient = {0};
  int divisor_degree = poly_degree(g);

  if (divisor_degree < 0)
  {
    fprintf(stderr, "polynomial division by zero\n");
    exit(EXIT_FAILURE);
  }
  while (poly_degree(f) >= divisor_degree)
  {
    int degree = poly_degree(f);
    int shift = degree - divisor_degree;
    unsigned short factor = (unsigned short)equ(g.x[divisor_degree], f.x[degree]);
    quotient.x[shift] ^= factor;
    for (int i = 0; i <= divisor_degree; i++)
      if (g.x[i] != 0)
        f.x[i + shift] ^= gf[mlt(fg[factor], fg[g.x[i]])];
  }
  return quotient;
}

vec poly_mod(vec f, vec g)
{
  int divisor_degree = poly_degree(g);

  if (divisor_degree < 0)
  {
    fprintf(stderr, "polynomial modulo zero\n");
    exit(EXIT_FAILURE);
  }
  while (poly_degree(f) >= divisor_degree)
  {
    int degree = poly_degree(f);
    int shift = degree - divisor_degree;
    unsigned short factor = (unsigned short)equ(g.x[divisor_degree], f.x[degree]);
    for (int i = 0; i <= divisor_degree; i++)
      if (g.x[i] != 0)
        f.x[i + shift] ^= gf[mlt(fg[factor], fg[g.x[i]])];
  }
  return f;
}

static unsigned short trace(vec f, unsigned short x)
{
  unsigned short result = 0;

  for (int i = 0; i < DEG; i++)
    if (f.x[i] != 0)
      result ^= gf[mlt(fg[f.x[i]], mltn(i, fg[x]))];
  return result;
}

vec poly_inverse(vec a, vec modulus)
{
  vec r0 = modulus, r1 = a;
  vec t0 = {0}, t1 = {0};
  t1.x[0] = 1;

  while (poly_degree(r1) >= 0)
  {
    vec quotient = poly_div(r0, r1);
    vec remainder = poly_mod(r0, r1);
    vec next = poly_add(t0, poly_mul(quotient, t1));
    r0 = r1;
    r1 = remainder;
    t0 = t1;
    t1 = next;
  }

  if (poly_degree(r0) != 0)
  {
    fprintf(stderr, "polynomial has no inverse modulo the supplied modulus\n");
    exit(EXIT_FAILURE);
  }
  return poly_mod(coeff(t0, r0.x[0]), modulus);
}

unsigned short
v2a(unsigned short a)
{
  int j;

  if (a == 0)
    return 0;

  // printf("aa=%d\n",a.a);
  for (j = 0; j < N; j++)
  {
    if (gf[j] == a && a > 0)
    {
      // printf("j==%d\n",j);
      return j - 1;
    }
  }
  fprintf(stderr, "Field element %u is not present in the field table.\n", a);
  exit(EXIT_FAILURE);
}

void printsage(vec a)
{
  int i, j;

  printf("poly=");
  for (i = 0; i < DEG; i++)
  {
    if (a.x[i] > 0)
    {
      j = v2a(a.x[i]);
      // printf("%d,==ba\n",b.a);
      // printf ("X**%d+", i); //for GF2
      printf("B('a^%d')*X**%d+", j, i); // for GF(2^m)
    }
  }
}

static vec poly_gcd_bounded(vec a, vec n, int limit)
{
  vec r0 = n, r1 = a;
  vec t0 = {0}, t1 = {0};
  t1.x[0] = 1;

  while (poly_degree(r1) > limit)
  {
    vec quotient = poly_div(r0, r1);
    vec remainder = poly_mod(r0, r1);
    vec next = poly_add(t0, poly_mul(quotient, t1));
    r0 = r1;
    r1 = remainder;
    t0 = t1;
    t1 = next;
  }
  return t1;
}

static vec poly_gcd(vec a, vec b)
{
  while (poly_degree(b) >= 0)
  {
    vec remainder = poly_mod(a, b);
    a = b;
    b = remainder;
  }
  if (poly_degree(a) < 0)
    return a;
  return coeff(a, a.x[poly_degree(a)]);
}

// 整数からベクトル型への変換
vec i2v(unsigned short n)
{
  vec v = {0};
  int i = 0;

  while (n > 0)
  {
    v.x[i++] = n % 2;
    n = (n >> 1);
  }

  return v;
}

// ベクトル型から整数への変換
unsigned short
v2i(vec v)
{
  unsigned int d = 0, i, e = 0;

  for (i = 0; i < deg(v) + 1; i++)
  {
    e = v.x[i];
    d ^= (e << i);
  }

  return d;
}

// 配列からベクトル表現の多項式へ変換する
vec Setvec(const unsigned short coefficients[], int n)
{
  int i;
  vec v = {0};

  for (i = 0; i < n; i++)
  {
    v.x[n - 1 - i] = coefficients[i];
  }

  return v;
}

static vec chen(vec f)
{
  vec errors = {0};
  int count = 0;

  for (int x = 0; x < N; x++)
  {
    unsigned short value = 0;
    for (int i = 0; i < DEG; i++)
      if (f.x[i] != 0)
        value ^= gf[mlt(mltn(i, fg[x]), fg[f.x[i]])];
    if (value == 0 && count < DEG)
      errors.x[count++] = (unsigned short)x;
  }
  return errors;
}

static int ben_or(vec f)
{
  return poly_degree(f) < 0 ? -1 : 0;
}

static vec setpol(unsigned short coefficients[], int n)
{
  return Setvec(coefficients, n);
}

// バイナリ型パリティチェック行列を生成する
static void bdet(const unsigned short matrix[N][K * E], MTX *check_matrix)
{
  int i, j, k, l;

  for (i = 0; i < N; i++)
  {
    for (j = 0; j < K; j++)
    {
      l = matrix[i][j];
      for (k = 0; k < E; k++)
      {
        check_matrix->x[i][j * E + k] = l % 2;
        l = (l >> 1);
      }
    }
  }

  for (i = 0; i < N; i++)
  {
    for (j = 0; j < E * K; j++)
      printf("%d,", check_matrix->x[i][j]);
    printf("\n");
  }
}

void toByte(const MTX *SH, unsigned short syndrome_matrix[N][K])
{
  int i, j, k;

  for (i = 0; i < N; i++)
  {
    printf("数え上げ %d\n", i);
    // #pragma omp parallel for
    for (j = 0; j < K; j++)
    {
    unsigned short x=0;
        k = j * E;
        x=0;
        for (int l = 0; l < E; l++){
          x ^= (SH->x[i][k + l]<<l);
          }
        syndrome_matrix[i][j] = x;
    }
  }
  printf("end of byte\n");
}

// 秘密置換を生成する
void Pgen(unsigned short permutation[N])
{
  unsigned int i;

  for (i = 0; i < N; i++)
    permutation[i] = i;
  random_shuffle(permutation, N);
}

// 有限体の元の平方を計算する
int isqrt(unsigned short u)
{
  int i;

  for (i = 0; i < N; i++)
  {
    if (gf[mlt(i, i)] == u)
      return i;
  }

  printf("来ちゃいけないところに来ました\n");
  exit(1);
}

// 多項式の平方を計算する
static vec osqrt(vec f, vec w)
{
  vec even = {0}, odd = {0}, h = {0}, r = {0}, ww = {0}, s = {0}, tmp = {0};

  for (int i = 0; i < DEG; i++)
  {
    if (f.x[i] == 0)
      continue;
    if ((i & 1) == 0)
    {
      even.x[i / 2] = gf[isqrt(f.x[i])];
    }
    else
    {
      odd.x[(i - 1) / 2] = gf[isqrt(f.x[i])];
    }
  }

  for (int i = 0; i < DEG; i++)
  {
    if (w.x[i] == 0)
      continue;
    if ((i & 1) == 0)
      h.x[i / 2] = gf[isqrt(w.x[i])];
    else
      r.x[(i - 1) / 2] = gf[isqrt(w.x[i])];
  }

  printpol(r);
  printf(" sqrt(g1)=======\n");
  if (poly_degree(r) > 0)
  {
    s = poly_inverse(r, w);
  }
  else if (poly_degree(r) <= 0)
  {
    printpol(r);
    printf(" deg(r)======0!\n");
    printpol(w);
    printf(" goppa======0\n");
    printpol(f);
    printf(" syn======0\n");
    if (leading_coefficient(r) > 0)
    {
      s.x[0] = gf[isqrt(leading_coefficient(r))];
      return poly_mod(poly_add(even, poly_mul(coeff(h, s.x[0]), odd)), w);
    }
    return even;
  }

  if (poly_degree(s) > 0)
    tmp = poly_mod(poly_mul(s, r), w);
  if (poly_degree(tmp) > 0)
  {
    fprintf(stderr, "square-root polynomial component is not invertible\n");
    exit(EXIT_FAILURE);
  }
  if (poly_degree(h) > 0 && poly_degree(s) > 0)
    ww = poly_mod(poly_mul(h, s), w);
  if (poly_degree(h) <= 0 || poly_degree(s) <= 0)
  {
    printpol(h);
    fprintf(stderr, "invalid polynomial square-root components\n");
    exit(EXIT_FAILURE);
  }

  if (poly_degree(ww) < 0)
  {
    fprintf(stderr, "polynomial square-root component is zero\n");
    exit(EXIT_FAILURE);
  }

  tmp = poly_mod(poly_mul(ww, ww), w);
  if (poly_degree(tmp) == 1)
  {
    printpol(ww);
    printf(" ww succsess!===========\n");
  }
  else
  {
    fprintf(stderr, "polynomial square-root verification failed\n");
    exit(EXIT_FAILURE);
  }

  h = ww;
  if (poly_degree(poly_mod(poly_mul(h, ww), w)) == 1)
    return poly_mod(poly_add(even, poly_mul(ww, odd)), w);
  fprintf(stderr, "unexpected polynomial square-root result\n");
  exit(EXIT_FAILURE);
}

// パターソンアルゴリズムでバイナリGoppa符号を復号する
vec patterson(vec w_vec, vec f_vec)
{
  vec w = w_vec, f = f_vec;
  vec g1 = {0}, ll = {0};
  int flg2 = 0;
  vec tt = {0}, ff = {0}, h = {0};
  vec v;
  vec r2 = {0}, b2 = {
                   0};
  // unsigned short g[K+1]={2,2,12,1,2,8,4,13,5,10,8,2,15,10,7,3,5};

  tt.x[1] = 1;

  ff = poly_inverse(f, w);
  b2 = poly_mod(poly_mul(ff, f), w);
  if (poly_degree((b2)) > 0)
  {
    printf("逆元が計算できません。error\n");
    // wait ();
    exit(1);
  }
  printf("locater==========\n");
  // exit(1);
  r2 = poly_add(ff, tt);
  printpol(r2);
  printf(" h+x==============\n");
  // wait();
  //   exit(1);
  g1 = osqrt(r2, w);
  printpol(g1);
  printf(" sqrt(h+x)=======\n");
  b2 = poly_mod(poly_mul(g1, g1), w);
  printpol(b2);
  printf(" g1*g1%%w===========\n");
  if (deg(b2) == 0)
  {
    printf("b2 is 0\n");
    exit(1);
  }
  if (lowest_coefficient(b2) != lowest_coefficient(r2))
  {
    printpol(w);
    printf(" w============\n");
    printpol(r2);
    printf(" r2============\n");
    printpol(g1);
    printf(" g1============\n");
    printf(" g1は平方ではありません。error");
    // wait ();
    exit(1);
  }
  printpol(w);
  printf(" ========goppa\n");
  printpol(g1);
  printf(" g1!=========\n");
  if (poly_degree(g1) < 0)
  {
    printf(" badkey=========\n");
    printf(" g1============\n");
    printf("平方が０になりました。 error\n");
    // wait ();
    exit(1);
  }
  h = poly_gcd(f, w);
  if (poly_degree((h)) > 0)
  {
    printsage(w);
    printf(" =goppa\n");
    printsage(f);
    printf(" =syn\n");
    printpol(h);
    printf(" =poly_gcd\n");

    // printf("f-ben_or=%d\n",ben_or(f));
    // printf("w-ben_or=%d\n",ben_or(w));
    printf(" s,wは互いに素じゃありません。\n");
    flg2 = 1;
    // wait ();
    exit(1);
    // goto label;
  }

  // exit(1);
  vec ppo = {0};

  // hh = xgcd(w, g1);
  ppo = poly_gcd_bounded(w, g1, T);
  printpol(ppo);
  printf(" =========ppo\n");
  printpol(ppo);
  printf(" =========ppo\n");
  // if(LT(ppo).a!=LT(ppo).a)
  // exit(1);
  printpol(ppo);
  printf(" =========ppo\n");
  ff = poly_mod(poly_mul(ppo, g1), w);
  printpol(ppo);
  printf(" alpha!=========\n");
  printpol(ff);
  printf(" beta!=========\n");
  printpol(w);
  printf(" goppa=========\n");
  printpol(f);
  printf(" syn=========\n");

  printpol(ppo);
  printf(" alpha!=========\n");
  // hh = xgcd(w, g1);
  ppo = poly_gcd_bounded(w, g1, T);
  //  h=poly_gcd_bounded(w,g1);
  ff = poly_mod(poly_mul(ppo, g1), w);
  ll = poly_add(poly_mul(ff, ff), poly_mul(tt, poly_mul(ppo, ppo)));
  v = chen(ll);
  if (v.x[K - 1] > 0)
  {
    // wait ();
    return v;
  }

  ppo = poly_gcd_bounded(w, g1, T);
  ff = poly_mod(poly_mul(ppo, g1), w);
  if (poly_degree((ff)) != K / 2)
  {

    printf("\nbefore h.d\n");
    ff = poly_mod(poly_mul(ppo, g1), w);
    printpol(ff);
    printf(" ==========beta!\n");
    printpol(ppo);
    printf(" alpha!=========\n");
    ll = poly_add(poly_mul(ff, ff), poly_mul(tt, poly_mul(ppo, ppo)));
    v = chen(ll);
    if (v.x[K - 1] > 0)
    {
      wait();
      return v;
    }
  }

  if (poly_degree((ff)) == 1)
  {
    ll = poly_add(poly_mul(ff, ff), poly_mul(tt, poly_mul(ppo, ppo))); // ff;
    printf("deg==1\n");
    exit(1);
  }

  printf("あっ、でる・・・！\n");
  printpol(ll);
  printf(" ll=================\n");
  printpol(f);
  printf(" syn=================\n");
  v = chen(ll);

  if (v.x[K - 1] > 0)
  {
    return v;
  }

  printf("flg2=%d\n", flg2);
  // exit(1);
  // ero2(v);
  return v;
}

static void van(unsigned short powers[K * 2][N])
{
  int i, j;

  printf("van der\n");

  for (i = 0; i < N; i++)
    powers[0][i] = 1;
  for (i = 1; i < K; i++)
  {
    for (j = 0; j < N; j++)
      powers[i][j] = gf[mltn(i, fg[j])];
  }
}

static void ogt(const unsigned short polynomial[K + 1], unsigned short table[K * 2][K * 2])
{
  int i, j;

  for (i = 0; i < K; i++)
  {
    for (j = 0; j < K - i; j++)
      table[i][j + i] = polynomial[j];
  }
}

/* input: in0, in1 in GF((2^m)^t)*/
/* output: out = in0*in1 */
void GF_mul(unsigned short *out, unsigned short *in0, unsigned short *in1)
{
  int i, j;

  unsigned short prod[K * 2 - 1] = {0};

  for (i = 0; i < K * 2 - 1; i++)
    prod[i] = 0;

  for (i = 0; i < K; i++)
  {
    for (j = 0; j < K; j++)
      prod[i + j] ^= gf[mlt(fg[in0[i]], fg[in1[j]])];
  }
  //

  for (i = (K - 1) * 2; i >= K; i--)
  {
    if (K == 512)
    {
      // GF(2^512) from sage
      prod[i - K + 8] ^= prod[i];
      prod[i - K + 5] ^= prod[i];
      prod[i - K + 2] ^= prod[i];
      prod[i - K + 0] ^= prod[i];
    }
    if (K == 256)
    {
      // GF(2^256) from sage
      prod[i - K + 10] ^= prod[i];
      prod[i - K + 5] ^= prod[i];
      prod[i - K + 2] ^= prod[i];
      prod[i - K + 0] ^= prod[i];
    }
    if (K == 128)
    {
      // 128
      prod[i - K + 7] ^= prod[i];
      prod[i - K + 2] ^= prod[i];
      prod[i - K + 1] ^= prod[i];
      prod[i - K + 0] ^= prod[i];
    }
             if (K == 96)
      {
        prod[i - K + 63] ^= prod[i];
        prod[i - K + 61] ^= prod[i];
        prod[i - K + 60] ^= prod[i];
        prod[i - K + 59] ^= prod[i];
        prod[i - K + 54] ^= prod[i];
        prod[i - K + 50] ^= prod[i];
        prod[i - K + 46] ^= prod[i];
        prod[i - K + 44] ^= prod[i];
        prod[i - K + 40] ^= prod[i];
        prod[i - K + 37] ^= prod[i];
        prod[i - K + 33] ^= prod[i];
        prod[i - K + 32] ^= prod[i];
        prod[i - K + 28] ^= prod[i];
        prod[i - K + 24] ^= prod[i];
        prod[i - K + 23] ^= prod[i];
        prod[i - K + 21] ^= prod[i];
        prod[i - K + 20] ^= prod[i];
        prod[i - K + 19] ^= prod[i];
        prod[i - K + 18] ^= prod[i];
        prod[i - K + 16] ^= prod[i];
        prod[i - K + 15] ^= prod[i];
        prod[i - K + 14] ^= prod[i];
        prod[i - K + 13] ^= prod[i];
        prod[i - K + 11] ^= prod[i];
        prod[i - K + 9] ^= prod[i];
        prod[i - K + 8] ^= prod[i];
        prod[i - K + 7] ^= prod[i];
        prod[i - K + 6] ^= prod[i];
        prod[i - K + 4] ^= prod[i];
        prod[i - K + 3] ^= prod[i];
        prod[i - K + 2] ^= prod[i];
        prod[i - K + 0] ^= prod[i];
      }
      if (K == 64)
      {
        prod[i - K + 33] ^= prod[i];
        prod[i - K + 30] ^= prod[i];
        prod[i - K + 26] ^= prod[i];
        prod[i - K + 25] ^= prod[i];
        prod[i - K + 24] ^= prod[i];
        prod[i - K + 23] ^= prod[i];
        prod[i - K + 22] ^= prod[i];
        prod[i - K + 21] ^= prod[i];
        prod[i - K + 20] ^= prod[i];
        prod[i - K + 18] ^= prod[i];
        prod[i - K + 13] ^= prod[i];
        prod[i - K + 12] ^= prod[i];
        prod[i - K + 11] ^= prod[i];
        prod[i - K + 10] ^= prod[i];
        prod[i - K + 7] ^= prod[i];
        prod[i - K + 5] ^= prod[i];
        prod[i - K + 4] ^= prod[i];
        prod[i - K + 2] ^= prod[i];
        prod[i - K + 1] ^= prod[i];
        prod[i - K + 0] ^= prod[i];
      }
    if (K == 32)
    {
      // 32
      prod[i - K + 15] ^= prod[i];
      prod[i - K + 9] ^= prod[i];
      prod[i - K + 7] ^= prod[i];
      prod[i - K + 4] ^= prod[i];
      prod[i - K + 3] ^= prod[i];
      prod[i - K + 0] ^= prod[i];
    }
    if (K == 16)
    {
      // 16
      prod[i - K + 5] ^= prod[i];
      prod[i - K + 3] ^= prod[i];
      prod[i - K + 2] ^= prod[i];
      prod[i - K + 0] ^= prod[i];
    }
    if (K == 8)
    {
      // 8
      prod[i - K + 4] ^= prod[i];
      prod[i - K + 3] ^= prod[i];
      prod[i - K + 2] ^= prod[i];
      prod[i - K + 0] ^= prod[i];
    }
    if (K == 4)
    {
      // 4
      prod[i - K + 1] ^= prod[i];
      prod[i - K + 0] ^= prod[i];
    }
  }

  for (i = 0; i < K; i++)
    out[i] = prod[i];
}

// #define NN 16
vec renritu(const MTX *input)
{
  unsigned short p, d;
  int i, j, k;
  MTX a = *input;
  vec v = {0};

  for (i = 0; i < K; i++)
  {
    p = a.x[i][i];

    for (j = 0; j < (K + 1); j++)
    {
      a.x[i][j] = gf[mlt(fg[a.x[i][j]], oinv(p))];
    }

    for (j = 0; j < K; j++)
    {
      if (i != j)
      {
        d = a.x[j][i];

        for (k = i; k < (K + 1); k++)
        {
          a.x[j][k] = a.x[j][k] ^ gf[mlt(fg[d], fg[a.x[i][k]])];
        }
      }
    }
  }
  for (i = 0; i < K; i++)
  {
    if (a.x[i][i] != 1)
    {
      for (j = 0; j < K + 1; j++)
        printf("a%d,", a.x[i][j]);
      printf("\n");
    }
  }
  printf("\n");
  vec x = {0};
  for (i = 0; i < K; i++)
  {
    v.x[i] = a.x[i][K];
    // v.x[128]=1;
    printf(" x%d = %d\n", i, v.x[i]);
    x.x[i + 1] = v.x[i];
  }
  x.x[0] = 1;

  vec pol = {0};
  pol = setpol(x.x, K + 1);
  printpol(pol);
  printf(" ==key\n");
  for (i = 0; i < N; i++)
  {
    if (trace(pol, i) == 0)
      printf("%d i=%d\n", i, fg[i] - 1);
  }

  return v;
}

/* input: f, element in GF((2^m)^t) */
/* output: out, minimal polynomial of f */
/* return: 0 for success and -1 for failure */
int mykey(unsigned short *out, vec x)
{
  unsigned short mat[K + 1][K] = {0};
  MTX a = {0};
  int i, j;

  // fill matrix

  mat[0][0] = 1;

  for (i = 1; i < K; i++)
    mat[0][i] = 0;

  for (i = 0; i < K; i++)
    mat[1][i] = x.x[i];

  for (j = 2; j <= K; j++)
  {
    GF_mul(mat[j], mat[j - 1], x.x);
  }
  // exit(1);
  //
  for (i = 0; i < K; i++)
  {
    for (j = 0; j < K + 1; j++)
    {
      a.x[i][j] = mat[j][i];
      printf("%d,", mat[j][i]);
    }
    printf("\n");
  }
  printf("\n");
  // exit(1);

  vec v = {0};
  v = renritu(&a);

  for (i = 0; i < K; i++)
  {
    out[i] = v.x[i];
    printf("%d,", out[i]);
  }
  printf("\n");
  return 0;
}

static int mkg(unsigned short matrix[N][K * E], vec *result)
{
  int i, j, k, l;
  vec w = {0};
  unsigned short tr[N] = {0};
  unsigned short ta[N] = {0};
  unsigned short ma[N][K * E] = {0};
  unsigned short powers[K * 2][N] = {0};
  unsigned short polynomial_table[K * 2][K * 2] = {0};

  // w=setpol(g,K+1);
  // printpol(w);
  // printf("\n");
  // exit(1);
  // w=poly_mul(f,f);

  // 既約性判定のためのBen-Orアルゴリズム。拡大体にも対応している。デフォルトでGF(8192)
  // 既約多項式しか使わない。

  const int max_polynomial_attempts = 10000;
  bool found_valid_polynomial = false;
  vec pp = {0}, tt = {0};

  for (int attempt = 0; attempt < max_polynomial_attempts; attempt++)
  {
    memset(pp.x, 0, sizeof(pp.x));
    memset(tt.x, 0, sizeof(tt.x));
    for (i = 0; i < K; i++)
      pp.x[i] = rand() % N;
    mykey(tt.x, pp);
    tt.x[K] = 1;
    l = ben_or(tt);
    if (l < 0)
      continue;

    w = setpol(tt.x, K + 1);
    found_valid_polynomial = true;
    for (i = 0; i < N; i++)
    {
      ta[i] = trace(w, i);
      if (ta[i] == 0)
      {
        printf("trace 0 @ %d\n", i);
        found_valid_polynomial = false;
        break;
      }
    }
    if (found_valid_polynomial)
      break;
  }

  if (!found_valid_polynomial)
  {
    fprintf(stderr, "failed to find a valid Goppa polynomial after %d attempts\n",
            max_polynomial_attempts);
    return -1;
  }
  for (i = 0; i < N; i++)
  {
    tr[i] = oinv(ta[i]);
    // printf("%d,", tr[i]);
  }

  printpol(w);
  printf(" =irreducible\n");
  printsage(w);
  printf("\n");
  // j++;
  // if(j>100)
  // exit(1);

  // 多項式を固定したい場合コメントアウトする。
  /*
  memset(ta, 0, sizeof(ta));
  w = setpol(g, K + 1);
  printpol(w);
  //printf(" =poly\n");
  for (i = 0; i < N; i++)
  {    ta[i] = trace(w, i);
    if (ta[i] == 0)
    {
      printf("trace 0 @ %d\n", i);
      fail = 1;
      break;
    }
  }
  printpol(w);
  printf("\n");
  printsage(w);
  printf("\n");
  printf("sagemath で既約性を検査してください！\n");
  */

  van(powers);
  ogt(tt.x, polynomial_table);

  printf("\nすげ、オレもうイキそ・・・\n");

  for (j = 0; j < N; j++)
  {
    for (i = 0; i < K; i++)
    {
      ma[j][i] = gf[mlt(fg[powers[i][j]], tr[j])];
    }
  }

  unsigned short s;
  for (i = 0; i < K; i++)
  {
    for (j = 0; j < N; j++)
    {
      s = 0;

      for (k = 0; k < K; k++)
        s ^= gf[mlt(fg[polynomial_table[k][i]], fg[ma[j][k]])];
       printf("%d,",s);
      matrix[j][i] = s;
    }
  }
  // printf("\n");
  // exit(1);

  for (j = 0; j < N; j++)
  {
    for (i = 0; i < K; i++)
      printf("%d,", matrix[j][i]);
    printf("\n");
  }
  // exit(1);
  // wait();

  *result = w;
  return 0;
}

// Niederreiter暗号の公開鍵を作る
vec pubkeygen(unsigned short permutation[N],
              MTX *inverse_scrambler,
              RNG *rng,
              unsigned short syndrome_matrix[N][K])
{
  int i, j;
  FILE *fp;
  unsigned short field_matrix[N][K * E] = {0};
  MTX check_matrix = {0}, scrambler = {0}, public_matrix = {0};
  MTX public_check_matrix = {0};
  vec w = {0};
  vec v = {0};

  if (mkg(field_matrix, &w) != 0)
    exit(EXIT_FAILURE);
  printpol(w);
  printf(" ==goppa polynomial\n");

  v = w;
  //fp = fopen("sk.key", "wb");
  //fwrite(g, 2, K + 1, fp);
  //fclose(fp);
  printpol(w);
  printf("\n");
  printsage(w);
  printf("\n");
  printf("sagemath で既約性を検査してください！\n");

  bdet(field_matrix, &check_matrix);

  Pgen(permutation);
  fp = fopen("P.key", "w");
  fclose(fp);
  if (random_full_rank_binary_matrix(&scrambler, rng) != 0)
  {
    fprintf(stderr, "failed to generate a full-rank scrambler matrix\n");
    exit(EXIT_FAILURE);
  }
  if (is_reg(&scrambler, inverse_scrambler, rng) != 0)
  {
    fprintf(stderr, "failed to invert generated full-rank scrambler matrix\n");
    exit(EXIT_FAILURE);
  }

  /*
    fp = fopen("inv_S.key", "wb");
  for(i=0;i<K*E;i++){
    for(j=0;j<K*E;j++)
      n[j]=inv_S.x[i][j];
      fwrite(n,1,K*E,fp);
  }
  fclose(fp);
  */
  /*
    for (i = 0; i < K * E; i++)
    {
      for (j = 0; j < K; j++)
      {
        memset(v.x, 0, sizeof(v.x));
        for (k = 0; k < E; k++)
          v.x[k] = inv_S.x[i][j * E + k];
        n[j] = v2i(v);
        //printf("%d,", n[j]);
      }
      //printf("\n");
    //  fwrite(n, 2, K, fp);
    }
  */

  //  exit(1);
  public_check_matrix = mulmat(&scrambler, &check_matrix, 1);
  for (i = 0; i < K * E; i++)
  {
    for (j = 0; j < N; j++)
      public_matrix.x[j][i] = public_check_matrix.x[permutation[j]][i];
  }
  toByte(&public_matrix, syndrome_matrix);

  return v;
}

static vec dec(const MTX *inverse_scrambler, unsigned short ss[])
{
  int i, j;
  vec v = {0};
  vec s = {0};
  unsigned ch[K * E] = {0};
  unsigned char h2o[K * E] = {0};

  for (i = 0; i < K; i++)
  {
    v = i2v(ss[i]);
    for (j = 0; j < E; j++)
      ch[i * E + j] = v.x[j];
  }
  for (i = 0; i < K * E; i++)
    printf("%d", ch[i]);
  printf("\n");

  unsigned short uk[K] = {0};

  for (i = 0; i < K * E; i++)
  {
    for (j = 0; j < K * E; j++)
      h2o[i] ^= (ch[j] & inverse_scrambler->x[i][j]);
  }
  // for (i = 0; i < K * E; i++)
  // printf("%d,", h2o[i]);
  // printf("\n");

  for (i = 0; i < K; i++)
  {
    memset(v.x, 0, sizeof(v.x));
    for (j = 0; j < E; j++)
      v.x[j] = h2o[i * E + j];
    uk[i] = v2i(v);
  }
  for (i = 0; i < K; i++)
    printf("%d,", uk[i]);
  printf("\n");
  //    exit(1);
  s = setpol(uk, K);

  return s;
}

int ero2(vec v, const unsigned short permutation[N])
{
  int i, count = 0;
  unsigned short ya[N] = {0}, xa[N] = {0};

  for (i = 0; i < T * 2; i++)
  {
    if (i == 0)
    {
      xa[v.x[i]] = 1;
      // printf("error position=%d %d う\n", i, v.x[i]);
      count++;
    }
    if (i > 0 && v.x[i] > 0)
    {
      xa[v.x[i]] = 1;
      // printf("error position=%d %d お\n", i, v.x[i]);
      count++;
    }
    if (i > 0 && v.x[i] == 0)
    {
      printf("baka %d %d\n", i, v.x[i]);
      printf("v.x[K-1]=%d\n", v.x[K - 1]);
      /*
        printpol(w);
        printf(" ============goppa\n");
        printsage(w);
        printf(" ============sage\n");
        printsage(f);
        printf(" ============syn\n");
        printpol(f);
        printf(" ==========synd\n");
        printf("{");
        for (k = 0; k < N; k++)
        {
          if (z1[k] > 0)
            printf("%d,", z1[k]);
        }
        printf("};\n");
        //AA++;
        //wait();
        */
      break;
      //
      // exit (1);
    }
    /*
      for (k = 0; k < N; k++)
      {
        if (z1[k] > 0)
        {
          if (k != v.x[cnt])
          {
            printf("%d,%d\n", k, v.x[cnt]);
            printsage(w);
            printf(" ========w\n");
            break;
            //exit(1);
          }
          cnt++;
        }
      }
*/
  }
  for (i = 0; i < N; i++)
    ya[i] = xa[permutation[i]];
  for (i = 0; i < N; i++)
  {
    if (ya[i] > 0 && i == 0)
    {
      printf("error position=%d う\n", i);
    }
    else if (ya[i] > 0)
    {
      printf("error position=%d お\n", i);
    }
  }
  // exit(1);

  if (count == T * 2)
  {
    printf("err=%dっ!! \n", count);
  }
  if (count < T * 2)
  {
    printf("error is too few\n");

    // memcpy (zz, z1, sizeof (zz));
    /*
      printf("{");
      for (i = 0; i < D; i++)
        printf("%d,", z1[i]);
      printf("};\n");
      printpol(w);
      printf(" =========goppa\n");
      printsage(w);
      printf(" =========sage\n");
      printsage(f);
      printf(" =========syn\n");
      printpol(f);
      printf(" ==========synd\n");
      */
    printf("へげえええーっ\n");
    // exit(1);
    return -1;
    //exit(1);
  }

  return count;
}

int mkerr(unsigned short *z1, int num, RNG *rng)
{
  unsigned short positions[N];
  int available = 0;

  if (num < 0 || num > N)
  {
    fprintf(stderr, "invalid error weight: %d (expected 0..%d)\n", num, N);
    return -1;
  }

  for (int i = 0; i < N; i++)
  {
    if (z1[i] == 0)
      positions[available++] = (unsigned short)i;
  }
  if (num > available)
  {
    fprintf(stderr, "cannot add %d errors; only %d zero positions remain\n", num, available);
    return -1;
  }

  for (int i = 0; i < num; i++)
  {
    int selected = i + (int)(xor128(rng) % (unsigned long)(available - i));
    unsigned short position = positions[selected];
    positions[selected] = positions[i];
    positions[i] = position;
    z1[position] = 1;
    printf("l=%u %d\n", position, i);
  }

  for (int i = 0; i < N; i++)
  {
    if (z1[i] > 0)
      printf("%d=%d\n", i, z1[i]);
  }
  return 0;
}

vec zin(const MTX *inverse_scrambler,
        const unsigned short syndrome_matrix[N][K],
        unsigned short zz[])
{
  int i, j;
  vec s = {0};
  unsigned short ss[K] = {0};

  for (i = 0; i < N; i++)
  {
    if (zz[i] > 0)
    {
      for (j = 0; j < K; j++)
      {
        ss[j] ^= syndrome_matrix[i][j];
      }
      // printf("\n");
    }
  }

  // 暗号化されたシンドロームを復元する
  s = dec(inverse_scrambler, ss);
  for (j = 0; j < K; j++)
    printf("%d,", ss[j]);
  printf(" ==ss\n");
  // exit(1);

  return s;
}
