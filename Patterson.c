// Niederreiter cryptosystem with Patterson decoding.

#include "common.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#define REG 64
extern void makeS();

#include <time.h>
static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec + t.tv_nsec*1e-9; }

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

// OP型からベクトル型への変換
vec o2v(OP f)
{
  vec a = {0};
  int i;

  // #pragma omp parallel for
  for (i = 0; i < DEG; i++)
  {
    if (f.t[i].a > 0 && f.t[i].n < DEG)
      a.x[f.t[i].n] = f.t[i].a;
  }

  return a;
}

// ベクトル型からOP型への変換
OP v2o(vec a)
{
  int i, j = 0;
  OP f = {0};

  // #pragma omp parallel for
  for (i = 0; i < DEG; i++)
  {
    if (a.x[i] > 0)
    {
      f.t[j].n = i;
      f.t[j++].a = a.x[i];
    }
  }

  return f;
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

// OP型を正規化する
OP conv(OP f)
{
  int i;
  int previous_degree = -1;
  bool has_ended = false;

  for (i = 0; i < DEG; i++)
  {
    if (f.t[i].a == 0)
    {
      has_ended = true;
      if (f.t[i].n != 0)
        break;
      continue;
    }

    if (has_ended || f.t[i].n >= DEG || f.t[i].n <= previous_degree)
      break;

    previous_degree = f.t[i].n;
  }

  if (i == DEG)
    return f;

  vec v = {0};
  OP g = {0};

  v = o2v(f);
  g = v2o(v);

  return g;
}

// 項の数
int terms(OP f)
{
  int i, count = 0;

  for (i = 0; i < DEG; i++)
  {
    if (f.t[i].a > 0)
      count++;
  }

  return count;
}


// 零多項式なら -1、それ以外は最高次の添字
static inline int deg(const vec *a)
{
    for (int i = DEG - 1; i >= 0; i--)
        if (a->x[i] != 0)
            return i;
    return -1;
}


// 多項式の次数(default)
int vdeg(vec a)
{
  int i, n = 0, flg = 0;

  // #pragma omp parallel for
  for (i = DEG-1; i >=0; i--)
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


// 多項式の次数(degのOP型)
int odeg(OP f)
{
  int i, j = 0;

  if (f.t[0].a == 0)
    return 0;

  // k=terms(f);
  for (i = DEG-1; i >= 0; i--)
  {
    if (j < f.t[i].n && f.t[i].a > 0)
      j = f.t[i].n;
  }

  return j;
}


// 多項式を表示する（OP型）
void oprintpol(OP f)
{
  int i, n;

  f = conv(f);
  n = odeg(f);
  printf("n=%d\n", n);
  // printf("terms=%d\n", terms(f));
  printf("deg=%d\n", odeg(f));

  for (i = n; i > -1; i--)
  {
    if (f.t[i].a > 0)
      printf("%ux^%u+", f.t[i].a, f.t[i].n);
  }
}

void op_print_raw(const OP f)
{
  puts("op_print_raw:");
  for (int i = 0; i < DEG; i++)
  {
    if (f.t[i].a > 0)
      printf("[%d] %ux^%u\n", i, f.t[i].a, f.t[i].n);
  }
}

bool op_verify(const OP f)
{
  bool end = false;
  unsigned short n_max = 0;
  for (int i = 0; i < DEG; i++)
  {
    if (end && (f.t[i].n != 0 || f.t[i].a != 0))
    {
      op_print_raw(f);
      printf("found data after end: i=%d\n", i);
      print_trace();
      fflush(stdout);
      return false;
    }
    if (f.t[i].a == 0)
    {
      end = true;
      continue;
    }
    if (f.t[i].n + 1 <= n_max)
    {
      op_print_raw(f);
      printf("found invalid order: i=%d\n", i);
      print_trace();
      fflush(stdout);
      return false;
    }
    n_max = f.t[i].n + 1;
  }
  return true;
}

OP norm(OP f)
{
  OP h = {0};
  int i;
  for (i = 0; i < DEG; i++)
  {
    if (f.t[i].a > 0)
    {
      //      h.t[f.t[i].n].n=f.t[i].n;
      h.t[f.t[i].n].a = f.t[i].a;
    }
  }

  // exit(1);

  return h;
}

vec vadd(vec a, vec b)
{
  vec c = {0};
  size_t len = (size_t)((deg(&a) > deg(&b) ? deg(&a) : deg(&b)) + 1);

  xor_u16(c.x, a.x, len);
  xor_u16(c.x, b.x, len);

  return c;
}

vec vmul(vec a, vec b)
{
  int i, j, k, l;
  vec c = {0};

  k = deg(&a);
  l = deg(&b);

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

// 20200816:正規化したいところだがうまく行かない
// 多項式の足し算
OP oadd(OP f, OP g)
{
  vec a = {0}, b = {0};
  OP h;

  for (int i = 0; i < DEG; i++)
  {
    if (f.t[i].a > 0 && f.t[i].n < DEG)
      a.x[f.t[i].n] = f.t[i].a;
    if (g.t[i].a > 0 && g.t[i].n < DEG)
      b.x[g.t[i].n] = g.t[i].a;
  }

  for (int i = 0; i < DEG; i++)
    a.x[i] ^= b.x[i];

  h = v2o(a);
  assert(op_verify(h));
  return h;
}

// 項の順序を降順に揃える
OP sort(OP f)
{
  oterm o = {0};
  int i, j, k;

  k = terms(f);
  for (i = 0; i < k + 1; ++i)
  {
    for (j = i + 1; j < k + 1; ++j)
    {
      if (f.t[i].n > f.t[j].n)
      {
        o = f.t[i];
        f.t[i] = f.t[j];
        f.t[j] = o;
      }
    }
  }

  return f;
}

// リーディングタームを抽出(OP型）
oterm oLT(OP f)
{
  int i, k, j;
  oterm s = {0};

  k = terms(f);
  s = f.t[0];
  for (i = 0; i < k + 1; i++)
  {
    // printf("a=%d %d\n",f.t[i].a,f.t[i].n);
    if (f.t[i].a > 0)
    {
      printf("in LT=%d %d\n", s.a, s.n);
      for (j = i; j < k + 1; j++)
      {
        if (s.n < f.t[j].n)
        {
          s.n = f.t[j].n;
          s.a = f.t[j].a;
        }

        //  else{
        // t=s;
        // }
      }
    }
  }
  //  exit(1);

  return s;
}

// 多項式を足し算する（OP型）
OP add(OP f, OP g)
{
  //  vec a={0},b={0},c={0};
  unsigned long long int i, j, n1 = 0, n2 = 0, m1 = 0, count = 0;
  OP h = {0};
  oterm o1 = {0}, o2 = {0};

  n1 = terms(f);
  printf("n1=%llu\n", n1);
  n2 = terms(g);
  printf("n2=%llu\n", n2);
  if (n1 > n2)
  {
  }

  oprintpol(f);
  printf(" fff==============\n");
  oprintpol(g);
  printf(" ggg==============\n");
  o1 = oLT(f);
  o2 = oLT(g);
  printf("LTadd==%d %d\n", o1.n, o2.n);
  m1 = n1 + n2;
  printf("m1=%llu\n", m1);
  // exit(1);

  for (i = 0; i < n1 + 1; i++)
  {
    for (j = 0; j < n2 + 1; j++)
    {
      if (f.t[i].n == g.t[j].n && g.t[j].a > 0 && f.t[i].a > 0)
      {
        o1 = oLT(f);
        o2 = oLT(g);
        printf("LT==%d %d\n", o1.n, o2.n);
        printf("f.n==%u %u %llu %llu\n", f.t[i].n, g.t[j].n, i, j);
        f.t[i].a = 0;
        g.t[j].a = 0;
      }
    }
  }
  for (i = 0; i < n2 + 1; i++)
  {
    if (g.t[i].a > 0)
    {
      h.t[count++] = g.t[i];
      g.t[i].a = 0;
    }
  }
  for (i = 0; i < n1 + 1; i++)
  {
    if (f.t[i].a > 0)
    {
      h.t[count++] = f.t[i];
      f.t[i].a = 0;
    }
  }

  h = sort(h);
  h = conv(h);
  if (odeg(h) > 0)
    oprintpol(h);
  printf(" addh==============\n");
  //   exit(1);

  return h;
}

// 多項式を項ずつ掛ける
OP oterml(OP f, oterm t)
{
  f = conv(f);
  assert(op_verify(f));
  int i;
  OP h = {0};

  // f=conv(f);
  // k = deg (o2v(f));
  for (i = 0; i < DEG; i++)
  {
    h.t[i].n = f.t[i].n + t.n;
    h.t[i].a = gf[mlt(fg[f.t[i].a], fg[t.a])];
  }

  h = conv(h);
  assert(op_verify(h));
  return h;
}

// 多項式の掛け算
OP omul(OP f, OP g)
{
  f = conv(f);
  g = conv(g);
  assert(op_verify(f));
  assert(op_verify(g));
  OP h = v2o(fmul(o2v(f), o2v(g)));
  assert(op_verify(h));
  return h;
}

// リーディグタームを抽出(default)
oterm LT(OP f)
{
  int i;
  oterm t = {0};

  // k = deg (o2v (f));
  for (i = 0; i < DEG; i++)
  {
    // printf("a=%d %d\n",f.t[i].a,f.t[i].n);
    if (f.t[i].a > 0)
    {
      t.n = f.t[i].n;
      t.a = f.t[i].a;
    }
  }

  return t;
}

// 多項式の最後の項を抽出
oterm LT2(OP f)
{
  oterm t = {0};

  t.n = f.t[0].n;
  t.a = f.t[0].a;

  return t;
}

// 多項式を単行式で割る
oterm LTdiv(OP f, oterm t)
{
  oterm tt = {0}, s = {
                      0};

  tt = LT(f);
  if (tt.n < t.n)
  {
    s.n = 0;
    s.a = 0;
  }
  else if (tt.n == t.n)
  {
    s.n = 0;
    s.a = equ(t.a, tt.a);
  }
  else if (tt.n > t.n)
  {
    s.n = tt.n - t.n;
    s.a = equ(t.a, tt.a);
    // printf("%u\n",s.a);
  }
  else if (t.n == 0 && t.a > 0)
  {
    s.a = gf[mlt(fg[tt.a], oinv(t.a))];
    s.n = tt.n;
  }

  return s;
}

// モニック多項式にする
OP coeff(OP f, unsigned short d)
{
  int i, k;

  f = conv(f);
  k = odeg((f)) + 1;
  for (i = 0; i < k; i++)
    f.t[i].a = gf[mlt(fg[f.t[i].a], oinv(d))];

  return f;
}

// 多項式を表示する(default)
void printpol(vec a)
{
  int i, n;

  n = deg(&a);

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

// 多項式の商を取る
OP odiv(OP f, OP g)
{

  f = conv(f);
  g = conv(g);
  assert(op_verify(f));
  assert(op_verify(g));
  int i = 0;
  OP h = {0}, tt = {0};
  oterm b = {0}, c = {0};

  if (LT(f).n == 0 && LT(g).a == 0)
  {
    printf("baka^\n");
    // return f;
    exit(1);
  }
  if (LT(g).a == 0)
  {
    print_trace();
    exit(1);
  }
  if (LT(g).n == 0 && LT(g).a > 1)
    return coeff(f, LT(g).a);

  b = LT(g);
  if (b.a == 1 && b.n == 0)
    return f;
  if (b.a == 0 && b.n == 0)
  {
    printf("baka in odiv\n");
    exit(1);
  }
  if (odeg((f)) < odeg((g)))
  {
    return f;
    //  a=LT(f);
  }

  i = 0;
  while (LT(f).n > 0 && LT(g).n > 0)
  {
    c = LTdiv(f, b);
    assert(c.n < DEG);
    tt.t[i] = c;
    i++;

    h = oterml(g, c);

    f = oadd(f, h);
    if (odeg((f)) == 0 || odeg((g)) == 0)
    {
      // printf ("blake2\n");
      break;
    }

    if (c.n == 0)
      break;
  }

  // tt は逆順に入ってるので入れ替える
  OP ret = {0};
  int tt_terms = terms(tt);
  for (i = 0; i < tt_terms; i++)
  {
    ret.t[i] = tt.t[tt_terms - i - 1];
  }
  ret = conv(ret);
  assert(op_verify(ret));
  return ret;
}



// 多項式のべき乗
OP opow(OP f, int n)
{
  int i;
  OP g = {0};

  g = f;

  for (i = 1; i < n; i++)
    g = omul(g, f);

  return g;
}

// 多項式の剰余を取る
OP omod(OP f, OP g)
{
  OP h = {0};
  oterm b = {0}, c = {0};

  // n = LT(g).n;

  //  assert (("baka^\n", LT (f).n != 0));

  //  assert (("baka(A)\n", LT (g).n != 0));

  if (LT(f).n < LT(g).n)
  {
    //    exit(1);
    return f;
  }

  // printf ("in omod\n");
  // exit(1);

  b = LT(g);
  // printpol(o2v(g));
  // printf("d %d\n",LT(g).n);

  // assert(("double baka\n", b.a != 0 && b.n != 0));
  while (LT(f).n > 0 && LT(g).n > 0)
  {

    c = LTdiv(f, b);
    h = oterml(g, c);
    f = oadd(f, h);
    if (odeg((f)) == 0 || odeg((g)) == 0)
    {
      //      printf("blake1\n");
      break;
    }

    if (c.n == 0 || b.n == 0)
      break;
  }

  return f;
}

// 多項式のべき乗余
OP opowmod(OP f, OP mod, int n)
{
  int i = 1;

  // printpol(o2v(mod));
  // printf(" =mod %d\n",LT(mod).n);
  // 繰り返し２乗法
  while (i < n + 1)
  {
    f = omod(omul(f, f), mod);
    i++;
  }

  return f;
}

// 多項式の代入値
static unsigned short trace(vec f, unsigned short x)
{
  unsigned short u = 0;
  int degree = deg(&f);

  for (int i = 0; i <= degree; i++)
    if (f.x[i] != 0)
      u ^= gf[mlt(fg[f.x[i]], mltn(i, fg[x]))];

  return u;
}

// invert of polynomial
OP inv(OP a, OP n)
{
  OP d = {0}, x = {0}, s = {0}, q = {0}, r = {0}, t = {0}, u = {0}, v = {0}, w = {0}, tt = {0}, gcd = {0}, tmp = {0};
  vec aa = {0}, bb = {0}, cc = {0};

  if (odeg((a)) > odeg((n)))
  {
    tmp = a;
    a = n;
    n = tmp;
    printf("baka_i\n");
    // exit (1);
  }
  if (LT(a).a == 0)
  {
    printf(" a ga 0\n");
    exit(1);
  }

  tt = n;

  d = n;
  x.t[0].a = 0;
  x.t[0].n = 0;
  s.t[0].a = 1;
  s.t[0].n = 0;
  while (odeg((a)) > 0)
  {
    if (odeg((a)) > 0)
      r = omod(d, a);
    if (LT(a).a == 0)
      break;
    if (LT(a).a > 0)
      q = odiv(d, a);

    d = a;
    a = r;
    t = oadd(x, omul(q, s));
    ////printpol (o2v (a));
    // printf ("\nin roop a==================%d\n", odeg ((a)));
    // printf ("\n");

    x = s;
    s = t;
  }
  // exit(1);
  //  if(LT(a).a>0){
  d = a;
  a = r;
  ////printpol (o2v (a));
  // printf ("\nin roop a|==================%d\n", odeg ((a)));
  // printf ("\n");

  x = s;
  s = t;

  ////printpol (o2v (d));
  // printf ("\nout1================\n");
  gcd = d; // $\gcd(a, n)$
  printpol(o2v(gcd));
  printf(" =========gcd\n");
  // exit(1);
  // printf ("\n");
  ////printpol (o2v (n));
  // printf ("\n");
  // printf ("out2===============\n");

  printf("before odiv\n");
  // w=tt;

  ////printpol (o2v (w));
  // printf ("\nw=======%d %d\n", b.a, b.n);
  // w=tt;
  aa = o2v(x);
  bb = o2v(n);
  cc = vadd(aa, bb);
  v = v2o(cc);
  w = tt;
  if (LT(v).n > 0 && LT(w).n > 0)
  {
    u = omod(v, w);
  }
  else
  {
    // printpol (o2v (v));
    printf(" v===========\n");
    // printpol (o2v (x));
    printf(" x==0?\n");
    // printpol (o2v (n));
    printf(" n==0?\n");

    exit(1);
  }
  // caution !!
  if (LT(u).a > 0 && LT(d).a > 0)
  {
    u = odiv(u, d);
  }

  if (LT(u).a == 0 || LT(d).a == 0)
  {
    printf("inv div u or d==0\n");
    // exit(1);
  }
  // u=coeff(u,d.t[0].a);
  ////printpol (o2v (u));
  // printf ("\nu==================\n");
  if (LT(u).a == 0)
  {
    printf("no return at u==0\n");
    exit(1);
  }

  return u;
}

unsigned short
v2a(oterm a)
{
  int j;

  if (a.a == 0)
    return 0;

  // printf("aa=%d\n",a.a);
  for (j = 0; j < N; j++)
  {
    if (gf[j] == a.a && a.a > 0)
    {
      // printf("j==%d\n",j);
      return j - 1;
    }
  }
  fprintf(stderr, "Field element %u is not present in the field table.\n", a.a);
  exit(EXIT_FAILURE);
}

void printsage(vec a)
{
  int i, j;
  oterm b;

  printf("poly=");
  for (i = 0; i < DEG; i++)
  {
    if (a.x[i] > 0)
    {
      b.a = a.x[i];
      b.n = i;
      j = v2a(b);
      // printf("%d,==ba\n",b.a);
      // printf ("X**%d+", i); //for GF2
      printf("B('a^%d')*X**%d+", j, i); // for GF(2^m)
    }
  }
}


// vec 多項式の剰余
vec vmod2(vec f, vec g)
{
    vec r = f;
    int dg = deg(&g);

    // ゼロ多項式で割ることはできない
    bool gzero = true;
    for (int i = 0; i < DEG; i++)
        if (g.x[i] != 0)
            gzero = false;

    if (gzero)
        exit(1);

    while (true)
    {
        int dr = deg(&r);

        bool rzero = true;
        for (int i = 0; i < DEG; i++)
            if (r.x[i] != 0)
                rzero = false;

        if (rzero || dr < dg)
            break;

        int shift = dr - dg;

        // 最高次係数の比（対数表現で除算）
        int q = gf[mlt(fg[r.x[dr]],
                    oinv(g.x[dg]))];

        // r ^= g * q * x^shift
        for (int j = 0; j <= dg; j++)
        {
            if (g.x[j] == 0)
                continue;

            int k = j + shift;
            int product = gf[mlt(fg[g.x[j]], fg[q])];

            r.x[k] = fg[r.x[k]] ^ product;

            // 通常表現から対数表現へ戻す
            r.x[k] = (r.x[k] == 0) ? 0 : fg[r.x[k]];
        }

        // 最高次項が消えていないなら異常
        if (deg(&r) >= dr)
        {
            printf("vmod: degree did not decrease\n");
            exit(1);
        }
    }

    return r;
}


// vec 多項式の剰余 (f mod g)
// 係数は通常表現（GF(2^E) の元そのもの）で保持する。
vec vmod(vec f, vec g)
{
    vec r = f;
    int i;

    // g の次数を求める（零多項式なら -1）
    int dg = -1;
    for (i = DEG - 1; i >= 0; i--)
    {
        if (g.x[i] != 0)
        {
            dg = i;
            break;
        }
    }
    if (dg < 0)
    {
        fprintf(stderr, "vmod: division by zero polynomial\n");
        exit(1);
    }

    // g の最高次係数の逆元（元のコードと同じく oinv の戻り値を mlt にそのまま渡す）
    const int inv_lead = oinv(g.x[dg]);

    // r の次数を求める（零多項式なら -1）
    int dr = -1;
    for (i = DEG - 1; i >= 0; i--)
    {
        if (r.x[i] != 0)
        {
            dr = i;
            break;
        }
    }

    // 零多項式（dr = -1）や dr < dg なら、そのまま剰余
    while (dr >= dg)
    {
        int shift = dr - dg;

        // 最高次係数の比 q = r[dr] / g[dg]
        int q = gf[mlt(fg[r.x[dr]], inv_lead)];

        // r ^= q * x^shift * g
        for (int j = 0; j <= dg; j++)
        {
            if (g.x[j] == 0)
                continue;
            r.x[j + shift] ^= gf[mlt(fg[g.x[j]], fg[q])];
        }

        // 最高次項が消えていないなら異常
        if (r.x[dr] != 0)
        {
            fprintf(stderr, "vmod: leading term was not eliminated\n");
            exit(1);
        }

        // 次の最高次へ（零多項式になれば dr = -1 でループ終了）
        while (dr >= 0 && r.x[dr] == 0)
            dr--;
    }

    return r;
}


vec vgcd2(vec xx, vec yy)
{
    vec tt = {0};

    while (1)
    {
        tt = o2v(omod(v2o(xx), v2o(yy)));
        xx = yy;
        yy = tt;
        if(LT(v2o(yy)).a==0)
        break;
    }

    return xx;
}

// gcd for patterson
OP zgcd(OP a, OP n, int P)
{
  OP d = {0}, x = {0}, s = {0}, q = {0}, r = {0}, t = {0}, rt = {0};

  if (odeg(a) > odeg(n))
  {
    rt = a;
    a = n;
    n = rt;
    printf("big is good\n");
    // exit (1);
  }
  if (LT(a).a == 0)
  {
    printf(" a ga 0\n");
    exit(1);
  }

  d = n;
  x.t[0].a = 0;
  x.t[0].n = 0;
  s.t[0].a = 1;
  s.t[0].n = 0;
  while (LT(a).n > P)
  {

    r = omod(d, a);
    q = odiv(d, a);

    d = a;
    a = r;
    t = oadd(x, omul(q, s));

    x = s;
    s = t;
  }

  d = a;
  a = r;

  x = s;
  s = t;

  printpol(o2v(x));
  printf(" =======x\n");
  printpol(o2v(a));
  printf(" =======a\n");
  printpol(o2v(s));
  printf(" =======s\n");
  printpol(o2v(r));
  printf(" =======r\n");

  return x;
}

// gcd
OP gcd(OP xx, OP yy)
{
  OP tt = {0}, tmp, h = {0};

  h.t[0].a = 1;
  h.t[0].n = 0;
  vec vtmp=o2v(xx);
  vec temp=o2v(yy);
  if (deg(&vtmp) < deg(&temp))
  {
    tmp = xx;
    xx = yy;
    yy = tmp;
  }
  tt = omod(xx, yy);
  while (odeg(tt) > 0)
  {
    xx = yy;
    yy = tt;
    if (odeg(yy) > 0)
      tt = omod(xx, yy);
    if (LT(tt).a == 0)
      return yy;
  }
  if (LT(tt).a == 0)
  {
    return yy;
  }
  else
  {
    return h;
  }
  //  return yy;
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

  for (i = 0; i < deg(&v) + 1; i++)
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

static OP bib(int i, OP d, const OP ww[T])
{
  int j;

  OP t = d;

  // #pragma omp parallel for
  for (j = 0; j < T; j++)
  {
    // #pragma omp critical
    if (i != j)
    {
      t = omul(t, ww[j]);
    }
  }

  return t;
}


// chen探索
static vec chen(vec f)
{
  vec e = {0};
  int count = 0, n, x = 0;
  unsigned short z;

  n = deg(&f);
  // exit(1);
  // #pragma omp parallel for private(i)
  for (x = 0; x < N; x++)
  {
    z = 0;
    // #pragma omp parallel for reduction (^:z)
    for (int i = 0; i <= n; i++)
      if (f.x[i] != 0)
        z ^= gf[mlt(mltn(i, fg[x]), fg[f.x[i]])];
    if (z == 0)
    {
      if (count < DEG)
      {
        e.x[count++] = x;
        printf("x=%d\n", x);
      }
    }
  }
  printpol(e);
  printf(" ==eee!\n");
  // exit(1);

  return e;
}


vec vgcd(vec xx, vec yy)
{
    vec tt = {0};

    while (1) {
        // yy がゼロ多項式なら xx が gcd
        int zero = 1;
        for (int i = 0; i < DEG; i++) {
            if (yy.x[i] != 0) {
                zero = 0;
                break;
            }
        }

        if (zero)
            break;

        tt = o2v(omod(v2o(xx), v2o(yy)));
        xx = yy;
        yy = tt;
    }

    return xx;
}


vec derivative(vec f) {
    vec df={0};
    int y = deg(&f);


    if(y>0){
    for (int i = 1; i <= y; i += 2) {
        df.x[i - 1] = f.x[i];
    }
  }  

    return df;
}

vec diff(vec a){
  vec b={0};

  b=derivative(a);

  b=(vgcd((a),(b)));
  

  return b;
}


// GF(2^m) then set m in this function.
int ben_or(OP f)
{
  int i, n;
  OP d = {0};
  OP *ww = calloc(T, sizeof(*ww));
  if (ww == NULL)
  {
    fprintf(stderr, "ben_or: failed to allocate polynomial workspace\n");
    return -1;
  }
  vec tmp=o2v(f);
  n = deg(&tmp);

  if (n == 0)
  {
    free(ww);
    return -1;
  }

#pragma omp parallel num_threads(omp_get_max_threads())
  {
#pragma omp for schedule(static)
    for (i = 0; i < n / 2 + 1; i++)
      (void)bib(i, d, ww);
  }

  free(ww);
  return 0;
}

// 配列の値を係数として多項式に設定する
OP setpol(unsigned short f[], int n)
{
  return v2o(Setvec(f, n));
}


/* ビット位置 a と b を交換する */
static unsigned short swap_bits(unsigned short v, unsigned int a, unsigned int b)
{
    unsigned int x = ((v >> a) ^ (v >> b)) & 1u;
    return (unsigned short)(v ^ ((x << a) | (x << b)));
}


// バイナリ型パリティチェック行列を生成する
static void bdet(const unsigned short matrix[N][K], MTX *check_matrix)
{
  int i, j, k, l;

  if (check_matrix->x == NULL || check_matrix->rows != N ||
      check_matrix->cols != (size_t)(K * E))
  {
    fprintf(stderr, "invalid check-matrix storage\n");
    exit(EXIT_FAILURE);
  }
  
  #pragma omp parallel for schedule(static)
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < K; j++)
    {
      l = matrix[i][j];
      for (k = 0; k < E; k++)
      {
        check_matrix->x[i][j * E + k] = l & 1;
        l = (l >> 1);
      }
    }
  }

}

static void toByte(const MTX *SH,
                   const unsigned short permutation[N],
                   unsigned short syndrome_matrix[N][K])
{
  int i, j, k;

  if (SH == NULL || SH->x == NULL || SH->rows < N ||
      SH->cols < (size_t)(K * E) ||
      permutation == NULL)
  {
    fprintf(stderr, "toByte: invalid matrix or permutation\n");
    exit(EXIT_FAILURE);
  }

  for (i = 0; i < N; i++)
  {
    if (permutation[i] >= SH->rows)
    {
      fprintf(stderr, "toByte: permutation index is outside the matrix\n");
      exit(EXIT_FAILURE);
    }
    // #pragma omp parallel for
    for (j = 0; j < K; j++)
    {
    unsigned short x=0;
        k = j * E;
        x=0;
        for (int l = 0; l < E; l++){
          x ^= (SH->x[permutation[i]][k + l] << l);
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


// sqrt(x) mod w を求める (w = g0^2 + x*g1^2 のとき sqrt(x) = g0/g1)
// 成功で0。鍵ごとに1回だけ呼ぶ。
int osqrt_x(OP w, OP *sx)
{
  OP h = {0}, r = {0}, s = {0}, ww = {0}, tmp = {0}, xx = {0};
  int i, j = 0, jj = 0;
  vec vtmp = o2v(w);
  int k = deg(&vtmp);

  for (i = 0; i < k + 1; i++)
  {
    if (w.t[i].a == 0) continue;
    if (w.t[i].n % 2 == 0)
    {
      h.t[j].a = gf[isqrt(w.t[i].a)];
      h.t[j++].n = w.t[i].n / 2;
    }
    else
    {
      r.t[jj].a = gf[isqrt(w.t[i].a)];
      r.t[jj++].n = (w.t[i].n - 1) / 2;
    }
  }
  if (jj == 0) return 1;            // 奇数次の項がない: g が平方式で square-free でない

  if (LT(r).n > 0)
  {
    s = inv(r, w);
    tmp = omod(omul(s, r), w);
    if (odeg(tmp) != 0 || LT(tmp).a != 1) return 2;   // 逆元が取れない
    ww = omod(omul(h, s), w);
  }
  else
  {
    // g1 が定数の場合 (元のコードの分岐をそのまま)
    s.t[0].a = gf[isqrt(LT(r).a)];
    s.t[0].n = 0;
    ww = coeff(h, s.t[0].a);
  }

  // 検証: ww^2 mod w が x と完全に一致するか
  xx.t[0].n = 1;
  xx.t[0].a = 1;
  tmp = oadd(omod(omul(ww, ww), w), xx);
  if (LT(tmp).a != 0) return 3;

  *sx = ww;
  return 0;
}

// f の平方根 mod w。sx は osqrt_x で作ったもの
OP osqrt_f(OP f, OP w, OP sx)
{
  OP even = {0}, odd = {0};
  int i, j = 0, jj = 0;
  vec vtmp = o2v(f);
  int k = deg(&vtmp);

  for (i = 0; i < k + 1; i++)
  {
    if (f.t[i].a == 0) continue;
    if (f.t[i].n % 2 == 0)
    {
      even.t[j].n = f.t[i].n / 2;
      even.t[j++].a = gf[isqrt(f.t[i].a)];
    }
    else
    {
      odd.t[jj].n = (f.t[i].n - 1) / 2;
      odd.t[jj++].a = gf[isqrt(f.t[i].a)];
    }
  }
  return omod(oadd(even, omul(sx, odd)), w);
}

// 既存の呼び出し用 (互換)。復号側は sx を使い回すように直すのが理想
// 多項式の平方を計算する
OP osqrt(OP f, OP w)
{
  OP sx = {0};
  if (osqrt_x(w, &sx) != 0) exit(1);
  return osqrt_f(f, w, sx);
}

// e: 本当の誤りベクトル、d: 復号結果 (どちらも長さ n, 誤り位置が非0)
// 戻り値: 本当の誤り T 個のうち、正しく見つかった個数
// *T_out: 本当の誤りの数, *extra: 誤って余分に出した位置の数
int count_correct_old(const vec *e, const vec *d, int n, int *T_out, int *extra)
{
    int t = 0, hit = 0, ex = 0;

    for (int i = 0; i < n; i++)
    {
        int ei = (e->x[i] != 0);
        int di = (d->x[i] != 0);

        if (ei) t++;
        if (ei && di) hit++;       // 正しく見つけた
        if (!ei && di) ex++;       // 誤って余分に出した
    }

    if (T_out) *T_out = t;
    if (extra) *extra = ex;
    return hit;
}

vec paloma(OP s, OP f)
{
    vec x = {0}, v = {0}, zero = {0};
    x.x[1] = 1;  v.x[0] = 1;
    OP one = v2o(v);

    if (odeg(s) < 0) return zero;          // 誤りなし

    OP ss = omod(oadd(one, omul(s, v2o(x))), f);

    OP g1 = gcd(s, f);
    OP g2 = gcd(ss, f);
    OP g12 = odiv(f, omul(g1, g2));

    OP s1 = odiv(s, g1);                   // 常に割る
    OP s2 = odiv(ss, g2);

    OP denom = omul(g2, s1);
    if (odeg(gcd(denom, g12)) > 0) {       // f が square-free なら起きないはず
        printf("denominator not invertible mod g12\n");
        exit(1);
    }

    OP u  = omod(omul(omul(g1, s2), inv(denom, g12)), g12);
    OP sp = osqrt(u, g12);

    /* 拡張ユークリッド: rb ≡ tb*sp (mod g12) を保つ */
    OP ra = g12, rb = sp;
    OP ta = {0}, tb = one;
    int Ba = K / 2 - odeg(g2);             // a' の次数の上限

    while (odeg(rb) > Ba) {                // 判定が先
        OP q = odiv(ra, rb);
        OP r = omod(ra, rb);
        OP t = oadd(ta, omul(q, tb));
        ra = rb;  rb = r;
        ta = tb;  tb = t;
    }
    OP a0 = rb, b0 = tb;

    OP a = omul(a0, g2);
    OP b = omul(b0, g1);
    OP sigma = oadd(omul(a, a), omul(omul(b, b), v2o(x)));

    int yt = odeg(sigma);                   // 誤りの個数 (K に固定しない)
    if (yt < 1 || yt > K) { printf("decode failed\n"); exit(1); }

    return chen(o2v(sigma));
}


// パターソンアルゴリズムでバイナリGoppa符号を復号する
vec patterson(vec w_vec, vec f_vec)
{
  OP w = v2o(w_vec), f = v2o(f_vec);
  OP g1 = {0}, ll = {0};
  int flg2 = 0;
  OP tt = {0}, ff = {0}, h = {0};
  vec v;
  OP r2 = {0}, b2 = {
                   0};
  // unsigned short g[K+1]={2,2,12,1,2,8,4,13,5,10,8,2,15,10,7,3,5};

  tt.t[0].n = 1;
  tt.t[0].a = 1;

  ff = inv(f, w);
  //  //printpol (o2v (ff));
  b2 = omod(omul(ff, f), w);
  if (odeg((b2)) > 0)
  {
    printf("逆元が計算できません。error\n");
    // wait ();
    exit(1);
  }
  printf("locater==========\n");
  // exit(1);
  r2 = oadd(ff, tt);
  printpol(o2v(r2));
  printf(" h+x==============\n");
  // wait();
  //   exit(1);
  g1 = osqrt(r2, w);
  printpol(o2v(g1));
  printf(" sqrt(h+x)=======\n");
  b2 = omod(omul(g1, g1), w);
  printpol(o2v(b2));
  printf(" g1*g1%%w===========\n");
  vec tmp=o2v(b2);
  if (deg(&tmp) == 0)
  {
    printf("b2 is 0\n");
    exit(1);
  }
  if (LT2(b2).a != LT2(r2).a)
  {
    printpol(o2v(w));
    printf(" w============\n");
    printpol(o2v(r2));
    printf(" r2============\n");
    printpol(o2v(g1));
    printf(" g1============\n");
    printf(" g1は平方ではありません。error");
    // wait ();
    exit(1);
  }
  printpol(o2v(w));
  printf(" ========goppa\n");
  printpol(o2v(g1));
  printf(" g1!=========\n");
  if (LT(g1).n == 0 && LT(g1).a == 0)
  {
    // printpol (o2v (w));
    printf(" badkey=========\n");
    // printpol (o2v (g1));
    printf(" g1============\n");
    printf("平方が０になりました。 error\n");
    // wait ();
    exit(1);
  }
  h = gcd(f, w);
  if (odeg((h)) > 0)
  {
    printsage(o2v(w));
    printf(" =goppa\n");
    printsage(o2v(f));
    printf(" =syn\n");
    printpol(o2v(h));
    printf(" =gcd\n");

    // printf("f-ben_or=%d\n",ben_or(f));
    // printf("w-ben_or=%d\n",ben_or(w));
    printf(" s,wは互いに素じゃありません。\n");
    flg2 = 1;
    // wait ();
    exit(1);
    // goto label;
  }

  // exit(1);
  OP ppo = {0};

  // hh = xgcd(w, g1);
  ppo = zgcd(w, g1, T);
  printpol(o2v(ppo));
  printf(" =========ppo\n");
  printpol(o2v(ppo));
  printf(" =========ppo\n");
  // if(LT(ppo).a!=LT(ppo).a)
  // exit(1);
  printpol(o2v(ppo));
  printf(" =========ppo\n");
  ff = omod(omul(ppo, g1), w);
  printpol(o2v(ppo));
  printf(" alpha!=========\n");
  printpol(o2v(ff));
  printf(" beta!=========\n");
  printpol(o2v(w));
  printf(" goppa=========\n");
  printpol(o2v(f));
  printf(" syn=========\n");

  printpol(o2v(ppo));
  printf(" alpha!=========\n");
  // hh = xgcd(w, g1);
  ppo = zgcd(w, g1, T);
  //  h=zgcd(w,g1);
  ff = omod(omul(ppo, g1), w);
  ll = oadd(omul(ff, ff), omul(tt, omul(ppo, ppo)));
  v = chen(o2v(ll));
  if (v.x[K - 1] > 0)
  {
    // wait ();
    return v;
  }

  ppo = zgcd(w, g1, T);
  ff = omod(omul(ppo, g1), w);
  if (odeg((ff)) != K / 2)
  {

    printf("\nbefore h.d\n");
    ff = omod(omul(ppo, g1), w);
    printpol(o2v(ff));
    printf(" ==========beta!\n");
    printpol(o2v(ppo));
    printf(" alpha!=========\n");
    ll = oadd(omul(ff, ff), omul(tt, omul(ppo, ppo)));
    v = chen(o2v(ll));
    if (v.x[K - 1] > 0)
    {
      wait();
      return v;
    }
  }

  if (odeg((ff)) == 1)
  {
    ll = oadd(omul(ff, ff), omul(tt, omul(ppo, ppo))); // ff;
    printf("deg==1\n");
    exit(1);
  }

  printf("あっ、でる・・・！\n");
  printpol(o2v(ll));
  printf(" ll=================\n");
  printpol(o2v(f));
  printf(" syn=================\n");
  v = chen(o2v(ll));

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
  MTX a = {0};
  vec v = {0};

  if (input == NULL || input->rows < (size_t)K ||
      input->cols < (size_t)(K + 1) || matrix_copy(&a, input) != 0)
  {
    fprintf(stderr, "renritu: invalid or unavailable matrix storage\n");
    return v;
  }

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

  vec pol = Setvec(x.x, K + 1);
  printpol(pol);
  printf(" ==key\n");
  for (i = 0; i < N; i++)
  {
    if (trace(pol, i) == 0)
      printf("%d i=%d\n", i, fg[i] - 1);
  }

  matrix_free(&a);
  return v;
}

static inline unsigned short eval(const unsigned short *c, unsigned short lx)
{
  unsigned short q = c[K];
  for (int i = K - 1; i >= 0; i--)
    q = gf[mlt(fg[q], lx)] ^ c[i];
  return q;
}


/* input: f, element in GF((2^m)^t) */
/* output: out, minimal polynomial of f */
/* return: 0 for success and -1 for failure */
int mykey(unsigned short *out, vec x)
{
  unsigned short (*mat)[K] = calloc((size_t)K + 1, sizeof(*mat));
  MTX a = {0};
  int i, j;

  if (out == NULL || mat == NULL || matrix_init(&a, K, K + 1) != 0)
  {
    free(mat);
    matrix_free(&a);
    fprintf(stderr, "mykey: failed to allocate matrix workspace\n");
    return -1;
  }

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
    }
  }
  // exit(1);

  vec v = {0};
  v = renritu(&a);

  for (i = 0; i < K; i++)
  {
    out[i] = v.x[i];
    printf("%d,", out[i]);
  }
  printf("\n");
  free(mat);
  matrix_free(&a);
  return 0;
}

// matrix[j][i] = (Σ_{m=0}^{i} g[m]·α_j^(i-m)) / g(α_j),  α_j = j
static void matrix_horner(const vec *polynomial,
                          const unsigned short *tr,   // tr[j] = oinv(g(j))
                          unsigned short matrix[N][K])
{
#pragma omp parallel for schedule(static)
  for (int j = 0; j < N; j++)
  {
    const unsigned short alpha = fg[j];
    unsigned short q = polynomial->x[0];
    matrix[j][0] = gf[mlt(fg[q], tr[j])];
    for (int i = 1; i < K; i++)
    {
      q = gf[mlt(alpha, fg[q])] ^ polynomial->x[i];
      matrix[j][i] = gf[mlt(fg[q], tr[j])];
    }
  }
}


int mkg(unsigned short matrix[N][K], vec *result)
{
  OP w = {0};
  vec pp = {0}, polynomial = {0};
  unsigned short *tr = calloc(N, sizeof(*tr));
  unsigned short *ta = calloc(N, sizeof(*ta));
  if (!tr || !ta) { free(tr); free(ta); return -1; }

bool found = false;
for (int attempt = 0; attempt < 100 && !found; attempt++)
{
  memset(pp.x, 0, sizeof(pp.x));
  for (int i = 0; i < K; i++) pp.x[i] = rand() % N;
  pp.x[K] = 1;
  w = setpol(pp.x, K + 1);
  polynomial = o2v(w);

  int has_root = 0;   /* 試行ごとに0に戻す */
#pragma omp parallel for schedule(static)
  for (int i = 0; i < N; i++)
  {
    int skip;
#pragma omp atomic read
    skip = has_root;
    if (skip) continue;
    ta[i] = i ? eval(polynomial.x, fg[i]) : polynomial.x[0];
    if (ta[i] == 0)
    {
#pragma omp atomic write
      has_root = 1;
    }
  }
  found = !has_root;
}
  if (!found) { free(tr); free(ta); return -1; }

  for (int i = 0; i < N; i++) tr[i] = oinv(ta[i]);
  matrix_horner(&pp, tr, matrix);
  *result = o2v(w);
  free(tr); free(ta);
  return 0;
}


// Niederreiter暗号の公開鍵を作る
vec pubkeygen(unsigned short permutation[N],
              MTX *inverse_scrambler,
              RNG *rng,
              unsigned short syndrome_matrix[N][K])
{
  FILE *fp;
  unsigned short (*field_matrix)[K] = calloc(N, sizeof(*field_matrix));
  MTX check_matrix = {0}, scrambler = {0}, public_check_matrix = {0};
  vec w = {0};
  double t0 = now();

  if (field_matrix == NULL ||
      //matrix_init(&check_matrix, N, K * E) != 0 ||
      matrix_init(&scrambler, K * E, K * E) != 0 ||
      matrix_init(inverse_scrambler, K * E, K * E) != 0)
  {
    fprintf(stderr, "pubkeygen: failed to allocate key-generation matrices\n");
    free(field_matrix);
    matrix_free(&check_matrix);
    matrix_free(&scrambler);
    matrix_free(&public_check_matrix);
    matrix_free(inverse_scrambler);
    exit(EXIT_FAILURE);
  }

  
  if (mkg(field_matrix, &w) != 0)
  {
    free(field_matrix);
    matrix_free(&check_matrix);
    matrix_free(&scrambler);
    matrix_free(&public_check_matrix);
    matrix_free(inverse_scrambler);
    exit(EXIT_FAILURE);
  }
  printpol(w);
  printf(" ==goppa polynomial\n");
  printf("Uh!-mkg: %.3f s\n", now() - t0); t0 = now();            // ← ② mkg の if ブロックの閉じ括弧の後

  //fp = fopen("sk.key", "wb");
  //fwrite(g, 2, K + 1, fp);
  //fclose(fp);
  printpol(w);
  printf("\n");
  printsage(w);
  printf("\n");
  printf("sagemath で既約性を検査してください！\n");

  //bdet(field_matrix, &check_matrix);
  //free(field_matrix);

  Pgen(permutation);
  //fp = fopen("P.key", "w");
  //fclose(fp);
  printf("Uh!-Pgen: %.3f s\n", now() - t0); t0 = now();
  /*
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
  */
  // random_full_rank_binary_matrix と is_reg の後
  printf("Uh!-scrambler: %.3f s\n", now() - t0); t0 = now();
  
  
  if (random_full_rank_binary_matrix_inv(&scrambler, inverse_scrambler, rng) != 0)
  {
    fprintf(stderr, "failed to generate the scrambler matrix and its inverse\n");
    exit(EXIT_FAILURE);
  }
  
  printf("Uh!-scrambler+inverse: %.3f s\n", now() - t0); t0 = now();

  //exit(1);

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
  if (pubkey_pack(&scrambler, field_matrix, permutation, syndrome_matrix) != 0)
  {
    fprintf(stderr, "pubkeygen: failed to build the public key\n");
    free(field_matrix);
    matrix_free(&check_matrix);
    matrix_free(&scrambler);
    matrix_free(&public_check_matrix);
    matrix_free(inverse_scrambler);
    exit(EXIT_FAILURE);
  }
  // pubkey_pack の後
  printf("Uh!-pubkey_pack: %.3f s\n", now() - t0);
  
  // 一致確認(OK と出たら、このブロックは削除してよい)
  /*
  {
    MTX cm = {0}, pc = {0};
    unsigned short (*ref)[K] = calloc(N, sizeof *ref);
    if (ref && matrix_init(&cm, N, K * E) == 0) {
      //bdet(field_matrix, &cm);
      mulmat(&scrambler, &cm, 1, &pc);
      //toByte(&pc, permutation, ref);
      printf("pubkey check: %s\n",
             memcmp(ref, syndrome_matrix, (size_t)N * sizeof *ref) == 0 ? "OK" : "MISMATCH");
    }
    matrix_free(&cm); matrix_free(&pc); free(ref);
  }
*/

  free(field_matrix);

  matrix_free(&check_matrix);
  matrix_free(&scrambler);
  matrix_free(&public_check_matrix);
  return w;
}


static vec dec(const MTX *inverse_scrambler, unsigned short ss[])
{
  int i, j;
  vec v = {0};
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
  return Setvec(uk, K);
}

int ero2(vec v, const unsigned short permutation[N],unsigned short *ya)
{
  int i, count = 0;
  unsigned short xa[N] = {0};

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
        printpol(o2v(w));
        printf(" ============goppa\n");
        printsage(o2v(w));
        printf(" ============sage\n");
        printsage(o2v(f));
        printf(" ============syn\n");
        printpol(o2v(f));
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
            printsage(o2v(w));
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
