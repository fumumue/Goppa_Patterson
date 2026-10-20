// Niederreiter cryptosystem with Patterson decoding.

#include "common.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#define REG 64
extern unsigned long xor128(void);
extern int mlt(int x, int y);
extern int mltn(int n, int a);
extern void makeS();

static unsigned short c[E * K + 1] = {0};
static unsigned short mat[N][K * E] = {0};
static unsigned short ma[N][K * E] = {0};
static unsigned short bm[N][K * E] = {0};
static unsigned short bm2[N][K * E] = {0};

// シンドロームのコピー
unsigned short P[N] = {0};
unsigned short inv_P[N] = {0};
unsigned short uu = 0;

// Goppa多項式
static unsigned short g[K + 1] = {1, 0, 0, 0, 1, 0, 1};
//{1,0,1,1};
static MTX BB = {0};
static MTX H = {0};

unsigned int AA = 0, B = 0; //, C = 0, A2 = 0;

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
  scanf("%d", &a);                           // fgets(line, LINESIZE, stdin); という手も
}

// OP型を正規化する
OP conv(OP f)
{
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


// 多項式の次数(degのOP型)
int odeg(OP f)
{
  int i, j = 0, k;

  if (f.t[0].a == 0)
    return 0;

  // k=terms(f);
  for (i = 0; i < DEG; i++)
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
  size_t len = (size_t)((deg(a) > deg(b) ? deg(a) : deg(b)) + 1);

  xor_u16(c.x, a.x, len);
  xor_u16(c.x, b.x, len);

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

// 20200816:正規化したいところだがうまく行かない
// 多項式の足し算
OP oadd(OP f, OP g)
{
  f = conv(f);
  g = conv(g);
  assert(op_verify(f));
  assert(op_verify(g));

  vec a = {0}, b = {0}, c = {0};
  int i, j, k, l = 0;
  OP h = {0}, f2 = {0}, g2 = {0};

  a = o2v(f);
  b = o2v(g);

  // k=deg(o2v(f));
  // l=deg(o2v(g));

  for (i = 0; i < DEG; i++)
  {
    c.x[i] = a.x[i] ^ b.x[i];
  }
  h = v2o(c);
  // h=conv(h);
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
  oterm o1 = {0}, o2 = {
                      0};

  n1 = terms(f);
  printf("n1=%d\n", n1);
  n2 = terms(g);
  printf("n2=%d\n", n2);
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
  printf("m1=%d\n", m1);
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
        printf("f.n==%d %d %d %d\n", f.t[i].n, g.t[j].n, i, j);
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
  int i, k, j;
  OP h = {0};
  vec test;
  unsigned int n;

  // f=conv(f);
  // k = deg (o2v(f));
  j = 0;
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
  int i, count = 0, k, l;
  oterm t = {0};
  OP h = {0}, e = {0}, r = {0};
  vec c = {0};

  k = odeg(f);
  l = odeg(g);
  if (l > k)
  {
    k = l;
  }

  for (i = 0; i < k + 1; i++)
  {
    t = g.t[i];
    e = oterml(f, t);
    h = oadd(h, e);
  }
  assert(op_verify(h));
  return h;
}

// リーディグタームを抽出(default)
oterm LT(OP f)
{
  int i, k;
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
  int i, k;
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
  int i, j, k;
  vec a, b;

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

  n = deg(a);

  // printf ("baka\n");
  assert(("baka\n", n >= 0));

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
  int i = 0, j, n, k;
  OP h = {0}, e = {0}, tt = {0};
  oterm a, b = {0}, c = {0};

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

  k = odeg(g);
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
  int i = 0, j, n, k;
  OP h = {0}, e = {
                  0};
  oterm a, b = {0}, c = {0};

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

  k = LT(g).n;
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
  int i = 1, j = 0, l;

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
unsigned short
trace(OP f, unsigned short x)
{
  int i, d;
  unsigned short u = 0;

  d = deg(o2v(f));

  for (i = 0; i < d + 1; i++)
  {
    u ^= gf[mlt(fg[f.t[i].a], mltn(f.t[i].n, fg[x]))];
  }

  return u;
}

// invert of polynomial
OP inv(OP a, OP n)
{
  OP d = {0}, x = {0}, s = {0}, q = {0}, r = {0}, t = {0}, u = {0}, v = {0}, w = {0}, tt = {0}, gcd = {0}, tmp = {0};
  oterm b = {0};
  vec vv = {0}, xx = {0}, aa = {0}, bb = {0}, cc = {0};

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

  b = LT(w);
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
  int i, j;

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
}

void printsage(vec a)
{
  int i, j, k;
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

// gcd for patterson
OP zgcd(OP a, OP n, int P)
{
  OP d = {0}, x = {0}, s = {0}, q = {0}, r = {0}, t = {0}, u = {0}, v = {0}, w = {0}, tt = {0}, gcd = {0}, rt = {0};
  oterm b = {0};
  vec vv = {0}, xx = {
                    0};

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

  tt = n;

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
  gcd = d; // $\gcd(a, n)$

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
  if (deg(o2v(xx)) < deg(o2v(yy)))
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

  for (i = 0; i < deg(v) + 1; i++)
  {
    e = v.x[i];
    d ^= (e << i);
  }

  return d;
}

// 配列からベクトル表現の多項式へ変換する
vec Setvec(int n)
{
  int i;
  vec v = {0};

  for (i = 0; i < n; i++)
  {
    v.x[n - 1 - i] = c[i];
  }

  return v;
}

OP ww[T] = {0};

OP bib(int i, OP d)
{
  int id, j;

  OP t[T] = {0};
  // omp_set_num_threads(omp_get_max_threads());
  id = omp_get_thread_num();
  t[id] = d;

  // #pragma omp parallel for
  for (j = 0; j < T; j++)
  {
    // #pragma omp critical
    if (i != j)
    {
      t[id] = omul(t[id], ww[j]);
    }
  }

  return t[id];
}

// chen探索
vec chen(OP f)
{
  vec e = {0};
  int i, count = 0, n, x = 0;
  unsigned short z;

  n = odeg((f));
  // exit(1);
  // #pragma omp parallel for private(i)
  for (x = 0; x < N; x++)
  {
    z = 0;
    // #pragma omp parallel for reduction (^:z)
    for (i = 0; i < n + 1; i++)
    {
      if (f.t[i].a > 0)
        z ^= gf[mlt(mltn(f.t[i].n, fg[x]), fg[f.t[i].a])];
    }
    if (z == 0)
    {
      e.x[count] = x;
      count++;
      printf("x=%d\n", x);
    }
  }
  printpol(e);
  printf(" ==eee!\n");
  // exit(1);

  return e;
}

// GF(2^m) then set m in this function.
int ben_or(OP f)
{
  int i, n;
  OP d = {0};
  n = deg(o2v(f));

  if (n == 0)
    return -1;

  OP t[K]={0};
#pragma omp parallel num_threads(omp_get_max_threads())
  {
#pragma omp for schedule(static)
    for (i = 0; i < n / 2 + 1; i++)
      t[i] = bib(i, d);
  }

  return 0;
}

// 配列の値を係数として多項式に設定する
OP setpol(unsigned short f[], int n)
{
  OP g;
  vec a;
  int i;

  memset(c, 0, sizeof(c));
  memcpy(c, f, 2 * n);
  a = Setvec(n);

  g = v2o(a);

  return g;
}

// バイナリ型パリティチェック行列を生成する
void bdet()
{
  int i, j, k, l;
  unsigned char dd[E * K] = {0};
  FILE *ff;

  // ff = fopen("Hb.key", "wb");

  for (i = 0; i < N; i++)
  {
    for (j = 0; j < K; j++)
    {
      l = mat[i][j];
      // #pragma omp parallel for
      for (k = 0; k < E; k++)
      {
        BB.x[i][j * E + k] = l % 2;
        l = (l >> 1);
      }
    }
  }

  for (i = 0; i < N; i++)
  {
    // #pragma omp parallel for
    for (j = 0; j < E * K; j++)
    {
      printf("%d,", BB.x[i][j]);
      // dd[j] = BH[j][i];
    }
    // fwrite(dd, 1, E * K, ff);
    printf("\n");
  }

  // fclose(ff);
}

unsigned short HH[N][K];

void toByte(const MTX *SH)
{
  vec v = {0};
  int i, j, k, cnt;

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
        HH[i][j] = x;//v2i(v);
      // printf("%d,", HH[i][j]);
      //= BH[j][i];
    }
    // fwrite(dd, 1, E * K, ff);
    // printf("\n");
  }
  printf("end of byte\n");
  // exit(1);
}

// 秘密置換を生成する
void Pgen()
{
  unsigned int i, j;
  FILE *fp;

  //fp = fopen("P.key", "wb");
  for (i = 0; i < N; i++)
    P[i] = i;
  random_shuffle(P, SIZE_OF_ARRAY(P));
  //  random_permutation(P);
  for (i = 0; i < N; i++)
    inv_P[P[i]] = i;
  //fwrite(P, 2, N, fp);
  //fclose(fp);

  // for (i = 0; i < N; i++)
  // printf ("%d,", inv_P[i]);
  // printf ("\n");

  // fp = fopen("inv_P.key", "wb");
  // fwrite(inv_P, 2, N, fp);
  // fclose(fp);
}

// 有限体の元の平方を計算する
int isqrt(unsigned short u)
{
  int i, j, k;

  for (i = 0; i < N; i++)
  {
    if (gf[mlt(i, i)] == u)
      return i;
  }

  printf("来ちゃいけないところに来ました\n");
  exit(1);
}

// 多項式の平方を計算する
OP osqrt(OP f, OP w)
{
  int i, j, k, jj, n, flg = 0;
  OP even = {0}, odd = {0}, h = {0}, r = {0}, ww = {0}, s = {0}, tmp = {0}, t = {0};
  oterm o = {0};
  vec v = {0};

  j = 0;
  jj = 0;
  k = deg(o2v(f));
  for (i = 0; i < k + 1; i++)
  {
    if (f.t[i].n % 2 == 0 && f.t[i].a > 0)
    {
      even.t[j].n = f.t[i].n / 2;
      even.t[j++].a = gf[isqrt(f.t[i].a)];
      printf("a=%d %d\n", f.t[i].a, i);
    }
    if (f.t[i].n % 2 == 1 && f.t[i].a > 0)
    {
      odd.t[jj].n = (f.t[i].n - 1) / 2;
      odd.t[jj++].a = gf[isqrt(f.t[i].a)];
      printf(" odd %d\n", i);
      flg = 1;
    }
  }

  k = deg(o2v(w));
  // printf ("%d\n", k);
  // exit(1);
  j = 0;
  jj = 0;
  for (i = 0; i < k + 1; i++)
  {
    if (w.t[i].n % 2 == 0 && w.t[i].a > 0)
    {
      h.t[j].a = gf[isqrt(w.t[i].a)];
      h.t[j++].n = w.t[i].n / 2;
      printf("h==%d %d\n", (w.t[i].n / 2), i);
    }
    if (w.t[i].n % 2 == 1 && w.t[i].a > 0)
    {
      r.t[jj].a = gf[isqrt(w.t[i].a)];
      r.t[jj++].n = (w.t[i].n - 1) / 2;
      printf("r=====%d %d\n", (w.t[i].n - 1) / 2, i);
    }
  }
  printpol(o2v(r));
  printf(" sqrt(g1)=======\n");

  //  exit(1);
  if (LT(r).n > 0)
  {
    s = inv(r, w);
  }
  else if (LT(r).n == 0)
  {
    printpol(o2v(r));
    printf(" deg(r)======0!\n");
    printpol(o2v(w));
    printf(" goppa======0\n");
    printpol(o2v(f));
    printf(" syn======0\n");
    if (LT(r).a > 0)
    {
      s.t[0].a = gf[isqrt(LT(r).a)];
      s.t[0].n = 0;
      return omod(oadd(even, omul(coeff(h, s.t[0].a), odd)), w);
    }
    return even;
    // s = inv (r, w);
    // wait ();
    // exit (1);
  }
  if (deg(o2v(s)) > 0)
    tmp = omod(omul(s, r), w);
  if (odeg((tmp)) > 0)
  {
    // printpol (o2v (tmp));
    printf(" r is not inv==========\n");
    wait();
    exit(1);
  }
  if (LT(h).n > 0 && odeg(s) > 0)
  {
    ww = omod(omul(h, s), w);
  }
  if (LT(h).n == 0 || odeg(s) == 0)
  {
    printpol(o2v(h));
    printf(" h=========0\n");
    exit(1);
  }

  if (LT(ww).n == 0 && LT(ww).a == 0)
  {
    printpol(o2v(s));
    printf(" s===========\n");
    printsage(o2v(w));
    printf(" w==============\n");
    printpol(o2v(r));
    printf(" r===========\n");
    printpol(o2v(h));
    printf(" h============\n");
    printpol(o2v(ww));
    printf(" ww==============\n");
    printf(" wwが0になりました。error\n");
    wait();
    // return ww;;
    exit(1);
  }

  tmp = omod(omul(ww, ww), w);
  if (LT(tmp).n == 1)
  {
    printpol(o2v(ww));
    printf(" ww succsess!===========\n");
  }
  else
  {
    // printpol (o2v (tmp));
    printf(" mod w^2==========\n");
    // printpol (o2v (ww));
    printf(" ww^2 failed!========\n");
    printpol(o2v(s));
    printf(" g1^-1==============\n");
    // printpol (o2v (w));
    printf(" w==============\n");
    // printpol (o2v (h));
    printf(" g0===========\n");
    // printpol (o2v (r));
    printf(" r===========\n");
    printf("この鍵では逆元が計算できません。error");
    wait();
    // return ww;
    exit(1);
  }

  //    exit(1);
  printpol(o2v(s));
  printf(" g1^-1=========\n");
  printpol(o2v(h));
  printf(" g0=========\n");
  // exit(1);
  printpol(o2v(ww));
  printf(" ww==========\n");
  //  exit(1);
  h = ww;
  if (odeg(omod(omul(h, ww), w)) == 1)
  {
    ww = h;
    h = omod(oadd(even, omul(ww, odd)), w);
    return h;
  }
  else if (LT(ww).a == 0)
  {
    printf("vaka\n");
    exit(1);
  }

  // //printpol(o2v(ww));
  printf(" 来ちゃだめなところに来ました\n");

  exit(1);
}

// パターソンアルゴリズムでバイナリGoppa符号を復号する
vec patterson(vec w_vec, vec f_vec)
{
  OP w = v2o(w_vec), f = v2o(f_vec);
  OP g1 = {0}, ll = {0}, op = {0};
  int i, j, k, l, c, n, count = 0, flg2 = 0;
  int flg, o1 = 0;
  OP tt = {0}, ff = {0}, h = {0};
  EX hh = {0}, opp = {0};
  vec v;
  oterm rr = {0};
  OP r2 = {0}, b2 = {
                   0};
  // unsigned short g[K+1]={2,2,12,1,2,8,4,13,5,10,8,2,15,10,7,3,5};

  tt.t[0].n = 1;
  tt.t[0].a = 1;
  o1 = 0;

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
  printf(" g1*g1%w===========\n");
  if (deg(o2v(b2)) == 0)
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
  flg = 0;
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
  v = chen(ll);
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
    flg = 1;
    printpol(o2v(ff));
    printf(" ==========beta!\n");
    printpol(o2v(ppo));
    printf(" alpha!=========\n");
    ll = oadd(omul(ff, ff), omul(tt, omul(ppo, ppo)));
    v = chen(ll);
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
  count = 0;
  printpol(o2v(ll));
  printf(" ll=================\n");
  printpol(o2v(f));
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

unsigned short vb[K * 2][N] = {0};
unsigned short gt[K * 2][K * 2] = {0};

void van()
{
  int i, j, k;

  printf("van der\n");

  for (i = 0; i < N; i++)
    vb[0][i] = 1;
  // #pragma omp parallel for private(i, j)
  for (i = 1; i < K; i++)
  {
    for (j = 0; j < N; j++)
      vb[i][j] = gf[mltn(i, fg[j])];
  }
}

void ogt(unsigned short g[K+1])
{
  int i, j, k;
  OP w = {0};
  unsigned short abc[N][K] = {0};

  // #pragma omp parallel for private(i, j)
  for (i = 0; i < K; i++)
  {
    for (j = 0; j < K - i; j++)
      gt[i][j + i] = g[j];
  }
}

unsigned short dd[N][N] = {0};

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

  OP pol = {0};
  pol = setpol(x.x, K + 1);
  printpol(o2v(pol));
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
  static MTX a = {0};
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

OP mkg()
{
  int i, j, k, l, ii = 0;
  OP w = {0}, f;
  unsigned short tr[N] = {0};
  unsigned short ta[N] = {0};

  unsigned short g1[3] = {1, 1, 3};

  // w=setpol(g,K+1);
  // printpol(o2v(w));
  // printf("\n");
  // exit(1);
  // w=omul(f,f);

  j = 0;
aa:

  // printf("\n");

  // 既約性判定のためのBen-Orアルゴリズム。拡大体にも対応している。デフォルトでGF(8192)
  // 既約多項式しか使わない。

  ii = 0;
  l = -1;
  vec pp = {0}, tt = {0};

  
          while (l < 0)
          {
              for (i = 0; i < K; i++)
                  pp.x[i] = rand() % N;
              mykey(tt.x, pp);
              tt.x[K] = 1;
              l = ben_or(v2o(tt));
              if (l == 0)
              {
                  printf("\n");
                  printsage(tt);
                  printf(" ==irr\n");
                 //exit(1);
              }
          }
          //w = v2o(tt);
  
  w = setpol(tt.x,K+1);

  // 多項式の値が0でないことを確認
  for (i = 0; i < N; i++)
  {
    ta[i] = trace(w, i);
    if (ta[i] == 0)
    {
      printf("trace 0 @ %d\n", i);
      // fail = 1;
      goto aa;
      //exit(1);
    }
  }
  for (i = 0; i < N; i++)
  {
    tr[i] = oinv(ta[i]);
    // printf("%d,", tr[i]);
  }

  printpol(o2v(w));
  printf(" =irreducible\n");
  printsage(o2v(w));
  printf("\n");
  // j++;
  // if(j>100)
  // exit(1);

  // 多項式を固定したい場合コメントアウトする。
  /*
  memset(ta, 0, sizeof(ta));
  w = setpol(g, K + 1);
  printpol(o2v(w));
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
  oprintpol(w);
  printf("\n");
  printsage(o2v(w));
  printf("\n");
  printf("sagemath で既約性を検査してください！\n");
  */

  van();
  //for(int i=0;i<K;i++)
  //g[i]=tt.x[i];
  ogt(tt.x);
  memset(mat, 0, sizeof(mat));

  // wait();

  // #pragma omp parallel for

  printf("\nすげ、オレもうイキそ・・・\n");
  // keygen(g);
  // exit(1);

  for (j = 0; j < N; j++)
  {
    for (i = 0; i < K; i++)
    {
      ma[j][i] = gf[mlt(fg[vb[i][j]], tr[j])];
    }
    // printf("tr[%d]=%d\n",j,tr[j]);
  }

  unsigned short s;
  // #pragma omp parallel for default(none) private(i, j, k, s) shared(mat, gt, ma, gf, fg)
  for (i = 0; i < K; i++)
  {
    for (j = 0; j < N; j++)
    {
      s = 0;

      for (k = 0; k < K; k++)
        s ^= gf[mlt(fg[gt[k][i]], fg[ma[j][k]])];
       printf("%d,",s);
      mat[j][i] = s;
    }
  }
  // printf("\n");
  // exit(1);

  for (j = 0; j < N; j++)
  {
    for (i = 0; i < K; i++)
      printf("%d,", mat[j][i]);
    printf("\n");
  }
  // exit(1);
  // wait();

  return w;
}

// Niederreiter暗号の公開鍵を作る
vec pubkeygen(void)
{
  int i, j, k, l;
  unsigned short n[K] = {0};
  FILE *fp;
  unsigned short dd[K] = {0};
  OP w = {0};
  vec v = {0};
  static MTX O, Q;

  w = mkg();
  printpol(o2v(w));
  printf(" ==goppa polynomial\n");

  v = o2v(w);
  //fp = fopen("sk.key", "wb");
  //fwrite(g, 2, K + 1, fp);
  //fclose(fp);
  oprintpol(w);
  printf("\n");
  printsage(o2v(w));
  printf("\n");
  printf("sagemath で既約性を検査してください！\n");

  bdet();
  //  toByte(BB);
  // exit(1);

  Pgen();
  fp = fopen("P.key", "w");
  //  fwrite(P, 2, N, fp);
  fclose(fp);
  // makeS();
  do
  {
    memset(Q.x, 0, sizeof(Q.x));
    memset(O.x, 0, sizeof(O.x));
    memset(S.x, 0, sizeof(S.x));
  
    for (i = 0; i < (K)*E; i++)
    {
      for (j = 0; j < (K)*E; j++)
        S.x[i][j] = xor128() % 2;
    }
  } while (is_reg(&S, &inv_S) == -1);

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
  H = mulmat(&S, &BB, 1);
  for (i = 0; i < K * E; i++)
  {
    for (j = 0; j < N; j++)
      O.x[j][i] = H.x[P[j]][i];
  }
  toByte(&O);

  return v;
}

OP dec(unsigned short ss[])
{
  int i, j, k;
  vec v = {0};
  OP s = {0};
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
      h2o[i] ^= (ch[j] & inv_S.x[i][j]);
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

int ero2(vec v)
{
  int i, j, count = 0;
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
    int cnt = 0;
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
            AA++;
            break;
            //exit(1);
          }
          cnt++;
        }
      }
*/
  }
  for (i = 0; i < N; i++)
    ya[i] = xa[P[i]];
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
    B++;
  }
  if (count < T * 2)
  {
    printf("error is too few\n");

    AA++;
    // memcpy (zz, z1, sizeof (zz));
    /*
      printf("{");
      for (i = 0; i < D; i++)
        printf("%d,", z1[i]);
      printf("};\n");
      printpol(o2v(w));
      printf(" =========goppa\n");
      printsage(o2v(w));
      printf(" =========sage\n");
      printsage(o2v(f));
      printf(" =========syn\n");
      printpol(o2v(f));
      printf(" ==========synd\n");
      */
    printf("へげえええーっ\n");
    // exit(1);
    return -1;
    //exit(1);
  }

  return count;
}

void mkerr(unsigned short *z1, int num)
{
  int i, j, l;
  j = 0;

  while (j < num)
  {
    l = xor128() % N;
    // printf ("l=%d\n", l);
    if (0 == z1[l])
    {
      z1[l] = 1;
      printf("l=%d %d\n", l, j);
      j++;
    }
  }

  for (i = 0; i < N; i++)
  {
    if (z1[i] > 0)
      printf("%d=%d\n", i, z1[i]);
  }
}

vec zin(unsigned short zz[])
{
  int i, j;
  OP s = {0};
  vec v = {0};
  unsigned short ss[K] = {0};

  for (i = 0; i < N; i++)
  {
    if (zz[i] > 0)
    {
      for (j = 0; j < K; j++)
      {
        ss[j] ^= HH[i][j];
        // printf("%d,", HH[i][j]);
      }
      // printf("\n");
    }
  }

  // 暗号化されたシンドロームを復元する
  s = dec(ss);
  for (j = 0; j < K; j++)
    printf("%d,", ss[j]);
  printf(" ==ss\n");
  // exit(1);

  return o2v(s);
}
