// Niederreiter cryptosystem with Patterson decoding.
//
// 整理・高速化版。外部から見える関数の名前・引数は変えていない。
//
//  - デバッグ出力は -DVERBOSE のときだけ出す(通常ビルドでは出ない)
//  - OP型の剰余・商・逆元・gcd は内部を vec で計算する(変換は入口と出口で1回だけ)
//  - equ / isqrt / v2a の O(N) 総当たりを O(1) の式に置き換え
//  - Chien探索は Horner 法
//  - patterson(): 同一入力の zgcd / chen の重複呼び出しを1回にし、
//                 sqrt(x) mod w を鍵ごとにキャッシュ
//  - dec(): ビット分解を直接のビット演算に変更

#include "common.h"
#include "systematic.h"
#include "sinv_pack.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <stdint.h>

#include <time.h>

// ---------------------------------------------------------------
// デバッグ出力 (make CPPFLAGS=-DVERBOSE で有効)
// ---------------------------------------------------------------
#ifdef VERBOSE
#define DBG(...) printf(__VA_ARGS__)
#else
#define DBG(...) ((void)0)
#endif

static inline void dbg_pol(const char *label, vec a)
{
#ifdef VERBOSE
  printpol(a);
  printf(" %s\n", label);
#else
  (void)label;
  (void)a;
#endif
}

static double now(void)
{
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec + t.tv_nsec * 1e-9;
}

// ===============================================================
// 有限体の基本演算
// ===============================================================

// 有限体の元の逆数 (戻り値は mlt にそのまま渡せる対数表現)
unsigned short oinv(unsigned short a)
{
  if (a == 0)
    return -1;

  return (N - fg[a]) % (N - 1) + 1;
}

// 有限体の元の逆数 (総当たり版。比較用に残してある)
unsigned short oinv2(unsigned short a)
{
  int i;

  if (a == 0)
    return -1;

  for (i = 0; i < N; i++)
  {
    if (gf[mlt(fg[a], i)] == 1)
      return (unsigned short)i;
  }

  fprintf(stderr, "no return\n");
  exit(1);
}

// a に何をかけたら b になるか (= b / a)
unsigned short equ(unsigned short a, unsigned short b)
{
  if (b == 0)
    return 0;
  if (a == 0)
    return N; // 元の総当たりが「見つからない」ときに返していた値
  return gf[mlt(fg[b], oinv(a))];
}

// 平方根 (戻り値は mlt にそのまま渡せる対数表現。gf[isqrt(u)] が sqrt(u))
// 2*(i-1) = log(u) (mod M-1) を解く。M-1 は奇数なので必ず解ける。
int isqrt(unsigned short u)
{
  if (u == 0)
    return 0;

  const int l = fg[u] - 1;
  return ((l + ((l & 1) ? (M - 1) : 0)) >> 1) + 1;
}

// ===============================================================
// OP型 <-> vec型
// ===============================================================

// OP型からベクトル型への変換
vec o2v(OP f)
{
  vec a = {0};
  int i;

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
  int a;
  printf(" (enter number and hit return) ");
  fflush(stdout);
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

  return v2o(o2v(f));
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

// 零多項式を 0 次として扱う次数 (odeg と同じ規約)
static inline int degz(const vec *a)
{
  const int d = deg(a);
  return d < 0 ? 0 : d;
}

// 多項式の次数(default)
int vdeg(vec a)
{
  int i, n = 0, flg = 0;

  for (i = DEG - 1; i >= 0; i--)
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

  for (i = DEG - 1; i >= 0; i--)
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
      h.t[f.t[i].n].a = f.t[i].a;
  }

  return h;
}

// ===============================================================
// vec型の多項式演算 (係数は GF(2^m) の元そのもの)
// ===============================================================

vec vadd(vec a, vec b)
{
  vec c = {0};
  size_t len = (size_t)((deg(&a) > deg(&b) ? deg(&a) : deg(&b)) + 1);

  xor_u16(c.x, a.x, len);
  xor_u16(c.x, b.x, len);

  return c;
}

// 注: 内側・外側とも「次数未満」(<) で回しているため最高次の項が抜ける。
//     テストとの互換のため挙動は変えていない (fmul が正しい掛け算)。
vec vmul(vec a, vec b)
{
  vec c = {0};
  int k = deg(&a);
  int l = deg(&b);

  for (int i = 0; i < k; i++)
  {
    if (a.x[i] == 0)
      continue;
    const int la = fg[a.x[i]];
    for (int j = 0; j < l; j++)
      c.x[i + j] ^= gf[mlt(la, fg[b.x[j]])];
  }
  return c;
}

vec vmula(vec a, vec b)
{
  return vmul(a, b);
}

// c = a + b (係数ごとの XOR。次数までしか回さない)
static inline vec gf_add(const vec *a, const vec *b)
{
  vec c = {0};
  const int da = deg(a), db = deg(b);
  const int n = da > db ? da : db;

  for (int i = 0; i <= n; i++)
    c.x[i] = a->x[i] ^ b->x[i];
  return c;
}

// c = a * b (DEG 次以上は切り捨て)
// 対数表現を内側ループの外で1回だけ引き、次数の範囲だけ回す
static vec gf_mul(const vec *a, const vec *b)
{
  vec c = {0};
  const int da = deg(a), db = deg(b);
  unsigned short lb[DEG];

  if (da < 0 || db < 0)
    return c;

  for (int j = 0; j <= db; j++)
    lb[j] = fg[b->x[j]];

  for (int i = 0; i <= da; i++)
  {
    if (a->x[i] == 0)
      continue;
    const int la = fg[a->x[i]];
    const int jmax = (db < DEG - 1 - i) ? db : DEG - 1 - i;
    for (int j = 0; j <= jmax; j++)
    {
      if (lb[j])
        c.x[i + j] ^= gf[mlt(la, lb[j])];
    }
  }
  return c;
}

vec fmul(vec a, vec b)
{
  return gf_mul(&a, &b);
}

// f = q*g + r を一度に求める。g は非零。q は NULL でもよい。
static void gf_divmod(const vec *f, const vec *g, vec *q, vec *r)
{
  const int dg = deg(g);

  if (dg < 0)
  {
    fprintf(stderr, "vmod: division by zero polynomial\n");
    exit(1);
  }

  *r = *f;
  if (q)
    memset(q->x, 0, sizeof q->x);

  int dr = deg(r);
  if (dr < dg)
    return;

  unsigned short lg[DEG];
  for (int j = 0; j <= dg; j++)
    lg[j] = fg[g->x[j]];

  const int inv_lead = oinv(g->x[dg]);

  while (dr >= dg)
  {
    const int qc = gf[mlt(fg[r->x[dr]], inv_lead)];
    const int lq = fg[qc];
    const int sh = dr - dg;

    if (q)
      q->x[sh] = qc;

    // r ^= qc * x^sh * g  (最高次の項はここで消える)
    for (int j = 0; j <= dg; j++)
    {
      if (lg[j])
        r->x[j + sh] ^= gf[mlt(lg[j], lq)];
    }

    // 次の最高次へ
    do
    {
      dr--;
    } while (dr >= 0 && r->x[dr] == 0);
  }
}

// vec 多項式の剰余 (f mod g)
vec vmod(vec f, vec g)
{
  vec r;
  gf_divmod(&f, &g, NULL, &r);
  return r;
}

// omod と同じ規約の剰余: g が定数/零、または deg f < deg g なら f をそのまま返す
static vec gf_omod(const vec *f, const vec *g)
{
  vec r;

  if (deg(g) <= 0 || deg(f) < deg(g))
    return *f;

  gf_divmod(f, g, NULL, &r);
  return r;
}

// ===============================================================
// OP型の多項式演算 (内部は vec で計算する)
// ===============================================================

// 多項式の足し算
OP oadd(OP f, OP g)
{
  vec a = o2v(f), b = o2v(g);
  vec c = gf_add(&a, &b);
  return v2o(c);
}

// 項の順序を昇順に揃える
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
    if (f.t[i].a > 0)
    {
      for (j = i; j < k + 1; j++)
      {
        if (s.n < f.t[j].n)
        {
          s.n = f.t[j].n;
          s.a = f.t[j].a;
        }
      }
    }
  }

  return s;
}

// 多項式を足し算する（OP型・旧実装。同じ次数の項は係数に関係なく消える点に注意）
OP add(OP f, OP g)
{
  unsigned long long int i, j, n1 = terms(f), n2 = terms(g), count = 0;
  OP h = {0};

  for (i = 0; i < n1 + 1; i++)
  {
    for (j = 0; j < n2 + 1; j++)
    {
      if (f.t[i].n == g.t[j].n && g.t[j].a > 0 && f.t[i].a > 0)
      {
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

  return h;
}

// 多項式を項ずつ掛ける
OP oterml(OP f, oterm t)
{
  f = conv(f);
  assert(op_verify(f));
  int i;
  OP h = {0};

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
  vec a = o2v(f), b = o2v(g);
  vec c = gf_mul(&a, &b);
  return v2o(c);
}

// リーディグタームを抽出(default)
oterm LT(OP f)
{
  oterm t = {0};

  for (int i = DEG - 1; i >= 0; i--)
  {
    if (f.t[i].a > 0)
    {
      t.n = f.t[i].n;
      t.a = f.t[i].a;
      break;
    }
  }

  return t;
}

// 多項式の最初の項を抽出
oterm LT2(OP f)
{
  oterm t = {0};

  t.n = f.t[0].n;
  t.a = f.t[0].a;

  return t;
}

// 多項式を単項式で割る
oterm LTdiv(OP f, oterm t)
{
  oterm tt = {0}, s = {0};

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
  else
  {
    s.n = tt.n - t.n;
    s.a = equ(t.a, tt.a);
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
  assert(n >= 0);

  for (i = n; i > -1; i--)
  {
    if (a.x[i] > 0)
    {
      printf("%u", a.x[i]);
      printf("x^%d", i);
      printf("+");
    }
  }
}

// 多項式の商を取る
// (元の挙動: g が零なら異常終了。g が定数なら f/g。deg f < deg g のときは f を返す)
OP odiv(OP f, OP g)
{
  vec vf = o2v(f), vg = o2v(g), q, r;
  const int dg = deg(&vg);

  if (dg < 0)
  {
    print_trace();
    exit(1);
  }

  if (dg == 0)
  {
    const unsigned short c = vg.x[0];
    if (c > 1)
    {
      const int ic = oinv(c);
      const int df = deg(&vf);
      for (int i = 0; i <= df; i++)
        vf.x[i] = gf[mlt(fg[vf.x[i]], ic)];
    }
    return v2o(vf);
  }

  if (degz(&vf) < dg)
    return v2o(vf);

  gf_divmod(&vf, &vg, &q, &r);
  return v2o(q);
}

// 多項式のべき乗
OP opow(OP f, int n)
{
  int i;
  OP g = f;

  for (i = 1; i < n; i++)
    g = omul(g, f);

  return g;
}

// 多項式の剰余を取る
OP omod(OP f, OP g)
{
  vec vf = o2v(f), vg = o2v(g), r;

  if (deg(&vg) <= 0 || deg(&vf) < deg(&vg))
    return f;

  gf_divmod(&vf, &vg, NULL, &r);
  return v2o(r);
}

// 多項式のべき乗余
OP opowmod(OP f, OP mod, int n)
{
  int i = 1;

  while (i < n + 1)
  {
    f = omod(omul(f, f), mod);
    i++;
  }

  return f;
}

// ===============================================================
// 逆元・gcd (EEA)
// ===============================================================

// 多項式の逆元 a^{-1} mod n (vec版)
static vec gf_inv(vec a, vec n)
{
  vec d, x = {0}, s = {0}, t, q, r = {0}, v, u;

  s.x[0] = 1;

  if (degz(&a) > degz(&n))
  {
    vec tmp = a;
    a = n;
    n = tmp;
  }
  if (deg(&a) < 0)
  {
    fprintf(stderr, " a ga 0\n");
    exit(1);
  }

  d = n;
  while (deg(&a) > 0)
  {
    gf_divmod(&d, &a, &q, &r);
    d = a;
    a = r;
    vec qs = gf_mul(&q, &s);
    t = gf_add(&x, &qs);
    x = s;
    s = t;
  }

  // ここで a は定数(= gcd)。x が係数
  d = a;
  x = s;

  v = gf_add(&x, &n);
  if (deg(&v) <= 0 || deg(&n) <= 0)
  {
    fprintf(stderr, "inv: invalid modulus\n");
    exit(1);
  }
  u = gf_omod(&v, &n);

  if (deg(&u) >= 0 && deg(&d) >= 0)
  {
    const unsigned short c = d.x[0];
    if (c > 1)
    {
      const int ic = oinv(c);
      const int du = deg(&u);
      for (int i = 0; i <= du; i++)
        u.x[i] = gf[mlt(fg[u.x[i]], ic)];
    }
  }

  if (deg(&u) < 0)
  {
    fprintf(stderr, "no return at u==0\n");
    exit(1);
  }

  return u;
}

// invert of polynomial
OP inv(OP a, OP n)
{
  return v2o(gf_inv(o2v(a), o2v(n)));
}

unsigned short v2a(oterm a)
{
  if (a.a == 0)
    return 0;

  return fg[a.a] - 1;
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
      printf("B('a^%d')*X**%d+", j, i); // for GF(2^m)
    }
  }
}

// vec 多項式の剰余 (旧実装。比較用に残してある)
vec vmod2(vec f, vec g)
{
  vec r = f;
  int dg = deg(&g);

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

    int q = gf[mlt(fg[r.x[dr]], oinv(g.x[dg]))];

    for (int j = 0; j <= dg; j++)
    {
      if (g.x[j] == 0)
        continue;

      int k = j + shift;
      int product = gf[mlt(fg[g.x[j]], fg[q])];

      r.x[k] = fg[r.x[k]] ^ product;
      r.x[k] = (r.x[k] == 0) ? 0 : fg[r.x[k]];
    }

    if (deg(&r) >= dr)
    {
      printf("vmod: degree did not decrease\n");
      exit(1);
    }
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
    if (LT(v2o(yy)).a == 0)
      break;
  }

  return xx;
}

// gcd for patterson (部分EEA: 余りの次数が P 以下になるまで)
OP zgcd(OP a_, OP n_, int P)
{
  vec a = o2v(a_), n = o2v(n_), x = {0}, s = {0}, t, q, r;

  s.x[0] = 1;

  if (degz(&a) > degz(&n))
  {
    vec tmp = a;
    a = n;
    n = tmp;
  }
  if (deg(&a) < 0)
  {
    fprintf(stderr, " a ga 0\n");
    exit(1);
  }

  vec d = n;
  while (deg(&a) > P)
  {
    gf_divmod(&d, &a, &q, &r);
    d = a;
    a = r;
    vec qs = gf_mul(&q, &s);
    t = gf_add(&x, &qs);
    x = s;
    s = t;
  }
  x = s;

  dbg_pol("=======x", x);

  return v2o(x);
}

// gcd
OP gcd(OP xx_, OP yy_)
{
  vec xx = o2v(xx_), yy = o2v(yy_), tt;
  OP h = {0};

  h.t[0].a = 1;
  h.t[0].n = 0;

  if (deg(&xx) < deg(&yy))
  {
    vec tmp = xx;
    xx = yy;
    yy = tmp;
  }

  tt = gf_omod(&xx, &yy);
  while (deg(&tt) > 0)
  {
    xx = yy;
    yy = tt;
    if (deg(&yy) > 0)
      tt = gf_omod(&xx, &yy);
    if (deg(&tt) < 0)
      return v2o(yy);
  }
  if (deg(&tt) < 0)
    return v2o(yy);

  return h;
}

// ===============================================================
// 変換ユーティリティ
// ===============================================================

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
unsigned short v2i(vec v)
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
    v.x[n - 1 - i] = coefficients[i];

  return v;
}

// 配列の値を係数として多項式に設定する
OP setpol(unsigned short f[], int n)
{
  return v2o(Setvec(f, n));
}

static OP bib(int i, OP d, const OP ww[T])
{
  int j;
  OP t = d;

  for (j = 0; j < T; j++)
  {
    if (i != j)
      t = omul(t, ww[j]);
  }

  return t;
}

// Chien探索 (Horner法。根を昇順に e.x[] に並べる)
static vec chen(vec f)
{
  vec e = {0};
  int count = 0;
  const int n = deg(&f);

  for (int x = 0; x < N; x++)
  {
    unsigned short z = 0;

    if (n >= 0)
    {
      const int lx = fg[x];
      z = f.x[n];
      for (int i = n - 1; i >= 0; i--)
        z = gf[mlt(fg[z], lx)] ^ f.x[i];
    }

    if (z == 0 && count < DEG)
    {
      e.x[count++] = x;
      DBG("x=%d\n", x);
    }
  }
  dbg_pol("==eee!", e);

  return e;
}

vec vgcd(vec xx, vec yy)
{
  vec tt = {0};

  while (1)
  {
    int zero = 1;
    for (int i = 0; i < DEG; i++)
    {
      if (yy.x[i] != 0)
      {
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

vec derivative(vec f)
{
  vec df = {0};
  int y = deg(&f);

  if (y > 0)
  {
    for (int i = 1; i <= y; i += 2)
      df.x[i - 1] = f.x[i];
  }

  return df;
}

vec diff(vec a)
{
  vec b = derivative(a);

  return vgcd(a, b);
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
  vec tmp = o2v(f);
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

// 秘密置換を生成する
void Pgen(unsigned short permutation[N])
{
  unsigned int i;

  for (i = 0; i < N; i++)
    permutation[i] = i;
  random_shuffle(permutation, N);
}

// ===============================================================
// 平方根 (w = g0^2 + x*g1^2 のとき sqrt(x) = g0/g1)
// ===============================================================

// sqrt(x) mod w を求める。成功で0。鍵ごとに1回だけ呼べばよい。
int osqrt_x(OP w, OP *sx)
{
  OP h = {0}, r = {0}, s = {0}, ww = {0}, tmp = {0}, xx = {0};
  int i, j = 0, jj = 0;
  vec vtmp = o2v(w);
  int k = deg(&vtmp);

  for (i = 0; i < k + 1; i++)
  {
    if (w.t[i].a == 0)
      continue;
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
  if (jj == 0)
    return 1; // 奇数次の項がない: g が平方式で square-free でない

  if (LT(r).n > 0)
  {
    s = inv(r, w);
    tmp = omod(omul(s, r), w);
    if (odeg(tmp) != 0 || LT(tmp).a != 1)
      return 2; // 逆元が取れない
    ww = omod(omul(h, s), w);
  }
  else
  {
    // g1 が定数の場合
    s.t[0].a = gf[isqrt(LT(r).a)];
    s.t[0].n = 0;
    ww = coeff(h, s.t[0].a);
  }

  // 検証: ww^2 mod w が x と完全に一致するか
  xx.t[0].n = 1;
  xx.t[0].a = 1;
  tmp = oadd(omod(omul(ww, ww), w), xx);
  if (LT(tmp).a != 0)
    return 3;

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
    if (f.t[i].a == 0)
      continue;
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

// 多項式の平方根を計算する (互換用。呼ぶたびに osqrt_x を計算する)
OP osqrt(OP f, OP w)
{
  OP sx = {0};
  if (osqrt_x(w, &sx) != 0)
    exit(1);
  return osqrt_f(f, w, sx);
}

// 同じ w が続けて来る間は osqrt_x の結果を使い回す
static int sqrt_x_cached(const vec *w_vec, OP w, OP *sx)
{
  static _Thread_local vec cached_w;
  static _Thread_local OP cached_sx;
  static _Thread_local bool valid = false;

  if (valid && memcmp(cached_w.x, w_vec->x, sizeof cached_w.x) == 0)
  {
    *sx = cached_sx;
    return 0;
  }

  OP tmp = {0};
  const int rc = osqrt_x(w, &tmp);
  if (rc != 0)
  {
    valid = false;
    return rc;
  }

  cached_w = *w_vec;
  cached_sx = tmp;
  valid = true;
  *sx = tmp;
  return 0;
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

    if (ei)
      t++;
    if (ei && di)
      hit++; // 正しく見つけた
    if (!ei && di)
      ex++; // 誤って余分に出した
  }

  if (T_out)
    *T_out = t;
  if (extra)
    *extra = ex;
  return hit;
}

// ===============================================================
// 復号
// ===============================================================

vec paloma(OP s, OP f)
{
  vec x = {0}, v = {0}, zero = {0};
  x.x[1] = 1;
  v.x[0] = 1;
  OP one = v2o(v);

  if (odeg(s) < 0)
    return zero; // 誤りなし

  OP ss = omod(oadd(one, omul(s, v2o(x))), f);

  OP g1 = gcd(s, f);
  OP g2 = gcd(ss, f);
  OP g12 = odiv(f, omul(g1, g2));

  OP s1 = odiv(s, g1); // 常に割る
  OP s2 = odiv(ss, g2);

  OP denom = omul(g2, s1);
  if (odeg(gcd(denom, g12)) > 0)
  { // f が square-free なら起きないはず
    fprintf(stderr, "denominator not invertible mod g12\n");
    exit(1);
  }

  OP u = omod(omul(omul(g1, s2), inv(denom, g12)), g12);
  OP sp = osqrt(u, g12);

  /* 拡張ユークリッド: rb ≡ tb*sp (mod g12) を保つ */
  OP ra = g12, rb = sp;
  OP ta = {0}, tb = one;
  int Ba = K / 2 - odeg(g2); // a' の次数の上限

  while (odeg(rb) > Ba)
  { // 判定が先
    OP q = odiv(ra, rb);
    OP r = omod(ra, rb);
    OP t = oadd(ta, omul(q, tb));
    ra = rb;
    rb = r;
    ta = tb;
    tb = t;
  }
  OP a0 = rb, b0 = tb;

  OP a = omul(a0, g2);
  OP b = omul(b0, g1);
  OP sigma = oadd(omul(a, a), omul(omul(b, b), v2o(x)));

  int yt = odeg(sigma); // 誤りの個数 (K に固定しない)
  if (yt < 1 || yt > K)
  {
    fprintf(stderr, "decode failed\n");
    exit(1);
  }

  return chen(o2v(sigma));
}

// パターソンアルゴリズムでバイナリGoppa符号を復号する
vec patterson(vec w_vec, vec f_vec)
{
  OP w = v2o(w_vec), f = v2o(f_vec);
  OP tt = {0};

  tt.t[0].n = 1;
  tt.t[0].a = 1;

  // f の逆元 ff を求めて、本当に逆元になっているか確認する
  OP ff = inv(f, w);
  OP b2 = omod(omul(ff, f), w);
  if (odeg(b2) > 0)
  {
    fprintf(stderr, "逆元が計算できません。error\n");
    exit(1);
  }

  // g1 = sqrt(ff + x) mod w
  OP r2 = oadd(ff, tt);
  OP sx = {0};
  if (sqrt_x_cached(&w_vec, w, &sx) != 0)
  {
    fprintf(stderr, "sqrt(x) mod w が計算できません。error\n");
    exit(1);
  }
  OP g1 = osqrt_f(r2, w, sx);
  dbg_pol("sqrt(h+x)=======", o2v(g1));

  b2 = omod(omul(g1, g1), w);
  vec tmp = o2v(b2);
  if (deg(&tmp) == 0)
  {
    fprintf(stderr, "b2 is 0\n");
    exit(1);
  }
  if (LT2(b2).a != LT2(r2).a)
  {
    fprintf(stderr, "g1は平方ではありません。error\n");
    exit(1);
  }
  if (LT(g1).n == 0 && LT(g1).a == 0)
  {
    fprintf(stderr, "平方が０になりました。 error\n");
    exit(1);
  }

  // s と w が互いに素であること
  OP h = gcd(f, w);
  if (odeg(h) > 0)
  {
    fprintf(stderr, "s,wは互いに素じゃありません。\n");
    exit(1);
  }

  // 位置多項式 l = beta^2 + x*alpha^2 (zgcd は1回だけ)
  OP ppo = zgcd(w, g1, T);
  ff = omod(omul(ppo, g1), w);
  OP ll = oadd(omul(ff, ff), omul(tt, omul(ppo, ppo)));
  dbg_pol("ll=================", o2v(ll));

  vec v = chen(o2v(ll));
  if (v.x[K - 1] > 0)
    return v;

  if (odeg(ff) == 1)
  {
    fprintf(stderr, "deg==1\n");
    exit(1);
  }

  return v;
}

// ===============================================================
// 鍵生成
// ===============================================================

/* input: in0, in1 in GF((2^m)^t)*/
/* output: out = in0*in1 */
void GF_mul(unsigned short *out, unsigned short *in0, unsigned short *in1)
{
  // 既約多項式 x^K + (taps の項) (sage で求めたもの)
  static const unsigned char taps_512[] = {8, 5, 2, 0};
  static const unsigned char taps_256[] = {10, 5, 2, 0};
  static const unsigned char taps_128[] = {7, 2, 1, 0};
  static const unsigned char taps_96[] = {63, 61, 60, 59, 54, 50, 46, 44, 40, 37, 33,
                                          32, 28, 24, 23, 21, 20, 19, 18, 16, 15, 14,
                                          13, 11, 9, 8, 7, 6, 4, 3, 2, 0};
  static const unsigned char taps_64[] = {33, 30, 26, 25, 24, 23, 22, 21, 20, 18,
                                          13, 12, 11, 10, 7, 5, 4, 2, 1, 0};
  static const unsigned char taps_32[] = {15, 9, 7, 4, 3, 0};
  static const unsigned char taps_16[] = {5, 3, 2, 0};
  static const unsigned char taps_8[] = {4, 3, 2, 0};
  static const unsigned char taps_4[] = {1, 0};

  const unsigned char *taps = NULL;
  int ntaps = 0;
#define SET_TAPS(kk, arr)        \
  if (K == (kk))                 \
  {                              \
    taps = (arr);                \
    ntaps = (int)sizeof(arr);    \
  }
  SET_TAPS(512, taps_512)
  SET_TAPS(256, taps_256)
  SET_TAPS(128, taps_128)
  SET_TAPS(96, taps_96)
  SET_TAPS(64, taps_64)
  SET_TAPS(32, taps_32)
  SET_TAPS(16, taps_16)
  SET_TAPS(8, taps_8)
  SET_TAPS(4, taps_4)
#undef SET_TAPS

  unsigned short prod[K * 2 - 1] = {0};
  unsigned short l1[K];
  int i, j;

  for (j = 0; j < K; j++)
    l1[j] = fg[in1[j]];

  for (i = 0; i < K; i++)
  {
    if (in0[i] == 0)
      continue;
    const int la = fg[in0[i]];
    for (j = 0; j < K; j++)
    {
      if (l1[j])
        prod[i + j] ^= gf[mlt(la, l1[j])];
    }
  }

  for (i = (K - 1) * 2; i >= K; i--)
  {
    if (prod[i] == 0)
      continue;
    for (int t = 0; t < ntaps; t++)
      prod[i - K + taps[t]] ^= prod[i];
  }

  for (i = 0; i < K; i++)
    out[i] = prod[i];
}

#ifdef VERBOSE
// 多項式の代入値 (renritu の確認表示用)
static unsigned short trace(vec f, unsigned short x)
{
  unsigned short u = 0;
  int degree = deg(&f);

  for (int i = 0; i <= degree; i++)
    if (f.x[i] != 0)
      u ^= gf[mlt(fg[f.x[i]], mltn(i, fg[x]))];

  return u;
}
#endif

// 連立方程式を解く (ガウス・ジョルダン法)
vec renritu(const MTX *input)
{
  MTX a = {0};
  vec v = {0};

  if (input == NULL || input->rows < (size_t)K ||
      input->cols < (size_t)(K + 1) || matrix_copy(&a, input) != 0)
  {
    fprintf(stderr, "renritu: invalid or unavailable matrix storage\n");
    return v;
  }

  for (int i = 0; i < K; i++)
  {
    const int ip = oinv(a.x[i][i]);

    for (int j = 0; j < K + 1; j++)
      a.x[i][j] = gf[mlt(fg[a.x[i][j]], ip)];

    for (int j = 0; j < K; j++)
    {
      if (i == j)
        continue;
      const int ld = fg[a.x[j][i]];
      if (ld == 0)
        continue;
      for (int k = i; k < K + 1; k++)
        a.x[j][k] ^= gf[mlt(ld, fg[a.x[i][k]])];
    }
  }

  for (int i = 0; i < K; i++)
    v.x[i] = a.x[i][K];

#ifdef VERBOSE
  {
    vec x = {0};
    for (int i = 0; i < K; i++)
    {
      printf(" x%d = %d\n", i, v.x[i]);
      x.x[i + 1] = v.x[i];
    }
    x.x[0] = 1;

    vec pol = Setvec(x.x, K + 1);
    printpol(pol);
    printf(" ==key\n");
    for (int i = 0; i < N; i++)
    {
      if (trace(pol, i) == 0)
        printf("%d i=%d\n", i, fg[i] - 1);
    }
  }
#endif

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

  mat[0][0] = 1;
  for (i = 0; i < K; i++)
    mat[1][i] = x.x[i];

  for (j = 2; j <= K; j++)
    GF_mul(mat[j], mat[j - 1], x.x);

  for (i = 0; i < K; i++)
  {
    for (j = 0; j < K + 1; j++)
      a.x[i][j] = mat[j][i];
  }

  vec v = renritu(&a);

  for (i = 0; i < K; i++)
    out[i] = v.x[i];

  free(mat);
  matrix_free(&a);
  return 0;
}

// matrix[j][i] = (Σ_{m=0}^{i} g[m]·α_j^(i-m)) / g(α_j),  α_j = j
static void matrix_horner(const vec *polynomial,
                          const unsigned short *tr, // tr[j] = oinv(g(j))
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
  if (!tr || !ta)
  {
    free(tr);
    free(ta);
    return -1;
  }

  bool found = false;
  for (int attempt = 0; attempt < 100 && !found; attempt++)
  {
    memset(pp.x, 0, sizeof(pp.x));
    for (int i = 0; i < K; i++)
      pp.x[i] = rand() % N;
    pp.x[K] = 1;
    w = setpol(pp.x, K + 1);
    polynomial = o2v(w);

    int has_root = 0; /* 試行ごとに0に戻す */
#pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++)
    {
      int skip;
#pragma omp atomic read
      skip = has_root;
      if (skip)
        continue;
      ta[i] = i ? eval(polynomial.x, fg[i]) : polynomial.x[0];
      if (ta[i] == 0)
      {
#pragma omp atomic write
        has_root = 1;
      }
    }
    found = !has_root;
  }
  if (!found)
  {
    free(tr);
    free(ta);
    return -1;
  }

  for (int i = 0; i < N; i++)
    tr[i] = oinv(ta[i]);
  matrix_horner(&pp, tr, matrix);
  *result = o2v(w);
  free(tr);
  free(ta);
  return 0;
}



/*
 * BM用シンドローム sp[0..2t-1] から パターソン用 S(x) の係数 S[0..t-1] を作る。
 *
 *   sp[j] = sum_i e_i * a_i^j / g(a_i)^2    (j = 0..2t-1)
 *   S(x)  = sum_i e_i / (x - a_i) mod g(x)
 *
 * g[0..t] は g(x) の係数 (g[t] = 1 のモニック)。
 * 標数2なので加算・減算はすべて XOR。
 */
void bm_to_patterson_syndrome(unsigned short *S, const unsigned short *sp, const unsigned short *g, int t)
{
    unsigned short s[t];                          /* g型シンドローム s_j = sum e_i a_i^j / g(a_i) */

    /* 段階1: s_j = sum_{k=0..t} g_k * sp[k+j]   (j = 0..t-1, 最大添字 2t-1) */
    for (int j = 0; j < t; j++) {
        unsigned short acc = 0;
        for (int k = 0; k <= t; k++)
            acc ^= gf[mlt(fg[g[k]], fg[sp[k + j]])];
        s[j] = acc;
    }

    /* 段階2: S[j] = sum_{k=j+1..t} g_k * s[k-1-j]   (j = 0..t-1) */
    for (int j = 0; j < t; j++) {
        unsigned short acc = 0;
        for (int k = j + 1; k <= t; k++)
            acc ^= gf[mlt(fg[g[k]], fg[s[k - 1 - j]])];
        S[j] = acc;
    }
}


// Niederreiter暗号の公開鍵を作る
vec pubkeygen(unsigned short permutation[N],
              MTX *inverse_scrambler,
              RNG *rng,
              unsigned short syndrome_matrix[N][K])
{
  unsigned short (*field_matrix)[K] = calloc(N, sizeof(*field_matrix));
  MTX scrambler = {0};
  vec w = {0};
  const char *err = NULL;
  double t0 = now();

  if (field_matrix == NULL ||
      matrix_init(&scrambler, K * E, K * E) != 0 ||
      matrix_init(inverse_scrambler, K * E, K * E) != 0)
  {
    err = "pubkeygen: failed to allocate key-generation matrices";
    goto fail;
  }

  if (mkg(field_matrix, &w) != 0)
  {
    err = "pubkeygen: mkg failed";
    goto fail;
  }
  dbg_pol("==goppa polynomial", w);
  printf("Uh!-mkg: %.3f s\n", now() - t0);
  t0 = now();

#ifdef VERBOSE
  printsage(w);
  printf("\nsagemath で既約性を検査してください！\n");
#endif

  Pgen(permutation);
  printf("Uh!-Pgen: %.3f s\n", now() - t0);
  t0 = now();

  if (random_full_rank_binary_matrix_inv(&scrambler, inverse_scrambler, rng) != 0)
  {
    err = "failed to generate the scrambler matrix and its inverse";
    goto fail;
  }
  printf("Uh!-scrambler+inverse: %.3f s\n", now() - t0);
  t0 = now();

  if (pubkey_pack(&scrambler, field_matrix, permutation, syndrome_matrix) != 0)
  {
    err = "pubkeygen: failed to build the public key";
    goto fail;
  }
  printf("Uh!-pubkey_pack: %.3f s\n", now() - t0);

  free(field_matrix);
  matrix_free(&scrambler);
  return w;

fail:
  fprintf(stderr, "%s\n", err);
  free(field_matrix);
  matrix_free(&scrambler);
  matrix_free(inverse_scrambler);
  exit(EXIT_FAILURE);
}

// ===============================================================
// 暗号化・復号の周辺
// ===============================================================
static vec dec(const MTX *inverse_scrambler, unsigned short ss[])
{
  unsigned short uk[K];
  sinv_mul(inverse_scrambler, ss, uk);
  return Setvec(uk, K);
}

// シンドローム ss に S^{-1} を掛けて、多項式(vec)に直す
static vec dec2(const MTX *inverse_scrambler, unsigned short ss[])
{
  unsigned char ch[K * E];
  unsigned char h2o[K * E];
  unsigned short uk[K];

  // 13ビット整数 -> ビット列
  for (int i = 0; i < K; i++)
    for (int j = 0; j < E; j++)
      ch[i * E + j] = (ss[i] >> j) & 1;

  // GF(2) の行列・ベクトル積 (各行の AND の XOR 和 = パリティ)
  for (int i = 0; i < K * E; i++)
  {
    unsigned acc = 0;
    for (int j = 0; j < K * E; j++)
      acc ^= ch[j] & inverse_scrambler->x[i][j];
    h2o[i] = acc & 1;
  }

  // ビット列 -> 13ビット整数
  for (int i = 0; i < K; i++)
  {
    unsigned short u = 0;
    for (int j = 0; j < E; j++)
      u |= (unsigned short)(h2o[i * E + j] << j);
    uk[i] = u;
  }

  return Setvec(uk, K);
}

int ero2(vec v, const unsigned short permutation[N], unsigned short *ya)
{
  int count = 0;
  unsigned short xa[N] = {0};

  for (int i = 0; i < T * 2; i++)
  {
    if (i > 0 && v.x[i] == 0)
    {
      DBG("baka %d %d\n", i, v.x[i]);
      break;
    }
    xa[v.x[i]] = 1;
    count++;
  }

  for (int i = 0; i < N; i++)
    ya[i] = xa[permutation[i]];

#ifdef VERBOSE
  for (int i = 0; i < N; i++)
  {
    if (ya[i] > 0)
      printf("error position=%d\n", i);
  }
  if (count <= T * 2)
    printf("err=%d\n", count);
#endif

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
    DBG("l=%u %d\n", position, i);
  }

  return 0;
}

vec zin(const MTX *inverse_scrambler,
        const unsigned short syndrome_matrix[N][K],
        unsigned short zz[])
{
  unsigned short ss[K] = {0};

  for (int i = 0; i < N; i++)
  {
    if (zz[i] > 0)
    {
      for (int j = 0; j < K; j++)
        ss[j] ^= syndrome_matrix[i][j];
    }
  }

  // 暗号化されたシンドロームを復元する
  return dec(inverse_scrambler, ss);
}
