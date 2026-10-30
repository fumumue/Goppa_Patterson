#include "sinv_pack.h"

#include <stdlib.h>
#include <string.h>

#define NB (K * E)          /* 行列の一辺 (1664) */
#define NW ((NB + 63) / 64) /* 1 行あたりの 64 ビット語の数 (26) */

typedef struct
{
  const MTX *owner;              /* どの MTX の写しか */
  const unsigned short *storage; /* owner->storage (差し替えの検出用) */
  uint64_t *w;                   /* NB 行 x NW 語 */
} SinvCache;

static SinvCache g;

void sinv_pack_invalidate(void)
{
  free(g.w);
  memset(&g, 0, sizeof g);
}

void sinv_pack_register(const MTX *Sinv, const uint64_t *inv, size_t n, size_t words)
{
  sinv_pack_invalidate();
  if (Sinv == NULL || inv == NULL || n != (size_t)NB || words != (size_t)NW)
    return;

  uint64_t *p = malloc((size_t)NB * NW * sizeof *p);
  if (p == NULL)
    return;
  memcpy(p, inv, (size_t)NB * NW * sizeof *p);

  g.owner = Sinv;
  g.storage = Sinv->storage;
  g.w = p;
}

/* 登録した写しが、いま渡された MTX と同じものかを確かめる。
 * ポインタの一致に加えて、固定の64か所を実際の MTX と突き合わせる
 * (別の鍵に差し替わっていれば、ほぼ確実にここで見つかる)。 */
static int cache_valid(const MTX *m)
{
  if (g.w == NULL || m == NULL || g.owner != m || g.storage != m->storage ||
      m->rows != (size_t)NB || m->cols != (size_t)NB)
    return 0;

  uint64_t s = UINT64_C(0x9E3779B97F4A7C15);
  for (int k = 0; k < 64; k++)
  {
    s = s * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
    const size_t r = (size_t)((s >> 33) % NB);
    const size_t c = (size_t)((s >> 11) % NB);
    const unsigned bit = (unsigned)((g.w[r * NW + c / 64] >> (c % 64)) & 1);
    if (bit != (m->x[r][c] != 0))
      return 0;
  }
  return 1;
}

/* 従来と同じ計算 (1ビット=1要素。パック版が使えないときの保険) */
static void mul_slow(const MTX *m, const unsigned short ss[K], unsigned short uk[K])
{
  unsigned char ch[NB], h2o[NB];

  for (int i = 0; i < K; i++)
    for (int j = 0; j < E; j++)
      ch[i * E + j] = (ss[i] >> j) & 1;

  for (int i = 0; i < NB; i++)
  {
    unsigned acc = 0;
    for (int j = 0; j < NB; j++)
      acc ^= ch[j] & m->x[i][j];
    h2o[i] = acc & 1;
  }

  for (int i = 0; i < K; i++)
  {
    unsigned short u = 0;
    for (int j = 0; j < E; j++)
      u |= (unsigned short)(h2o[i * E + j] << j);
    uk[i] = u;
  }
}

/* ビットパック版: 各行 26 語との AND を XOR でまとめ、パリティを取る */
static void mul_packed(const uint64_t *w, const unsigned short ss[K], unsigned short uk[K])
{
  uint64_t chp[NW] = {0};

  for (int i = 0; i < K; i++)
  {
    const unsigned p = (unsigned)i * E;
    const uint64_t v = (uint64_t)(ss[i] & ((1u << E) - 1));
    chp[p >> 6] |= v << (p & 63);
    if ((p & 63) + E > 64)
      chp[(p >> 6) + 1] |= v >> (64 - (p & 63));
  }

  for (int k = 0; k < K; k++)
  {
    unsigned short u = 0;
    for (int j = 0; j < E; j++)
    {
      const uint64_t *row = w + (size_t)(k * E + j) * NW;
      uint64_t acc = 0;
      for (size_t x = 0; x < NW; x++)
        acc ^= row[x] & chp[x];
      u |= (unsigned short)(__builtin_parityll(acc) << j);
    }
    uk[k] = u;
  }
}

int sinv_mul(const MTX *Sinv, const unsigned short ss[K], unsigned short uk[K])
{
  if (cache_valid(Sinv))
  {
    mul_packed(g.w, ss, uk);
    return 1;
  }
  mul_slow(Sinv, ss, uk);
  return 0;
}
