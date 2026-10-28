#include "common.h"

static unsigned char packed_bit(const uint64_t *row, size_t bit)
{
  return (unsigned char)((row[bit / 64] >> (bit % 64)) & 1U);
}

int matrix_init(MTX *matrix, size_t rows, size_t cols)
{
  if (matrix == NULL || rows == 0 || cols == 0 ||
      rows > SIZE_MAX / cols ||
      rows * cols > SIZE_MAX / sizeof(*matrix->storage))
    return -1;

  MTX allocated = {0};
  allocated.x = calloc(rows, sizeof(*allocated.x));
  allocated.storage = calloc(rows * cols, sizeof(*allocated.storage));
  if (allocated.x == NULL || allocated.storage == NULL)
  {
    free(allocated.x);
    free(allocated.storage);
    return -1;
  }

  allocated.rows = rows;
  allocated.cols = cols;
  for (size_t row = 0; row < rows; row++)
    allocated.x[row] = allocated.storage + row * cols;

  matrix_free(matrix);
  *matrix = allocated;
  return 0;
}

void matrix_free(MTX *matrix)
{
  if (matrix == NULL)
    return;
  free(matrix->x);
  free(matrix->storage);
  memset(matrix, 0, sizeof(*matrix));
}

int matrix_copy(MTX *destination, const MTX *source)
{
  if (destination == NULL || source == NULL ||
      source->x == NULL || source->storage == NULL)
    return -1;
  if (destination == source)
    return 0;
  if (matrix_init(destination, source->rows, source->cols) != 0)
    return -1;
  memcpy(destination->storage, source->storage,
         source->rows * source->cols * sizeof(*source->storage));
  return 0;
}


#ifndef INV_WB
#define INV_WB 4   /* 後退代入の語ブロック幅。4と8を比べてください */
#endif

/* シードから64ビットずつ乱数を作る(速度検証用。暗号用ではありません) */
static inline uint64_t splitmix64(uint64_t *s)
{
  uint64_t z = (*s += UINT64_C(0x9E3779B97F4A7C15));
  z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
  return z ^ (z >> 31);
}

// S と S^{-1} を同時に作る。従属な行はその行だけ作り直す。
int random_full_rank_binary_matrix_inv(MTX *S, MTX *Sinv, RNG *rng)
{
  if (S == NULL || Sinv == NULL || rng == NULL ||
      S->x == NULL || Sinv->x == NULL ||
      S->storage == NULL || Sinv->storage == NULL ||
      S->rows == 0 || S->rows != S->cols ||
      Sinv->rows != S->rows || Sinv->cols != S->cols)
  {
    fprintf(stderr, "invalid matrix arguments for scrambler generation\n");
    return -1;
  }

  const size_t n = S->rows;
  const size_t words = (n + 63) / 64;
  if (n > SIZE_MAX / words || n * words > SIZE_MAX / sizeof(uint64_t))
    return -1;
  const uint64_t lastmask = (n % 64) ? ((UINT64_C(1) << (n % 64)) - 1)
                                     : ~UINT64_C(0);

  uint64_t *basis = calloc(n * words, sizeof *basis); // basis[p]: 最下位ビットが p の行
  uint64_t *track = calloc(n * words, sizeof *track); // basis[p] = track[p] に立つ S の行の和
  uint64_t *inv   = calloc(n * words, sizeof *inv);   // S^{-1} (ビットパック)
  uint64_t *spack = calloc(n * words, sizeof *spack); // 採用した S の行(ビットパック)
  uint64_t *umask = calloc(words, sizeof *umask);     // ピボットとして使用済みのビット
  uint64_t *cand  = malloc(words * sizeof *cand);
  uint64_t *red   = malloc(words * sizeof *red);
  uint64_t *trk   = malloc(words * sizeof *trk);
  int rc = -1;

  /* 元のコードと同じ「最下位ビット」から64ビットのシードを作る */
  uint64_t sm = 0;
  for (int i = 0; i < 64; i++)
    sm = (sm << 1) | (uint64_t)(xor128(rng) & 1UL);

  if (!basis || !track || !inv || !spack || !umask || !cand || !red || !trk)
  {
    fprintf(stderr, "failed to allocate scrambler generation workspace\n");
    goto done;
  }

  for (size_t row = 0; row < n; row++)
  {
    const size_t trk_words = row / 64 + 1;   // 追跡ベクトルは 0..row のビットだけ使う
    bool accepted = false;

    for (int attempt = 0; attempt < 4096 && !accepted; attempt++)
    {
#ifdef SCR_RNG_ORIG
      /* 切り分け用: 元のコードと同じ 1 ビットずつの生成 */
      memset(cand, 0, words * sizeof *cand);
      for (size_t col = 0; col < n; col++)
        cand[col / 64] |= (uint64_t)(xor128(rng) & 1UL) << (col % 64);
#else
      for (size_t w = 0; w < words; w++)
        cand[w] = splitmix64(&sm);
      cand[words - 1] &= lastmask;
#endif

      memcpy(red, cand, words * sizeof *red);
      memset(trk, 0, trk_words * sizeof *trk);
      trk[row / 64] |= UINT64_C(1) << (row % 64);

      // 立っているピボットビットを下位から消す(basis[p] は p より下が0)
      for (size_t w = 0; w < words; w++)
      {
        uint64_t bits;
        while ((bits = red[w] & umask[w]) != 0)
        {
          const size_t p = w * 64 + (size_t)__builtin_ctzll(bits);
          const uint64_t *b = basis + p * words;
          const uint64_t *t = track + p * words;
          for (size_t x = w; x < words; x++)
            red[x] ^= b[x];
          for (size_t x = 0; x < trk_words; x++)
            trk[x] ^= t[x];
        }
      }

      size_t q = n;
      for (size_t w = 0; w < words; w++)
        if (red[w] != 0)
        {
          q = w * 64 + (size_t)__builtin_ctzll(red[w]);
          break;
        }
      if (q == n)
        continue;                             // 従属: この行だけ作り直す

      memcpy(basis + q * words, red, words * sizeof *red);
      memcpy(track + q * words, trk, trk_words * sizeof *trk);
      umask[q / 64] |= UINT64_C(1) << (q % 64);
      memcpy(spack + row * words, cand, words * sizeof *cand);
      accepted = true;
    }

    if (!accepted)
    {
      fprintf(stderr, "scrambler: no independent row at row %zu of %zu\n", row, n);
      goto done;
    }
  }

  // B = T·S  (B は対角が1の上三角)。S^{-1} = B^{-1}·T を後退代入で求める。
  //   inv[p] = track[p] ^ XOR_{q>p, B[p][q]=1} inv[q]
  // 語ごとに独立なので、語のブロック単位で並列化する。
  {
    const size_t WB = INV_WB;
    const size_t nblk = (words + WB - 1) / WB;

#pragma omp parallel for schedule(static)
    for (size_t blk = 0; blk < nblk; blk++)
    {
      const size_t w0 = blk * WB;
      const size_t w1 = (w0 + WB < words) ? w0 + WB : words;

      for (size_t p = n; p-- > 0;)
      {
        uint64_t *xp = inv + p * words;
        const uint64_t *bp = basis + p * words;
        for (size_t w = w0; w < w1; w++)
          xp[w] = track[p * words + w];

        for (size_t wi = p / 64; wi < words; wi++)
        {
          uint64_t bits = bp[wi];
          if (wi == p / 64)
            bits &= ~((UINT64_C(2) << (p % 64)) - 1);   // q > p だけ残す
          while (bits)
          {
            size_t q = wi * 64 + (size_t)__builtin_ctzll(bits);
            bits &= bits - 1;
            const uint64_t *xq = inv + q * words;
            for (size_t w = w0; w < w1; w++)
              xp[w] ^= xq[w];
          }
        }
      }
    }
  }

  // S と S^{-1} を 1 バイト形式へ展開(行ごとに独立なので並列)
#pragma omp parallel for schedule(static)
  for (size_t r = 0; r < n; r++)
  {
    const uint64_t *sp = spack + r * words;
    const uint64_t *ip = inv + r * words;
    for (size_t w = 0; w < words; w++)
    {
      const size_t c0 = w * 64;
      const size_t cn = (n - c0 < 64) ? (n - c0) : 64;
      const uint64_t sv = sp[w], iv = ip[w];
      for (size_t b = 0; b < cn; b++)
      {
        S->x[r][c0 + b]    = (sv >> b) & 1;
        Sinv->x[r][c0 + b] = (iv >> b) & 1;
      }
    }
  }
  rc = 0;

done:
  if (rc != 0 && S != NULL && S->storage != NULL)
    memset(S->storage, 0, S->rows * S->cols * sizeof(*S->storage));
  free(basis); free(track); free(inv); free(spack); free(umask);
  free(cand);  free(red);   free(trk);
  return rc;
}
// S と S^{-1} を同時に作る。従属な行はその行だけ作り直す。
int random_full_rank_binary_matrix_inv2(MTX *S, MTX *Sinv, RNG *rng)
{ 
  if (S == NULL || Sinv == NULL || rng == NULL ||
      S->x == NULL || Sinv->x == NULL ||
      S->storage == NULL || Sinv->storage == NULL ||
      S->rows == 0 || S->rows != S->cols ||
      Sinv->rows != S->rows || Sinv->cols != S->cols)
  {
    fprintf(stderr, "invalid matrix arguments for scrambler generation\n");
    return -1;
  }

  const size_t n = S->rows;
  const size_t words = (n + 63) / 64;
  if (n > SIZE_MAX / words || n * words > SIZE_MAX / sizeof(uint64_t))
    return -1;

  uint64_t *basis = calloc(n * words, sizeof *basis); // basis[p]: 最下位ビットが p の行
  uint64_t *track = calloc(n * words, sizeof *track); // basis[p] = track[p] に立つ S の行の和
  uint64_t *inv   = calloc(n * words, sizeof *inv);   // S^{-1} (ビットパック)
  uint64_t *cand  = malloc(words * sizeof *cand);
  uint64_t *red   = malloc(words * sizeof *red);
  uint64_t *trk   = malloc(words * sizeof *trk);
  unsigned char *used = calloc(n, 1);
  int rc = -1;

  if (!basis || !track || !inv || !cand || !red || !trk || !used)
  {
    fprintf(stderr, "failed to allocate scrambler generation workspace\n");
    goto done;
  }

  for (size_t row = 0; row < n; row++)
  {
    const size_t trk_words = row / 64 + 1;   // 追跡ベクトルは 0..row のビットだけ使う
    bool accepted = false;

    for (int attempt = 0; attempt < 4096 && !accepted; attempt++)
    {
      memset(cand, 0, words * sizeof *cand);
      for (size_t col = 0; col < n; col++)
        cand[col / 64] |= (uint64_t)(xor128(rng) & 1UL) << (col % 64);

      memcpy(red, cand, words * sizeof *red);
      memset(trk, 0, trk_words * sizeof *trk);
      trk[row / 64] |= UINT64_C(1) << (row % 64);

      for (size_t p = 0; p < n; p++)
      {
        if (!used[p] || !packed_bit(red, p))
          continue;
        const uint64_t *b = basis + p * words;
        const uint64_t *t = track + p * words;
        for (size_t w = p / 64; w < words; w++)
          red[w] ^= b[w];
        for (size_t w = 0; w < trk_words; w++)
          trk[w] ^= t[w];
      }

      size_t q = n;
      for (size_t w = 0; w < words; w++)
        if (red[w] != 0)
        {
          q = w * 64 + (size_t)__builtin_ctzll(red[w]);
          break;
        }
      if (q == n)
        continue;                             // 従属: この行だけ作り直す

      memcpy(basis + q * words, red, words * sizeof *red);
      memcpy(track + q * words, trk, trk_words * sizeof *trk);
      used[q] = 1;

      for (size_t col = 0; col < n; col++)
        S->x[row][col] = (cand[col / 64] >> (col % 64)) & 1;
      accepted = true;
    }

    if (!accepted)
    {
      fprintf(stderr, "scrambler: no independent row at row %zu of %zu\n", row, n);
      goto done;
    }
  }

  // B = T·S  (B は対角が1の上三角)。S^{-1} = B^{-1}·T を後退代入で求める。
  //   inv[p] = track[p] ^ XOR_{q>p, B[p][q]=1} inv[q]
  // 語ごとに独立なので、語のブロック単位で並列化する。
  {
    const size_t WB = 8;
    const size_t nblk = (words + WB - 1) / WB;

#pragma omp parallel for schedule(static)
    for (size_t blk = 0; blk < nblk; blk++)
    {
      const size_t w0 = blk * WB;
      const size_t w1 = (w0 + WB < words) ? w0 + WB : words;

      for (size_t p = n; p-- > 0;)
      {
        uint64_t *xp = inv + p * words;
        const uint64_t *bp = basis + p * words;
        for (size_t w = w0; w < w1; w++)
          xp[w] = track[p * words + w];

        for (size_t wi = p / 64; wi < words; wi++)
        {
          uint64_t bits = bp[wi];
          if (wi == p / 64)
            bits &= ~((UINT64_C(2) << (p % 64)) - 1);   // q > p だけ残す
          while (bits)
          {
            size_t q = wi * 64 + (size_t)__builtin_ctzll(bits);
            bits &= bits - 1;
            const uint64_t *xq = inv + q * words;
            for (size_t w = w0; w < w1; w++)
              xp[w] ^= xq[w];
          }
        }
      }
    }
  }

  for (size_t r = 0; r < n; r++)
    for (size_t c = 0; c < n; c++)
      Sinv->x[r][c] = (inv[r * words + c / 64] >> (c % 64)) & 1;
  rc = 0;

done:
  if (rc != 0 && S != NULL && S->storage != NULL)
    memset(S->storage, 0, S->rows * S->cols * sizeof(*S->storage));
  free(basis); free(track); free(inv);
  free(cand);  free(red);   free(trk);
  free(used);
  return rc;
}


int random_full_rank_binary_matrix(MTX *matrix, RNG *rng)
{
  /* Keep one reduced basis row per pivot, packed to avoid cubic bit loops. */
  if (matrix == NULL || matrix->x == NULL || matrix->storage == NULL ||
      matrix->rows == 0 || matrix->rows != matrix->cols || rng == NULL)
  {
    fprintf(stderr, "invalid or unavailable storage for binary scrambler matrix\n");
    return -1;
  }
  size_t n = matrix->rows;
  size_t words = (n + 63) / 64;
  if (n > SIZE_MAX / words || n * words > SIZE_MAX / sizeof(uint64_t))
    return -1;
  uint64_t *basis = calloc(n * words, sizeof(*basis));
  unsigned char *used = calloc(n, sizeof(*used));
  uint64_t *candidate = malloc(words * sizeof(*candidate));
  uint64_t *reduced = malloc(words * sizeof(*reduced));
  if (basis == NULL || used == NULL || candidate == NULL || reduced == NULL)
  {
    free(basis);
    free(used);
    free(candidate);
    free(reduced);
    fprintf(stderr, "failed to allocate binary scrambler generation workspace\n");
    return -1;
  }
  memset(matrix->storage, 0, matrix->rows * matrix->cols * sizeof(*matrix->storage));

  for (size_t row = 0; row < n; row++)
  {
    bool accepted = false;
    for (int attempt = 0; attempt < 4096 && !accepted; attempt++)
    {
      memset(candidate, 0, words * sizeof(*candidate));
      for (size_t w = 0; w < words; w++)
        candidate[w] = (uint64_t)xor128(rng);
      if (n % 64)
        candidate[words - 1] &= (UINT64_C(1) << (n % 64)) - 1;
      for (size_t col = 0; col < n; col++)
        candidate[col / 64] |= (uint64_t)(xor128(rng) & 1UL) << (col % 64);
      memcpy(reduced, candidate, words * sizeof(*reduced));

      for (size_t pivot = 0; pivot < n; pivot++)
      {
        if (!used[pivot])
          continue;
        if (packed_bit(reduced, pivot))
        {
          const uint64_t *basis_row = basis + pivot * words;
          for (size_t word = pivot / 64; word < words; word++)
            reduced[word] ^= basis_row[word];
        }
      }

      size_t pivot = n;
      for (size_t word = 0; word < words; word++)
      {
        uint64_t bits = reduced[word];
        if (bits != 0)
        {
          pivot = word * 64 + (size_t)__builtin_ctzll(bits);
          break;
        }
      }

      if (pivot == n)
        continue;

      memcpy(basis + pivot * words, reduced, words * sizeof(*reduced));
      used[pivot] = true;

      for (size_t col = 0; col < n; col++)
        matrix->x[row][col] = packed_bit(reduced, col);

      accepted = true;
    }

    if (!accepted)
    {
      memset(matrix->storage, 0,
             matrix->rows * matrix->cols * sizeof(*matrix->storage));
      free(basis);
      free(used);
      free(candidate);
      free(reduced);
      return -1;
    }
  }

  free(basis);
  free(used);
  free(candidate);
  free(reduced);
  return 0;
}


static int check_inverse(const MTX *S, const MTX *Sinv, RNG *rng)
{
  size_t n = S->rows;
  unsigned short *v = malloc(n * sizeof *v), *y = malloc(n * sizeof *y);
  int ok = (v && y);
  for (int t = 0; ok && t < 8; t++)
  {
    for (size_t i = 0; i < n; i++) v[i] = xor128(rng) & 1UL;
    for (size_t i = 0; i < n; i++)
    {
      unsigned short s = 0;
      for (size_t j = 0; j < n; j++) s ^= Sinv->x[i][j] & v[j];
      y[i] = s;
    }
    for (size_t i = 0; i < n; i++)
    {
      unsigned short s = 0;
      for (size_t j = 0; j < n; j++) s ^= S->x[i][j] & y[j];
      if (s != v[i]) { ok = 0; break; }
    }
  }
  free(v); free(y);
  return ok ? 0 : -1;
}


int is_reg(const MTX *cc_input, MTX *R, RNG *rng)
{
  (void)rng;
  if (cc_input == NULL || R == NULL || cc_input->x == NULL ||
      cc_input->rows == 0 || cc_input->rows != cc_input->cols ||
      R->x == NULL || R->rows != cc_input->rows ||
      R->cols != cc_input->cols)
  {
    fprintf(stderr, "invalid matrix dimensions for binary matrix inversion\n");
    return -1;
  }
  size_t n = cc_input->rows;
  size_t words = (n + 63) / 64;
  if (n > SIZE_MAX / words || n * words > SIZE_MAX / sizeof(uint64_t))
    return -1;
  uint64_t *reduced = calloc(n * words, sizeof(*reduced));
  uint64_t *inverse = calloc(n * words, sizeof(*inverse));
  if (reduced == NULL || inverse == NULL)
  {
    free(reduced);
    free(inverse);
    fprintf(stderr, "failed to allocate binary matrix inversion workspace\n");
    return -1;
  }
  memset(R->storage, 0, R->rows * R->cols * sizeof(*R->storage));

  for (size_t row = 0; row < n; row++)
  {
    for (size_t col = 0; col < n; col++)
    {
      unsigned short value = cc_input->x[row][col];
      if (value > 1)
      {
        free(reduced);
        free(inverse);
        return -1;
      }
      if (value != 0)
        reduced[row * words + col / 64] |= UINT64_C(1) << (col % 64);
    }
    inverse[row * words + row / 64] |= UINT64_C(1) << (row % 64);
  }

  for (size_t col = 0; col < n; col++)
  {
    size_t pivot = col;
    while (pivot < n && !packed_bit(reduced + pivot * words, col))
      pivot++;
    if (pivot == n)
    {
      free(reduced);
      free(inverse);
      return -1;
    }

    if (pivot != col)
    {
      for (size_t word = 0; word < words; word++)
      {
        uint64_t value = reduced[col * words + word];
        reduced[col * words + word] = reduced[pivot * words + word];
        reduced[pivot * words + word] = value;
        value = inverse[col * words + word];
        inverse[col * words + word] = inverse[pivot * words + word];
        inverse[pivot * words + word] = value;
      }
    }

    
    //#pragma omp parallel for schedule(static)
    for (size_t row = 0; row < n; row++)
    {
      if (row == col || !packed_bit(reduced + row * words, col))
        continue;

      for (size_t word = col / 64; word < words; word++)
        reduced[row * words + word] ^= reduced[col * words + word];
      for (size_t word = 0; word < words; word++)
        inverse[row * words + word] ^= inverse[col * words + word];
      /*
      for (size_t word = 0; word < words; word++)
      {
        reduced[row * words + word] ^= reduced[col * words + word];
        inverse[row * words + word] ^= inverse[col * words + word];
      }
      */
    }
  }

  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      R->x[row][col] = packed_bit(inverse + row * words, col);

  free(reduced);
  free(inverse);
  return 0;
}
