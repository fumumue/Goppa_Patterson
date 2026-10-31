#include "common.h"

#define F (K * E)

int random_full_rank_binary_matrix(MTX *matrix, RNG *rng)
{
  /*
   * 0/1 のスクランブラー行列 S を、GF(2) 上で必ずフルランクにする。
   *
   * 既存の実装では「乱数で行を作ってから独立性を判定」しており、
   * 一次従属になった場合にそのまま落ちていた。ここでは、
   * 一次従属の候補は捨てて、新しい候補を再生成し、
   * 独立な基底ベクトルを pivot 列ごとに積み上げる。
   *
   * これで最終的に F = K * E の独立な行が揃い、
   * rank F の正則行列が 1 回の生成ループで作れる。
   */
  unsigned char basis[F][F] = {0};
  bool used[F] = {false};
  memset(matrix->x, 0, sizeof(matrix->x));

  for (int row = 0; row < F; row++)
  {
    bool accepted = false;
    for (int attempt = 0; attempt < 4096 && !accepted; attempt++)
    {
      unsigned char candidate[F] = {0};
      for (int col = 0; col < F; col++)
        candidate[col] = (unsigned char)(xor128(rng) & 1UL);

      unsigned char reduced[F];
      memcpy(reduced, candidate, sizeof(reduced));

      for (int pivot = 0; pivot < F; pivot++)
      {
        if (!used[pivot])
          continue;
        if (reduced[pivot])
        {
          for (int col = pivot; col < F; col++)
            reduced[col] ^= basis[pivot][col];
        }
      }

      int pivot = -1;
      for (int col = 0; col < F; col++)
      {
        if (reduced[col] != 0)
        {
          pivot = col;
          break;
        }
      }

      if (pivot < 0)
        continue;

      memcpy(basis[pivot], reduced, sizeof(reduced));
      used[pivot] = true;

      for (int col = 0; col < F; col++)
        matrix->x[row][col] = reduced[col];

      accepted = true;
    }

    if (!accepted)
    {
      memset(matrix->x, 0, sizeof(matrix->x));
      return -1;
    }
  }

  return 0;
}

int is_reg(const MTX *cc_input, MTX *R, RNG *rng)
{
  (void)rng;
  MTX reduced = *cc_input;
  unsigned char inverse[F][F] = {0};

  memset(R, 0, sizeof(*R));

  for (int row = 0; row < F; row++)
  {
    for (int col = 0; col < F; col++)
    {
      if (reduced.x[row][col] > 1)
        return -1;
      inverse[row][col] = (unsigned char)(row == col);
    }
  }

  for (int col = 0; col < F; col++)
  {
    int pivot = col;
    while (pivot < F && reduced.x[pivot][col] == 0)
      pivot++;
    if (pivot == F)
      return -1;

    if (pivot != col)
    {
      for (int j = 0; j < F; j++)
      {
        unsigned short matrix_value = reduced.x[col][j];
        reduced.x[col][j] = reduced.x[pivot][j];
        reduced.x[pivot][j] = matrix_value;

        unsigned char inverse_value = inverse[col][j];
        inverse[col][j] = inverse[pivot][j];
        inverse[pivot][j] = inverse_value;
      }
    }

    for (int row = 0; row < F; row++)
    {
      if (row == col || reduced.x[row][col] == 0)
        continue;

      for (int j = 0; j < F; j++)
      {
        reduced.x[row][j] ^= reduced.x[col][j];
        inverse[row][j] ^= inverse[col][j];
      }
    }
  }

  for (int row = 0; row < F; row++)
    for (int col = 0; col < F; col++)
      R->x[row][col] = inverse[row][col];

  return 0;
}
