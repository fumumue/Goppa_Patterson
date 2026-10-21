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
