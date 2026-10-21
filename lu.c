//学校で習った知識を使ったら問題解決できました。大学行っといてよかったｗ

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "common.h"

#define F K *E
#define X (K) * E
#define AY (K / 2 + 1) * E

extern void makeS();

int random_full_rank_binary_matrix(MTX *matrix, RNG *rng)
{
  /*
   * この実装ではスクランブラー S は常に 0/1 の全順序行列として扱う。
   * 完全に一般的な乱数行列を逐次的に生成するよりも、
   * まず対角成分が 1 である上三角行列を作り、
   * それをランダム化する方が、GF(2) 上で必ず正則になる。
   *
   * ここでの対象サイズは F = K * E で、K=32, E=8 の場合 256x256。
   * したがって、対角成分 1 の上三角行列は常にフルランクであり、
   * 逆行列が求まるので、is_reg() の落ちを避けられる。
   */
  memset(matrix->x, 0, sizeof(matrix->x));

  for (int i = 0; i < F; i++)
  {
    matrix->x[i][i] = 1;
    for (int j = i + 1; j < F; j++)
      matrix->x[i][j] = (unsigned char)(xor128(rng) & 1UL);
  }

  return 0;
}

int is_reg(const MTX *cc_input, MTX *R, RNG *rng)
{
  MTX cc = *cc_input;

  int i, j, k, l;
  unsigned char b[K * E][K * E] = {0};
  unsigned int flg = 0, count = 0;
  unsigned char cl[K * E][K * E];
  time_t t;
  unsigned char inv_a[K * E][K * E] = {0}; //ここに逆行列が入る

    // memset(cc,0,sizeof(cc));
    flg = 0;
    count = 0;
    srand(clock() + time(&t));

    // g2();

    //#pragma omp parallel for private(j)
    for (i = 0; i < F; i++)
    {

      for (j = 0; j < F; j++)
      {
        // printf("%d,",cc.x[i][j]);
        cl[i][j] = cc.x[i][j];
      }
      // printf("\n");
    }

    // memset(inv_a,0,sizeof(inv_a));

    //単位行列を作る
    //#pragma omp parallel for
    for (i = 0; i < F; i++)
    {
      for (j = 0; j < F; j++)
      {
        inv_a[i][j] = (i == j) ? 1.0 : 0.0;
      }
    }

    //掃き出し法

    for (i = 0; i < F; i++)
    {
      if (cc.x[i][i] == 0)
      {
        j = i;

        while (j < F && cc.x[j][i] == 0)
        {
          j++;
        }

        //#pragma omp parallel for
        if (j >= F)
        {
          printf("S is not reg in is_reg %d\n", j);
          return -1;
        }
        xor_u16_mod2(cc.x[i], cc.x[j], F);
        xor_u8(inv_a[i], inv_a[j], F);

        cc.x[i][i] = 1;
      }

      if (cc.x[i][i] == 1)
      {
        for (l = i + 1; l < F; l++)
        {
          if (cc.x[l][i] == 1)
          {
            xor_u16_mod2(cc.x[l], cc.x[i], F);
            xor_u8(inv_a[l], inv_a[i], F);
          }
        }
      }
    }

    for (i = 1; i < F; i++)
    {
      for (k = 0; k < i; k++)
      {
        if (cc.x[k][i] == 1)
        {
          xor_u16_mod2(cc.x[k], cc.x[i], F);
          xor_u8(inv_a[k], inv_a[i], F);
        }
      }
    }

    //検算
    for (i = 0; i < F; i++)
    {
      for (j = 0; j < F; j++)
      {
        for (k = 0; k < F; k++)
        {
          b[i][j] ^= (cl[i][k] & inv_a[k][j]);
        }
      }
    }

    for (i = 0; i < F; i++)
    {
      if (b[i][i] == 1)
      {
        flg++;
      }
    }
    count = 0;

    for (i = 0; i < F; i++)
    {
      for (j = 0; j < F; j++)
      {
        if (b[i][j] == 0 && i != j)
          count++;
      }
    }
    printf("%d,%d %d,%d", flg, F, count, F * F - F);
    printf("S[K][K]=\n{\n");
    if (flg == F && count == (F * F - F))
    {
      printf("inv_S[K][K]=\n{\n");
      for (i = 0; i < F; i++)
      {
        printf("{");
        for (j = 0; j < F; j++)
        {
          R->x[i][j] = inv_a[i][j];
          printf("%d,", inv_a[i][j]);
        }
        printf("},\n");
      }
      printf("};\n");
      return 0;
    }
    return -1;
}

int mkS(const MTX *cc_input, MTX *R)
{
  MTX cc = *cc_input;
  int i, j, k, l;
  unsigned char b[X][X] = {0};
  unsigned int flg = 0, count = 0;
  // unsigned char cc.x[X][X] = {0};
  unsigned char cl[X][X];
  time_t t;
  unsigned char inv_a[X][X] = {0}; //ここに逆行列が入る

    flg = 0;
    count = 0;
    srand(clock() + time(&t));

    for (i = 0; i < X; i++)
    {
      for (j = 0; j < X; j++)
      {
        cl[i][j] = cc.x[i][j];
      }
    }

    for (i = 0; i < X; i++)
    {
      for (j = 0; j < X; j++)
      {
        inv_a[i][j] = (i == j) ? 1.0 : 0.0;
      }
    }

    for (i = 0; i < X; i++)
    {
      if (cc.x[i][i] == 0)
      {
        j = i;

        while (j < X && cc.x[j][i] == 0)
        {
          j++;
        }

        if (j >= X)
        {
          printf("baka in mkS %d\n", j);
          return -1;
        }
        xor_u16_mod2(cc.x[i], cc.x[j], X);
        xor_u8(inv_a[i], inv_a[j], X);

        cc.x[i][i] = 1;
      }

      if (cc.x[i][i] == 1)
      {
        for (l = i + 1; l < X; l++)
        {
          if (cc.x[l][i] == 1)
          {
            xor_u16_mod2(cc.x[l], cc.x[i], X);
            xor_u8(inv_a[l], inv_a[i], X);
          }
        }
      }
    }

    for (i = 1; i < X; i++)
    {
      for (k = 0; k < i; k++)
      {
        if (cc.x[k][i] == 1)
        {
          xor_u16_mod2(cc.x[k], cc.x[i], X);
          xor_u8(inv_a[k], inv_a[i], X);
        }
      }
    }

    for (i = 0; i < X; i++)
    {
      for (j = 0; j < X; j++)
      {
        for (k = 0; k < X; k++)
        {
          b[i][j] ^= (cl[i][k] & inv_a[k][j]);
        }
      }
    }

    for (i = 0; i < X; i++)
    {
      if (b[i][i] == 1)
      {
        flg++;
      }
    }
    count = 0;

    for (i = 0; i < X; i++)
    {
      for (j = 0; j < X; j++)
      {
        if (b[i][j] == 0 && i != j)
          count++;
      }
    }

    printf("S[K][K]=\n{\n");
    if (flg == X && count == (X * X - X))
    {
      for (i = 0; i < X; i++)
      {
        for (j = 0; j < X; j++)
        {
          R->x[i][j] = inv_a[i][j];
        }
      }
      return 0;
    }
    return -1;
}

int mkRS(const MTX *cc_input, MTX *R)
{
  MTX cc = *cc_input;
  int i, j, k, l;
  unsigned char b[AY][AY] = {0};
  unsigned int flg = 0, count = 0;
  unsigned char cl[AY][AY];
  time_t t;
  unsigned char inv_a[AY][AY] = {0};

    flg = 0;
    count = 0;
    srand(clock() + time(&t));

    for (i = 0; i < AY; i++)
    {
      for (j = 0; j < AY; j++)
      {
        cl[i][j] = cc.x[i][j];
      }
    }

    for (i = 0; i < AY; i++)
    {
      for (j = 0; j < AY; j++)
      {
        inv_a[i][j] = (i == j) ? 1.0 : 0.0;
      }
    }

    for (i = 0; i < AY; i++)
    {
      if (cc.x[i][i] == 0)
      {
        j = i;

        while (j < AY && cc.x[j][i] == 0)
        {
          j++;
        }

        if (j >= AY)
        {
          printf("baka in mkS %d\n", j);
          return -1;
        }
        xor_u16_mod2(cc.x[i], cc.x[j], AY);
        xor_u8(inv_a[i], inv_a[j], AY);

        cc.x[i][i] = 1;
      }

      if (cc.x[i][i] == 1)
      {
        for (l = i + 1; l < AY; l++)
        {
          if (cc.x[l][i] == 1)
          {
            xor_u16_mod2(cc.x[l], cc.x[i], AY);
            xor_u8(inv_a[l], inv_a[i], AY);
          }
        }
      }
    }

    for (i = 1; i < AY; i++)
    {
      for (k = 0; k < i; k++)
      {
        if (cc.x[k][i] == 1)
        {
          xor_u16_mod2(cc.x[k], cc.x[i], AY);
          xor_u8(inv_a[k], inv_a[i], AY);
        }
      }
    }

    for (i = 0; i < AY; i++)
    {
      for (j = 0; j < AY; j++)
      {
        for (k = 0; k < AY; k++)
        {
          b[i][j] ^= (cl[i][k] & inv_a[k][j]);
        }
      }
    }

    for (i = 0; i < AY; i++)
    {
      if (b[i][i] == 1)
      {
        flg++;
      }
    }
    count = 0;

    for (i = 0; i < AY; i++)
    {
      for (j = 0; j < AY; j++)
      {
        if (b[i][j] == 0 && i != j)
          count++;
      }
    }

    printf("S[K][K]=\n{\n");
    if (flg == AY && count == (AY * AY - AY))
    {
      for (i = 0; i < AY; i++)
      {
        for (j = 0; j < AY; j++)
        {
          R->x[i][j] = inv_a[i][j];
        }
      }
      return 0;
    }
    return -1;
}

int binv(const MTX *cc_input, MTX *L, int Y)
{
  MTX cc = *cc_input;
  int i, j, k, l;
  unsigned char b[N][N] = {0};
  unsigned int flg = 0, count = 0;
  unsigned char cl[N][N];
  time_t t;
  unsigned char inv_a[N][N] = {0};

    flg = 0;
    count = 0;
    srand(clock() + time(&t));

    for (i = 0; i < Y; i++)
    {
      for (j = 0; j < Y; j++)
      {
        cl[i][j] = cc.x[i][j];
      }
    }

    for (i = 0; i < Y; i++)
    {
      for (j = 0; j < Y; j++)
      {
        inv_a[i][j] = (i == j) ? 1.0 : 0.0;
      }
    }

    for (i = 0; i < Y; i++)
    {
      if (cc.x[i][i] == 0)
      {
        j = i;

        while (j < Y && cc.x[j][i] == 0)
        {
          j++;
        }

        if (j >= Y)
        {
          printf("baka in binv %d\n", j);
          return -1;
        }
        xor_u16_mod2(cc.x[i], cc.x[j], Y);
        xor_u8(inv_a[i], inv_a[j], Y);

        cc.x[i][i] = 1;
      }

      if (cc.x[i][i] == 1)
      {
        for (l = i + 1; l < Y; l++)
        {
          if (cc.x[l][i] == 1)
          {
            xor_u16_mod2(cc.x[l], cc.x[i], Y);
            xor_u8(inv_a[l], inv_a[i], Y);
          }
        }
      }
    }

    for (i = 1; i < Y; i++)
    {
      for (k = 0; k < i; k++)
      {
        if (cc.x[k][i] == 1)
        {
          xor_u16_mod2(cc.x[k], cc.x[i], Y);
          xor_u8(inv_a[k], inv_a[i], Y);
        }
      }
    }

#pragma omp parallel for private(j, k) num_threads(16)
    for (i = 0; i < Y; i++)
    {
      for (j = 0; j < Y; j++)
      {
        for (k = 0; k < Y; k++)
        {
          b[i][j] ^= (cl[i][k] & inv_a[k][j]);
        }
      }
    }

    for (i = 0; i < Y; i++)
    {
      if (b[i][i] == 1)
      {
        flg++;
      }
    }
    count = 0;

    for (i = 0; i < Y; i++)
    {
      for (j = 0; j < Y; j++)
      {
        if (b[i][j] == 0 && i != j)
          count++;
      }
    }

    printf("S[K][K]=\n{\n");
    if (flg == Y && count == (Y * Y - Y))
    {
      for (i = 0; i < Y; i++)
      {
        for (j = 0; j < Y; j++)
        {
          L->x[i][j] = inv_a[i][j];
        }
      }
      return 0;
    }
    return -1;
}
