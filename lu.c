// 学校で習った知識を使ったら問題解決できました。大学行っといてよかったｗ

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define F K *E
#define X (K) * E
#define AY (K / 2 + 1) * E
#define MATRIX_SIZE 16

// GF(2)上の元を表すデータ型
typedef unsigned short gf2;


static MTX inv_S = {0};
static MTX S = {0};
static MTX SS = {0};

extern void makeS();

// GF(2)上の乗算
unsigned short gf2_mul(unsigned short a, unsigned short b) {
    unsigned short result = 0;
    while (b) {
        if (b & 1) {
            result ^= a;
        }
        b >>= 1;
        a <<= 1;
    }
    return result;
}

// GF(2)上での単位行列の生成
void gf2_identity(gf2 matrix[MATRIX_SIZE][MATRIX_SIZE]) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            matrix[i][j] = (i == j) ? 1 : 0;
        }
    }
}

// 行列式を計算する関数（GF(2)上では単純に対角線上の要素の積）
unsigned short determinant(unsigned short matrix[MATRIX_SIZE][MATRIX_SIZE]) {
    unsigned short det = 1;
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        det = gf2_mul(det, matrix[i][i]); // 対角線上の要素の積
    }
    return det;
}


// GF(2)上でランダムな行列を生成する関数
void generate_random_regular_matrix(gf2 matrix[MATRIX_SIZE][MATRIX_SIZE]) {
    srand(time(NULL));
    gf2_identity(matrix); // 単位行列を初期化

    // ランダムな正則行列を生成する
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            matrix[i][j] = rand() % 2; // 0または1のランダムな値を設定
        }
    }

    // 正則性の確認
    // 逆行列が存在しない場合は、再度行列を生成する
    while (1) {
        // 生成した行列が正則かどうかを確認
        // 逆行列が存在するかどうかは、行列式が0でないかどうかで判断できる
        // この例では、行列式が0であれば再度行列を生成する
        if (determinant(matrix) != 0) {
            break;
        } else {
            // 行列を再生成
            for (int i = 0; i < MATRIX_SIZE; ++i) {
                for (int j = 0; j < MATRIX_SIZE; ++j) {
                    matrix[i][j] = rand() % 2;
                }
            }
        }
    }
}

int is_reg(MTX cc, MTX *R)
{
  unsigned char inv_a[F][F] = {{0}}; // 逆行列を格納する配列
  unsigned char cl[F][F];            // 行列 cc のコピー
  unsigned char b[F][F] = {{0}};     // 検算用の一時的な行列
  unsigned int flg = 0, count = 0;


  while (flg != F)
  {
    // 行列 cc をコピー
    for (int i = 0; i < F; i++)
    {
      for (int j = 0; j < F; j++)
      {
        cl[i][j] = cc.x[i][j];
      }
    }

    // 単位行列を作成
    for (int i = 0; i < F; i++)
    {
      for (int j = 0; j < F; j++)
      {
        inv_a[i][j] = (i == j) ? 1 : 0;
      }
    }

    // 掃き出し法
    for (int i = 0; i < F; i++)
    {
      if (cc.x[i][i] == 0)
      {
        int j = i + 1;
        while (j < F && cc.x[j][i] == 0)
        {
          j++;
        }
        if (j == F)
        {
          printf("S is not regular.\n");
          return -1;
        }
        for (int k = 0; k < F; k++)
        {
          cc.x[i][k] ^= cc.x[j][k] % 2;
          inv_a[i][k] ^= inv_a[j][k];
        }
        cc.x[i][i] = 1;
      }
      if (cc.x[i][i] == 1)
      {
        for (int l = i + 1; l < F; l++)
        {
          if (cc.x[l][i] == 1)
          {
            for (int k = 0; k < F; k++)
            {
              cc.x[l][k] ^= cc.x[i][k] % 2;
              inv_a[l][k] ^= inv_a[i][k];
            }
          }
        }
      }
    }

    for (int i = 1; i < F; i++)
    {
      for (int k = 0; k < i; k++)
      {
        if (cc.x[k][i] == 1)
        {
          for (int j = 0; j < F; j++)
          {
            cc.x[k][j] ^= cc.x[i][j] % 2;
            inv_a[k][j] ^= inv_a[i][j];
          }
        }
      }
    }

    // 検算
    for (int i = 0; i < F; i++)
    {
      for (int j = 0; j < F; j++)
      {
        for (int k = 0; k < F; k++)
        {
          b[i][j] ^= (cl[i][k] & inv_a[k][j]);
        }
      }
    }

    for (int i = 0; i < F; i++)
    {
      if (b[i][i] == 1)
      {
        flg++;
      }
    }

    count = 0;
    for (int i = 0; i < F; i++)
    {
      for (int j = 0; j < F; j++)
      {
        printf("%d,",b[i][j]);
        if (b[i][j] == 0 && i != j)
        {
          count++;
        }
      }
      printf("\n");
    }
    printf("\n%d,%d %d,%d\n", flg, F, count, F * F - F);

    if (flg == F && count == (F * F - F))
    {
      printf("inv_S[K][K]=\n{\n");
      for (int i = 0; i < F; i++)
      {
        printf("{");
        for (int j = 0; j < F; j++)
        {
          R->x[i][j] = inv_a[i][j];
          printf("%d,", inv_a[i][j]);
        }
        printf("},\n");
      }
      printf("};\n");
      printf("S is regular\n");
      return 0;
    }
    return -1;
  }

  return -1;
}

int mkS(MTX cc, MTX *R)
{
  int i, j, k, l;
  unsigned char b[X][X] = {0};
  unsigned char dd[X] = {0};
  unsigned int flg = 0, count = 0;
  // unsigned char cc.x[X][X] = {0};
  unsigned char cl[X][X];
  time_t t;
  FILE *fq;
  unsigned char inv_a[X][X] = {0}; // ここに逆行列が入る
  unsigned char buf;               // 一時的なデータを蓄える
  int n = K * E;                   // 配列の次数

  // while(flg!=F || count!=F*F-F)
  // while(count!=F*F-F)
  while (flg != X)
  {
  labo:
    // memset(cc,0,sizeof(cc));
    flg = 0;
    count = 0;
    srand(clock() + time(&t));

    // g2();
    /*
    for (i = 0; i < X; i++)
    {
      for (j = 0; j < X; j++)
        cc.x[i][j] = xor128() % 2;
    }
    */
    printf("end of g2\n");
    // exit(1);

    // #pragma omp parallel for private(j)
    for (i = 0; i < X; i++)
    {

      for (j = 0; j < X; j++)
      {
        // printf("%d,",cc.x[i][j]);
        cl[i][j] = cc.x[i][j];
        dd[j] = cc.x[i][j];
      }
      // printf("\n");
    }

    // memset(inv_a,0,sizeof(inv_a));

    // 単位行列を作る
    // #pragma omp parallel for private(j)
    for (i = 0; i < X; i++)
    {
      for (j = 0; j < X; j++)
      {
        inv_a[i][j] = (i == j) ? 1.0 : 0.0;
      }
    }

    // 掃き出し法

    for (i = 0; i < X; i++)
    {
      if (cc.x[i][i] == 0)
      {
        j = i;
        /*
  cc.x[i][i]=1;
  for(k=i+1;k<F;k++)
    cc.x[i][k]^=rand()%2;
  //printf("i=%d\n",i);
  */

        while (cc.x[j][i] == 0 && j < X)
        {
          j++;
          // buf=cc.x[j++][i];
        }

        //  cc.x[i][i]=1;
        //  printf("j=%d\n",j);

        //  exit(1);
        // #pragma omp parallel for
        if (j >= X)
        {
          printf("baka in mkS %d\n", j);
          // exit(1);
          return -1;
          // goto labo;
        }
        for (k = 0; k < X; k++)
        {
          cc.x[i][k] ^= cc.x[j][k] % 2;
          inv_a[i][k] ^= inv_a[j][k];
        }

        cc.x[i][i] = 1;
      }
      //  exit(1);

      if (cc.x[i][i] == 1)
      {
        for (l = i + 1; l < X; l++)
        {
          if (cc.x[l][i] == 1)
          {
            // #pragma omp parallel for
            for (k = 0; k < X; k++)
            {
              cc.x[l][k] ^= cc.x[i][k] % 2;
              inv_a[l][k] ^= inv_a[i][k];
            }
          }
        }

        // printf("@%d\n",i);
      }
      // printf("@i=%d\n",i);
    }

    //  exit(1);
    // #pragma omp parallel for private(j,k)
    for (i = 1; i < X; i++)
    {
      for (k = 0; k < i; k++)
      {
        if (cc.x[k][i] == 1)
        {
          for (j = 0; j < X; j++)
          {
            // if(a[k][i]==1){
            cc.x[k][j] ^= cc.x[i][j] % 2;
            inv_a[k][j] ^= inv_a[i][j];
            //}
          }
        }
      }
    }

    /*
        //逆行列を出力
        for (i = 0; i < F; i++)
        {
          for (j = 0; j < F; j++)
          {
            printf("a %d,", inv_a[i][j]);
          }
          printf("\n");
        }
    */
    // exit(1);

    // 検算
    // #pragma omp parallel for private(j, k) num_threads(16)
    for (i = 0; i < X; i++)
    {
      // #pragma omp parallel num_threads(8) //private(j,k)
      {
        for (j = 0; j < X; j++)
        {
          l = 0;
          // #pragma omp parallel for reduction(^:l)
          for (k = 0; k < X; k++)
          {
            b[i][j] ^= (cl[i][k] & inv_a[k][j]);
            // l^=(cl[i][k]&inv_a[k][j]);
          }
          // b[i][j]=l;
        }
      }
    }

    for (i = 0; i < X; i++)
    {
      //   printf("%d",b[i][i]);
      // printf("==\n");
      if (b[i][i] == 1)
      {
        // printf("baka");
        //    exit(1);
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

    // if(cl[0][0]>0)
    //   goto labo;
    //
    printf("S[K][K]=\n{\n");
    if (flg == X && count == (X * X - X))
    // if(flg==F)
    {
      for (i = 0; i < X; i++)
      {
        // printf("{");
        for (j = 0; j < X; j++)
        {
          //
          dd[j] = cl[i][j];
          S.x[i][j] = cl[i][j];
          printf("%d,", S.x[i][j]);
        }

        printf("},\n");
      }
      printf("};\n");

      printf("inv_S[K][K]=\n{\n");
      for (i = 0; i < X; i++)
      {
        // printf("{");
        for (j = 0; j < X; j++)
        {
          dd[j] = inv_a[i][j];
          R->x[i][j] = inv_a[i][j];
          // printf("%d,", inv_S.w[i][j]);
        }
        // printf("},\n");
      }
      // printf("};\n");

      /*
            for (i = 0; i < F; i++)
            {
              for (j = 0; j < F; j++)
                printf("%d, ", b[i][j]);
              printf("\n");
            }
            //  exit(1);
            */
      return 0;
    }
    return -1;
  }
  return -1;
}

int mkRS(MTX cc, MTX *R)
{
  int i, j, k, l;
  unsigned char b[AY][AY] = {0};
  unsigned char dd[AY] = {0};
  unsigned int flg = 0, count = 0;
  // unsigned char cc.x[X][X] = {0};
  unsigned char cl[AY][AY];
  time_t t;
  FILE *fq;
  unsigned char inv_a[AY][AY] = {0}; // ここに逆行列が入る
  unsigned char buf;                 // 一時的なデータを蓄える
  int n = K * E;                     // 配列の次数

  // while(flg!=F || count!=F*F-F)
  // while(count!=F*F-F)
  while (flg != AY)
  {
  labo:
    // memset(cc,0,sizeof(cc));
    flg = 0;
    count = 0;
    srand(clock() + time(&t));

    // g2();
    /*
    for (i = 0; i < X; i++)
    {
      for (j = 0; j < X; j++)
        cc.x[i][j] = xor128() % 2;
    }
    */
    printf("end of g2\n");
    // exit(1);

    // #pragma omp parallel for private(j)
    for (i = 0; i < AY; i++)
    {

      for (j = 0; j < AY; j++)
      {
        // printf("%d,",cc.x[i][j]);
        cl[i][j] = cc.x[i][j];
        dd[j] = cc.x[i][j];
      }
      // printf("\n");
    }

    // memset(inv_a,0,sizeof(inv_a));

    // 単位行列を作る
    // #pragma omp parallel for private(j)
    for (i = 0; i < AY; i++)
    {
      for (j = 0; j < AY; j++)
      {
        inv_a[i][j] = (i == j) ? 1.0 : 0.0;
      }
    }

    // 掃き出し法

    for (i = 0; i < AY; i++)
    {
      if (cc.x[i][i] == 0)
      {
        j = i;
        /*
  cc.x[i][i]=1;
  for(k=i+1;k<F;k++)
    cc.x[i][k]^=rand()%2;
  //printf("i=%d\n",i);
  */

        while (cc.x[j][i] == 0 && j < AY)
        {
          j++;
          // buf=cc.x[j++][i];
        }

        //  cc.x[i][i]=1;
        //  printf("j=%d\n",j);

        //  exit(1);
        // #pragma omp parallel for
        if (j >= AY)
        {
          printf("baka in mkS %d\n", j);
          // exit(1);
          return -1;
          // goto labo;
        }
        for (k = 0; k < AY; k++)
        {
          cc.x[i][k] ^= cc.x[j][k] % 2;
          inv_a[i][k] ^= inv_a[j][k];
        }

        cc.x[i][i] = 1;
      }
      //  exit(1);

      if (cc.x[i][i] == 1)
      {
        for (l = i + 1; l < AY; l++)
        {
          if (cc.x[l][i] == 1)
          {
            // #pragma omp parallel for
            for (k = 0; k < AY; k++)
            {
              cc.x[l][k] ^= cc.x[i][k] % 2;
              inv_a[l][k] ^= inv_a[i][k];
            }
          }
        }

        // printf("@%d\n",i);
      }
      // printf("@i=%d\n",i);
    }

    //  exit(1);
    // #pragma omp parallel for private(j,k)
    for (i = 1; i < AY; i++)
    {
      for (k = 0; k < i; k++)
      {
        if (cc.x[k][i] == 1)
        {
          for (j = 0; j < AY; j++)
          {
            // if(a[k][i]==1){
            cc.x[k][j] ^= cc.x[i][j] % 2;
            inv_a[k][j] ^= inv_a[i][j];
            //}
          }
        }
      }
    }

    /*
        //逆行列を出力
        for (i = 0; i < F; i++)
        {
          for (j = 0; j < F; j++)
          {
            printf("a %d,", inv_a[i][j]);
          }
          printf("\n");
        }
    */
    // exit(1);

    // 検算
    // #pragma omp parallel for private(j, k) num_threads(16)
    for (i = 0; i < AY; i++)
    {
      // #pragma omp parallel num_threads(8) //private(j,k)
      {
        for (j = 0; j < AY; j++)
        {
          l = 0;
          // #pragma omp parallel for reduction(^:l)
          for (k = 0; k < AY; k++)
          {
            b[i][j] ^= (cl[i][k] & inv_a[k][j]);
            // l^=(cl[i][k]&inv_a[k][j]);
          }
          // b[i][j]=l;
        }
      }
    }

    for (i = 0; i < AY; i++)
    {
      //   printf("%d",b[i][i]);
      // printf("==\n");
      if (b[i][i] == 1)
      {
        // printf("baka");
        //    exit(1);
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

    // if(cl[0][0]>0)
    //   goto labo;
    //
    printf("S[K][K]=\n{\n");
    if (flg == AY && count == (AY * AY - AY))
    // if(flg==F)
    {
      for (i = 0; i < AY; i++)
      {
        // printf("{");
        for (j = 0; j < AY; j++)
        {
          //
          dd[j] = cl[i][j];
          S.x[i][j] = cl[i][j];
          printf("%d,", S.x[i][j]);
        }

        printf("},\n");
      }
      printf("};\n");

      printf("inv_S[K][K]=\n{\n");
      for (i = 0; i < AY; i++)
      {
        // printf("{");
        for (j = 0; j < AY; j++)
        {
          dd[j] = inv_a[i][j];
          R->x[i][j] = inv_a[i][j];
          // printf("%d,", inv_S.w[i][j]);
        }
        // printf("},\n");
      }
      // printf("};\n");

      /*
            for (i = 0; i < F; i++)
            {
              for (j = 0; j < F; j++)
                printf("%d, ", b[i][j]);
              printf("\n");
            }
            //  exit(1);
            */
      return 0;
    }
    return -1;
  }
  return -1;
}

int binv(MTX cc, MTX *L, int Y)
{
  int i, j, k, l;
  unsigned char b[N][N] = {0};
  unsigned char dd[N] = {0};
  unsigned int flg = 0, count = 0;
  // unsigned char cc.x[X][X] = {0};
  unsigned char cl[N][N];
  time_t t;
  FILE *fq;
  unsigned char inv_a[N][N] = {0}; // ここに逆行列が入る
  unsigned char buf;               // 一時的なデータを蓄える
  int n = K * E;                   // 配列の次数

  // while(flg!=F || count!=F*F-F)
  // while(count!=F*F-F)
  while (flg != Y)
  {
  labo:
    // memset(cc,0,sizeof(cc));
    flg = 0;
    count = 0;
    srand(clock() + time(&t));

    // g2();
    /*
    for (i = 0; i < Y; i++)
    {
      for (j = 0; j < Y; j++)
        cc.x[i][j] = xor128() % 2;
    }
    */
    printf("end of g2\n");
    // exit(1);

#pragma omp parallel for private(j)
    for (i = 0; i < Y; i++)
    {

      for (j = 0; j < Y; j++)
      {
        // printf("%d,",cc.x[i][j]);
        cl[i][j] = cc.x[i][j];
        dd[j] = cc.x[i][j];
      }
      // printf("\n");
    }

// memset(inv_a,0,sizeof(inv_a));

// 単位行列を作る
#pragma omp parallel for private(j)
    for (i = 0; i < Y; i++)
    {
      for (j = 0; j < Y; j++)
      {
        inv_a[i][j] = (i == j) ? 1.0 : 0.0;
      }
    }

    // 掃き出し法

    for (i = 0; i < Y; i++)
    {
      if (cc.x[i][i] == 0)
      {
        j = i;
        /*
  cc.x[i][i]=1;
  for(k=i+1;k<F;k++)
    cc.x[i][k]^=rand()%2;
  //printf("i=%d\n",i);
  */

        while (cc.x[j][i] == 0 && j < Y)
        {
          j++;
          // buf=cc.x[j++][i];
        }

        //  cc.x[i][i]=1;
        //  printf("j=%d\n",j);

        //  exit(1);
        // #pragma omp parallel for
        if (j >= Y)
        {
          printf("baka in binv %d\n", j);
          // exit(1);
          return -1;
          // goto labo;
        }
        for (k = 0; k < Y; k++)
        {
          cc.x[i][k] ^= cc.x[j][k] % 2;
          inv_a[i][k] ^= inv_a[j][k];
        }

        cc.x[i][i] = 1;
      }
      //  exit(1);

      if (cc.x[i][i] == 1)
      {
        for (l = i + 1; l < Y; l++)
        {
          if (cc.x[l][i] == 1)
          {
            // #pragma omp parallel for
            for (k = 0; k < Y; k++)
            {
              cc.x[l][k] ^= cc.x[i][k] % 2;
              inv_a[l][k] ^= inv_a[i][k];
            }
          }
        }

        // printf("@%d\n",i);
      }
      // printf("@i=%d\n",i);
    }

    //  exit(1);
    // #pragma omp parallel for private(j,k)
    for (i = 1; i < Y; i++)
    {
      for (k = 0; k < i; k++)
      {
        if (cc.x[k][i] == 1)
        {
          for (j = 0; j < Y; j++)
          {
            // if(a[k][i]==1){
            cc.x[k][j] ^= cc.x[i][j] % 2;
            inv_a[k][j] ^= inv_a[i][j];
            //}
          }
        }
      }
    }

    /*
        //逆行列を出力
        for (i = 0; i < F; i++)
        {
          for (j = 0; j < F; j++)
          {
            printf("a %d,", inv_a[i][j]);
          }
          printf("\n");
        }
    */
    // exit(1);

    // 検算
#pragma omp parallel for private(j, k) num_threads(16)
    for (i = 0; i < Y; i++)
    {
      // #pragma omp parallel num_threads(8) //private(j,k)
      {
        for (j = 0; j < Y; j++)
        {
          l = 0;
          // #pragma omp parallel for reduction(^:l)
          for (k = 0; k < Y; k++)
          {
            b[i][j] ^= (cl[i][k] & inv_a[k][j]);
            // l^=(cl[i][k]&inv_a[k][j]);
          }
          // b[i][j]=l;
        }
      }
    }

    for (i = 0; i < Y; i++)
    {
      //   printf("%d",b[i][i]);
      // printf("==\n");
      if (b[i][i] == 1)
      {
        // printf("baka");
        //    exit(1);
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

    // if(cl[0][0]>0)
    //   goto labo;
    //
    printf("S[K][K]=\n{\n");
    if (flg == Y && count == (Y * Y - Y))
    // if(flg==F)
    {
      for (i = 0; i < Y; i++)
      {
        // printf("{");
        for (j = 0; j < Y; j++)
        {
          //
          dd[j] = cl[i][j];
          S.x[i][j] = cl[i][j];
          printf("%d,", S.x[i][j]);
        }

        printf("},\n");
      }
      printf("};\n");

      printf("inv_S[K][K]=\n{\n");
      for (i = 0; i < Y; i++)
      {
        // printf("{");
        for (j = 0; j < Y; j++)
        {
          dd[j] = inv_a[i][j];
          L->x[i][j] = inv_a[i][j];
          // printf("%d,", inv_S.w[i][j]);
        }
        // printf("},\n");
      }
      // printf("};\n");

      /*
            for (i = 0; i < F; i++)
            {
              for (j = 0; j < F; j++)
                printf("%d, ", b[i][j]);
              printf("\n");
            }
            //  exit(1);
            */
      return 0;
    }
    return -1;
  }

  /*
  fq=fopen("S.key","wb");
  for(i=0;i<F;i++){
    for(j=0;j<F;j++)
      dd[j]=cl[i][j];
    fwrite(dd,1,n,fq);

  }
  fclose(fq);
  fq=fopen("inv_S.key","wb");
  for(i=0;i<F;i++){
    for(j=0;j<F;j++)
      dd[j]=inv_a[i][j];
    fwrite(dd,1,n,fq);
  }
  fclose(fq);
*/

  // free(b);
  return -1;
}
