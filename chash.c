#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>


//#define O 8  //素数の和(3*5=15 の位数を持つ)
#define MAX 2 // 素数表の先頭から何個素数を足すか
#define NN 8  //置換配列の次元

unsigned char p[2] = {3, 5}; // 8
// unsigned char p[2]={2,3}; //5

int x[5][NN] = {0};

unsigned char aa[64]; //={ 148, 246, 52, 251, 16, 194, 72, 150, 249, 23, 90, 107, 151, 42, 154, 124, 48, 58, 30, 24, 42, 33, 38, 10, 115, 41, 164, 16, 33, 32, 252, 143, 86, 175, 8, 132, 103, 231, 95, 190, 61, 29, 215, 75, 251, 248, 72, 48, 224, 200, 147, 93, 112, 25, 227, 223, 206, 137, 51, 88, 109, 214, 17, 172};

#define str_length 128
#define password_length 256

#define SIZE_OF_ARRAY(array) (sizeof(array) / sizeof(array[0]))
#define SWAP(type, a, b) \
  {                      \
    type work = a;       \
    a = b;               \
    b = work;            \
  }

#define ROTL32(X, B) rotl32((X), (B))
static inline uint32_t
rotl32(const uint32_t x, const int b)
{
  return (x << b) | (x >> (32 - b));
}

#define ROTL64(X, B) rotl64((X), (B))
static inline uint64_t
rotl64(const uint64_t x, const int b)
{
  return (x << b) | (x >> (64 - b));
}

char password[password_length + 1];

int desc(const void *a, const void *b)
{
  return *(int *)b - *(int *)a;
}

int compareInt(const void *a, const void *b)
{
  int aNum = *(int *)a;
  int bNum = *(int *)b;

  return aNum - bNum;
}

void printArray(const int *array, size_t size)
{
  for (size_t i = 0; i < size; ++i)
  {
    printf("%d ", array[i]);
  }
  printf("\n");
}

void mkcycle()
{
  int i, j, r, rr, vv, pko, flg, flg2, l, n, ll, jj, flg3, cnt2, v, cnt3 = 0, ii;
  unsigned short cnt, k, count = 0, kk, count2 = 0;
  ;
  unsigned long long int o;

  // p は素数表
  //{2,13,17};
  n = 0;
  for (i = 0; i < MAX; i++)
  {
    n += p[i];
  }

  // printf("%d\n",n);
  //    exit(1);

  // scanf("%llu",&o);
  // srand(o);

a3:
  for (j = 0; j < 3; j++)
  {

    cnt = 0;

    count2 = 0;
    for (i = 0; i < NN; i++)
      x[j][i] = -1;

    // printf("aa");
    cnt = 0;
    pko = 0;

    for (i = 0; i < MAX; i++)
    {
      count = 0;
      flg = 0;

      while (count < p[i])
      {

        // printf("cx=%d %d\n",count,i);
        if (i == NN)
        {
          cnt3 = 0;
          for (ii = 0; ii < NN; ii++)
          {
            if (x[j][ii] == -1)
              cnt3++;
          }
          // printf("-1,%d\n%d",cnt3);
          if (cnt3 == 0)
            break;
        }

        if (count == 0)
        {
          do
          {
            // printf("-----------------------\n");
            kk = rand() % NN;
          } while (x[j][kk] >= 0);

          pko = kk;
          // printf("pko=%d\n",kk);
          flg = 1;
        }

        if (x[j][kk] == -1)
        {

          do
          {
            // printf("aa");
            flg2 = 0;
            do
            {
              // printf("bb");
              k = rand() % NN;
            } while (k == pko);
            for (l = 0; l < NN; l++)
            {
              if (x[j][l] == k)
                flg2 = 1;
            }
            if (flg2 == 0 && kk != k && x[j][kk] == -1)
              x[j][kk] = k;

          } while (kk == k || flg2 == 1 || k == pko);
        }
        else
        {
          do
          {
            kk = rand() % NN;
            // printf("cc");
          } while (kk == pko);

          do
          {
            flg2 = 0;

            k = rand() % NN;
            for (l = 0; l < NN; l++)
            {
              if (x[j][l] == k)
                flg2 = 1;
            }
            if (flg2 == 0 && kk != k && x[j][kk] == -1 && k != pko)
            {
              x[j][kk] = k;
              count++;
            }
            if (count == p[i] - 1)
              break;
            // printf("ee");
          } while (flg2 == 1 || kk == k || k == pko);
        }

        if (flg2 == 0 && kk != k)
        {
          flg2 = 0;
          if (x[j][kk] == -1)
          {
            for (l = 0; l < NN; l++)
            {
              if (x[j][l] == k)
                flg2 = 1;
            }
            if (flg2 == 0)
            {
              x[j][kk] = k;
              count++;
            }
          }
          while (x[j][kk] > -1)
            kk = rand() % NN;
          if (k != pko)
            kk = k;
          count++;
          cnt++;
        }
        else if (flg2 == 1 || kk == k)
        {
          kk = rand() % NN;

          do
          {
            flg2 = 0;
            // printf("ff");
            k = rand() % NN;
            for (ll = 0; ll < NN; ll++)
            {
              if (x[j][ll] == k)
                flg2 = 1;
            }
            while (flg2 == 1 || k == pko || kk == k)
              ;

            if (flg2 == 0 && k != pko && kk != k)
            {
              x[j][kk] = k;
              count++;
              kk = k;
            }

          } while (flg2 == 1);
        }
        // printf("loop=%d,%d\n",k,kk);

        if (count == p[i] - 1 && pko != kk)
        {
          flg2 = 0;
          for (l = 0; l < NN; l++)
          {
            if (x[j][l] == pko)
              flg2 = 1;
          }
          if (flg2 == 0)
          {
            x[j][kk] = pko;
            count++;
            cnt++;
          }
        }

        // printf("count=%d\n",cnt);
        // printf("i=%d\n",i);
        if (count - 1 == p[i])
        {
          x[j][k] = pko;
          cnt++;
        }
        // printf("%d\n",cnt3);
      }
    }

    // printf("cnt=%d\n",cnt);

    for (i = 0; i < NN; i++)
    {
      if (x[j][i] == -1)
        goto a3;
    }
  }
}

int chkp(int *x)
{
  int i, j, k, count = 1, cnt[NN] = {0};
  unsigned char c[NN] = {0};
  int cnt2 = 0;

  /*
  for(i=0;i<N;i++)
  printf("%d,",x[i]);
  printf("\n");
  */
  i = 0;
  while (cnt2 < N)
  {
    k = x[i];
    i = k;
    // printf("k=%d\n",k);
    while (1)
    {
      j = x[i];

      // printf("j=%d\n",j);
      i = j;
      // if(j==k)
      // break;

      if (c[j] == 0)
      {
        c[j] = count;
        cnt[count]++;
        cnt2++;
      }
      if (i == k)
      {
        i = 0;
        // if(cnt2!=N)
        count++;
        while (c[i] != 0)
          i++;
        // if(count==N)
        break;
      }
    }
    // break;
  }

  // for(i=0;i<count;i++)
  // cnt2+=cnt[i];
  /*
  printf("cycle1=%d\n",cnt2);
  for(i=0;i<N;i++)
  printf("%d,",c[i]);
  printf("\n");
  */
  qsort(cnt, SIZE_OF_ARRAY(cnt), sizeof(int), compareInt);
  // printArray(cnt, N);
  i = 0;
  while (cnt[i] == 0)
    i++;
  int i2 = 0;
  for (j = 0; j < MAX; j++)
  {
    if (cnt[i + j] == p[j])
      i2++;
  }
  if (cnt2 == N && i2 == MAX)
  {
    // exit(1);
    // printf("@\n");
    return 1;
  }
  return 0;
}

/*
    Fisher-Yates shuffle による方法
    配列の要素をランダムシャッフルする
*/
void random_shuffle(unsigned short *array, size_t size)
{
  for (size_t i = size; i > 1; --i)
  {
    size_t a = i - 1;
    size_t b = rand() % i;
    SWAP(int, array[a], array[b]);
  }
}

/*
    配列の要素を出力
*/
void print_array(const unsigned short *array, size_t size)
{
  for (size_t i = 0; i < size; ++i)
  {
    printf("%d ", array[i]);
  }
  printf("\n");
}

unsigned long xor128(void)
{
  unsigned int a = 0;

  static unsigned long x = 123456789, y = 362436069, z = 521288629, w = 88675123;
  unsigned long t;

  a = rand();
  t = x ^ (a << 11);
  a = y;
  y = z;
  z = w;
  return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
}

void seed(void)
{
  /*
   * 変数宣言
   */
  char str[str_length + 1];
  time_t t;
  int i, j, k, rnd;

  /*
   * 乱数の初期化
   */
  srand(clock() + time(&t));

  /*
   * 乱数生成とパスワードの生成
   */
  for (i = 0; i < str_length; i++)
  {
    for (j = 0; j < 2; j++)
    {
      k = i * 2 + j;
      do
      {
        rnd = rand();
        password[k] = str[i] + rnd;
      } while (!isalnum(password[k]));
    }
  }

  /*
   * NULL文字の挿入
   */
  password[password_length] = '\0';

  /*
   * パスワードの出力
   */
  //    printf("生成パスワード：%s",password);

  return;
}

int mlt(int x, int y)
{

  if (x == 0 || y == 0)
    return 0;

  return ((x + y - 2) % (M - 1)) + 1;
}

int mlt2(int n, int x)
{
  int i, j;

  if (n == 0)
    return 1;
  i = x;
  for (j = 0; j < n - 1; j++)
    i = mlt(i, x);

  return i;
}

unsigned short tabloog[N] = {0};
unsigned short tabloof[N] = {0};

int mktbl()
{
  int i, j;

  j = 2;
  tabloog[0] = 0;
  tabloog[1] = 1;
  tabloof[0] = 0;
  tabloof[1] = 1;
  for (i = 1; i < N; i = i * 2)
  {
    tabloog[j] = gf[i];
    tabloof[j++] = fg[i];
  }
}


int mltn(int n, int x) {
    int ret = 1;
    while (n > 0) {
        if (n & 1) ret = mlt(ret , x) ;  // n の最下位bitが 1 ならば x^(2^i) をかける
        x = mlt(x , x);
        n >>= 1;  // n を1bit 左にずらす
    }
    return ret;
}



/*
 * S-box transformation table
 */
static const unsigned char s_box[256] = {
    // 0     1     2     3     4     5     6     7     8     9     a     b     c     d     e     f
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,  // 0
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,  // 1
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,  // 2
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,  // 3
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,  // 4
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,  // 5
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,  // 6
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,  // 7
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,  // 8
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,  // 9
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,  // a
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,  // b
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,  // c
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,  // d
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,  // e
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16}; // f

/*
 * Inverse S-box transformation table
 */
static const unsigned char inv_s_box[256] = {
    // 0     1     2     3     4     5     6     7     8     9     a     b     c     d     e     f
    0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb,  // 0
    0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87, 0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb,  // 1
    0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,  // 2
    0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25,  // 3
    0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92,  // 4
    0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,  // 5
    0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06,  // 6
    0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02, 0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b,  // 7
    0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,  // 8
    0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e,  // 9
    0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89, 0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b,  // a
    0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,  // b
    0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f,  // c
    0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d, 0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef,  // d
    0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,  // e
    0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d}; // f

void rp2(unsigned char *a)
{
  int i, j, x;
  time_t t;

  //  srand(clock() + time(&t));

  for (i = 0; i < K; i++)
  {
    a[i] = i;
  }
  for (i = 0; i < K - 2; i++)
  {
    // rand from i+1 to F-1
    j = (rand() % (K - 1 - i)) + i + 1;

    // swap a[i] and a[j]
    x = a[j];
    a[j] = a[i];
    a[i] = x;
  }
  if (a[K - 1] == K - 1)
  {
    a[K - 1] = a[K - 2];
    a[K - 2] = K - 1;
  }
}

int print_uint128(__uint128_t n)
{
  if (n == 0)
    return printf("0\n");

  char str[40] = {0};              // log10(1 << 128) + '\0'
  char *s = str + sizeof(str) - 1; // start at the end
  while (n != 0)
  {
    if (s == str)
      return -1; // never happens

    *--s = "0123456789"[n % 10]; // save last digit
    n /= 10;                     // drop it
  }
  return printf("%s", s);
}

//ハッシュ関数本体
arrayul
chash()
{
  int i; //, j = 0;
  arrayul n;
  arrayul vw;
  unsigned char key[K] = {0}; //{148, 246, 52, 251, 16, 194, 72, 150, 249, 23, 90, 107, 151, 42, 154, 124}; //, 48, 58, 30, 24, 42, 33, 38, 10, 115, 41, 164, 16, 33, 32, 252, 143, 86, 175, 8, 132, 103, 231, 95, 190, 61, 29, 215, 75, 251, 248, 72, 48, 224, 200, 147, 93, 112, 25, 227, 223, 206, 137, 51, 88, 109, 214, 17, 172};
  unsigned char z[K];
  unsigned char x0[K] = {0};
  unsigned char inv_x[K] = {0};
  unsigned char x1[K] = {0};
  unsigned char tmp[K] = {0};
  unsigned short u = 0;
  time_t t;

  rp2(x0);
  rp2(x1);

  n.x[0] = 65535; //(unsigned int)time(&t);
  n.x[1] = 65535; //(unsigned long)clock();
  for (i = 0; i < K; i++)
    key[i] = n.d[i];

  for (i = 0; i < K; i++)
    inv_x[x0[i]] = i;
  // key[i]=rand()%2;//salt[i];
  // key[0]=1;

  int count = 0;

  // memset(f, 0, sizeof(f));
  while (count < 8)
  {
    //バッファを埋める回数だけ回す
    for (i = 0; i < K; i++)
      z[i] = x0[x1[inv_x[i]]];

    for (i = 0; i < K; i++)
      tmp[i] += s_box[key[z[i]]];
    for (i = 0; i < K; i++)
    {
      key[i] ^= inv_s_box[ROTL8(tmp[i], 3)];
      // printf("%d,",key[i]);
    }
    // printf("\n");
    memcpy(x1, z, sizeof(x1));

    count++;
  }

  int uu = 0;
  for (i = 0; i < K; i++)
  {
    // vw.x[0] ^= key[i];
    // vw.x[0]=(vw.x[0]<<1);
    vw.d[i] = key[i];
  }
  // vw.x[0]=(vw.x[0]>>1);

  return vw;
}
