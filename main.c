#include "common.h"
#include <sys/random.h>


unsigned short oinv(unsigned short a);
unsigned short oinv2(unsigned short a);

static int verify_decode_result(vec decoded,
                                const unsigned short permutation[N],
                                const unsigned short expected[N],
                                int *decoded_count,
                                int *expected_count,
                                int *matching_count)
{
  unsigned char decoded_positions[N] = {0};
  int valid = 1;

  *decoded_count = 0;
  *expected_count = 0;
  *matching_count = 0;

  for (int i = 0; i < T * 2; i++)
  {
    unsigned short position = decoded.x[i];
    if (i > 0 && position == 0)
      break;
    if (position >= N || decoded_positions[position])
    {
      valid = 0;
      continue;
    }
    decoded_positions[position] = 1;
    (*decoded_count)++;
  }

  for (int i = 0; i < N; i++)
  {
    int was_expected = expected[i] != 0;
    int was_decoded = decoded_positions[permutation[i]] != 0;
    *expected_count += was_expected;
    *matching_count += was_expected && was_decoded;
    if (was_expected != was_decoded)
      valid = 0;
  }

  return valid;
}

// 復号したエラー d と、正解の z1 を比べる。誤り位置が完全に一致なら PASS
static int check_decode(const unsigned short *d, const unsigned short *z1, int n)
{
    int ok = 1;
    for (int i = 0; i < n; i++)
    {
        if ((d[i] != 0) != (z1[i] != 0)) { ok = 0; break; }
    }

    printf("decode check: %s\n", ok ? "PASS" : "FAIL");

    printf("decoded error positions:");
    for (int i = 0; i < n; i++)
        if (d[i] != 0) printf(" %d", i);
    printf("\n");

    return ok;
}


// e: 本当の誤り (長さ n), d: 復号結果 (長さ n)。誤り位置が非0
// 戻り値: 本当の誤りのうち、正しく見つかった個数
static int count_correct_old(const unsigned short *e, const unsigned short *d,
                         int n, int *T_out, int *extra)
{
    int t = 0, hit = 0, ex = 0;

    for (int i = 0; i < n; i++)
    {
        int ei = (e[i] != 0);
        int di = (d[i] != 0);

        if (ei) t++;
        if (ei && di) hit++;
        if (!ei && di) ex++;
    }

    if (T_out) *T_out = t;
    if (extra) *extra = ex;
    return hit;
}


int main(void)
{
  unsigned short *z1 = calloc(N, sizeof(*z1));
  unsigned short *permutation = malloc((size_t)N * sizeof(*permutation));
  unsigned short (*syndrome_matrix)[K] = calloc(N, sizeof(*syndrome_matrix));
  MTX inverse_scrambler = {0};
  RNG rng = {123456789, 362436069, 521288629, 88675123};
  vec f = {0}, w = {0}, v = {0};

  srand(clock());
  
  if (z1 == NULL || permutation == NULL || syndrome_matrix == NULL)
  {
    fprintf(stderr, "main: failed to allocate input/output buffers\n");
    free(z1);
    free(permutation);
    free(syndrome_matrix);
    return EXIT_FAILURE;
  }

  if (K > N)
  {
    printf("configuration error! K is too big K\n");
    exit(1);
  }
  gen_gf(E, N, 1);

  //for (int i = 0; i < N; i++)
  //  printf("%d %d\n", oinv(i), oinv2(i));

  w = pubkeygen(permutation, &inverse_scrambler, &rng, syndrome_matrix);

  memset(z1, 0, (size_t)N * sizeof(*z1));
  //if (mkerr(z1, rand()%K+1, &rng) != 0)
  if(mkerr(z1,K,&rng) != 0)
  {
    free(z1);
    free(permutation);
    free(syndrome_matrix);
    matrix_free(&inverse_scrambler);
    return EXIT_FAILURE;
  }

  f = zin(&inverse_scrambler, syndrome_matrix, z1);
  printpol(f);
  printf(" ==syndrome\n");


  unsigned short ya[N] = {0};

//v = paloma(v2o(f), v2o(w));
 v = patterson(w, f);

int reported = ero2(v, permutation, ya);   // ya に元の並びの誤りを書き込む

int expected = 0, matched = 0, extra = 0;
for (int i = 0; i < N; i++)
{
    int e = (z1[i] != 0);
    int d = (ya[i] != 0);
    expected += e;
    matched  += (e && d);
    extra    += (!e && d);
}


int ok = (matched == expected && extra == 0 && reported == expected),count=0;

if (ok)
{
    printf("PASS\n");
    for(int i = 0; i < N; i++)
    {
        if (z1[i] != 0 && ya[i] != 0){
            printf("i= %d %d\n", i,count++);
        }
    }
}
else
{
    printf("FAIL matched %d/%d, extra %d, reported %d\n",
           matched, expected, extra, reported);
    for (int i = 0; i < N; i++)
    {
        if ((ya[i] != 0) != (z1[i] != 0))
            printf("%s i=%d\n", z1[i] ? "missed" : "extra ", i);
    }
}
  matrix_free(&inverse_scrambler);
  free(z1);
  free(permutation);
  free(syndrome_matrix);


return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
