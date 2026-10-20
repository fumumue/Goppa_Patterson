#include "common.h"

unsigned short oinv(unsigned short a);
unsigned short oinv2(unsigned short a);

int main(void)
{
  unsigned short z1[N] = {0}; //{1,0,1,1,1,0,0,0,0,0,1,1,1,0,0,1};
  unsigned short permutation[N] = {0};
  unsigned short syndrome_matrix[N][K] = {0};
  MTX inverse_scrambler = {0};
  RNG rng = {123456789, 362436069, 521288629, 88675123};
  vec f = {0}, w = {0}, v = {0};

  srand(clock());
  if (K > N)
  {
    printf("configuration error! K is too big K\n");
    exit(1);
  }
  gen_gf(E, N, 1);

  for (int i = 0; i < N; i++)
    printf("%d %d\n", oinv(i), oinv2(i));

  w = pubkeygen(permutation, &inverse_scrambler, &rng, syndrome_matrix);

  memset(z1, 0, sizeof(z1));
  if (mkerr(z1, K, &rng) != 0)
    return EXIT_FAILURE;

  f = zin(&inverse_scrambler, syndrome_matrix, z1);
  printpol(f);
  printf(" ==syndrome\n");

  v = patterson(w, f);

  if (ero2(v, permutation) < 0)
  {
    for (int i = 0; i < N; i++)
    {
      if (z1[i] > 0)
        printf("i=%d\n", i);
    }
  }

  return 0;
}
