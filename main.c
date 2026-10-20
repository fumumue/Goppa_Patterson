#include "common.h"

void mkerr(unsigned short *z1, int num);
int ero2(vec v);
unsigned short oinv(unsigned short a);
unsigned short oinv2(unsigned short a);

int main(void)
{
  unsigned short z1[N] = {0}; //{1,0,1,1,1,0,0,0,0,0,1,1,1,0,0,1};
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

  w = pubkeygen();

  while (1)
  {
    memset(z1, 0, sizeof(z1));
    mkerr(z1, K);

    f = zin(z1);
    printpol(f);
    printf(" ==syndrome\n");

    v = patterson(w, f);

    if (ero2(v) < 0)
    {
      for (int i = 0; i < N; i++)
      {
        if (z1[i] > 0)
          printf("i=%d\n", i);
      }
    }

    break;
  }

  return 0;
}
