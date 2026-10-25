#include "common.h"

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

  for (int i = 0; i < N; i++)
    printf("%d %d\n", oinv(i), oinv2(i));

  w = pubkeygen(permutation, &inverse_scrambler, &rng, syndrome_matrix);

  memset(z1, 0, (size_t)N * sizeof(*z1));
  if (mkerr(z1, K, &rng) != 0)
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


  v=paloma(v2o(f),v2o(w));
  //v = patterson(w, f);

  int decoded_count = 0, expected_count = 0, matching_count = 0;
  int decode_matches = verify_decode_result(v, permutation, z1,
                                            &decoded_count, &expected_count,
                                            &matching_count);
  int reported_count = ero2(v, permutation);
  decode_matches = decode_matches && reported_count == decoded_count;
  printf("Patterson decode verification: %s "
         "(matched %d/%d error positions; decoded %d, expected %d)\n",
         decode_matches ? "PASS" : "FAIL",
         matching_count, expected_count, decoded_count, expected_count);

  if (reported_count < 0 || !decode_matches)
  {
    for (int i = 0; i < N; i++)
    {
      if (z1[i] > 0)
        printf("i=%d\n", i);
    }
  }

  matrix_free(&inverse_scrambler);
  free(z1);
  free(permutation);
  free(syndrome_matrix);
  return decode_matches ? EXIT_SUCCESS : EXIT_FAILURE;
}
