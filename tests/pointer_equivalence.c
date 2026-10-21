#include "../common.h"

static int expect_vec(const char *name, vec actual, const unsigned short expected[DEG])
{
  int i;

  for (i = 0; i < DEG; i++)
  {
    if (actual.x[i] != expected[i])
    {
      fprintf(stderr, "%s[%d]: expected %u, got %u\n",
              name, i, expected[i], actual.x[i]);
      return 1;
    }
  }

  return 0;
}

int main(void)
{
  vec a = {0}, b = {0}, expected = {0}, expected_product = {0};
  vec expected_vector_product = {0};
  unsigned short xor16[67], xor16_expected[67], xor16_mod2[67], xor16_mod2_expected[67];
  uint8_t xor8[67], xor8_expected[67];
  MTX full_rank_matrix = {0};
  MTX singular_matrix = {0}, inverse_matrix = {0};
  unsigned short rank_rows[DEG][DEG] = {0};
  RNG rng = {123456789, 362436069, 521288629, 88675123};
  int failed = 0;

  gen_gf(E, N, 1);
  if (random_full_rank_binary_matrix(&full_rank_matrix, &rng) != 0)
  {
    fprintf(stderr, "full-rank matrix generation failed\n");
    return 1;
  }

  for (int i = 0; i < DEG; i++)
  {
    for (int j = 0; j < DEG; j++)
    {
      if (full_rank_matrix.x[i][j] > 1)
      {
        fprintf(stderr, "generated matrix contains a non-binary value\n");
        failed = 1;
      }
      rank_rows[i][j] = full_rank_matrix.x[i][j];
    }
  }

  int rank = 0;
  for (int col = 0; col < DEG; col++)
  {
    int pivot = rank;
    while (pivot < DEG && rank_rows[pivot][col] == 0)
      pivot++;
    if (pivot == DEG)
      continue;

    for (int j = 0; j < DEG; j++)
    {
      unsigned short tmp = rank_rows[rank][j];
      rank_rows[rank][j] = rank_rows[pivot][j];
      rank_rows[pivot][j] = tmp;
    }
    for (int row = rank + 1; row < DEG; row++)
    {
      if (rank_rows[row][col] == 1)
      {
        for (int j = col; j < DEG; j++)
          rank_rows[row][j] ^= rank_rows[rank][j];
      }
    }
    rank++;
  }
  if (rank != DEG)
  {
    fprintf(stderr, "generated binary matrix rank: expected %d, got %d\n", DEG, rank);
    failed = 1;
  }
  if (is_reg(&full_rank_matrix, &inverse_matrix, &rng) != 0)
  {
    fprintf(stderr, "is_reg failed to invert a full-rank matrix\n");
    failed = 1;
  }
  else
  {
    for (int row = 0; row < DEG; row++)
    {
      for (int col = 0; col < DEG; col++)
      {
        unsigned short product = 0;
        for (int k = 0; k < DEG; k++)
          product ^= full_rank_matrix.x[row][k] & inverse_matrix.x[k][col];
        if (product != (unsigned short)(row == col))
        {
          fprintf(stderr, "matrix inverse mismatch at [%d][%d]\n", row, col);
          failed = 1;
        }
      }
    }
  }
  if (is_reg(&singular_matrix, &inverse_matrix, &rng) == 0)
  {
    fprintf(stderr, "is_reg accepted a singular matrix\n");
    failed = 1;
  }

  unsigned short invalid_error[N] = {0};
  if (mkerr(invalid_error, N + 1, &rng) == 0)
  {
    fprintf(stderr, "mkerr accepted an impossible error weight\n");
    failed = 1;
  }
  invalid_error[0] = 1;
  if (mkerr(invalid_error, N, &rng) == 0)
  {
    fprintf(stderr, "mkerr accepted an impossible additional error weight\n");
    failed = 1;
  }

  unsigned short generated_errors[N] = {0};
  RNG zero_rng = {0};
  if (mkerr(generated_errors, K, &zero_rng) != 0)
  {
    fprintf(stderr, "mkerr failed with a zero-initialized RNG state\n");
    failed = 1;
  }
  else
  {
    int error_weight = 0;
    for (int i = 0; i < N; i++)
      error_weight += generated_errors[i] != 0;
    if (error_weight != K)
    {
      fprintf(stderr, "generated error weight: expected %d, got %d\n", K, error_weight);
      failed = 1;
    }
  }

  for (int i = 0; i < 67; i++)
  {
    xor16[i] = xor16_expected[i] = (unsigned short)(i * 13);
    xor16_mod2[i] = xor16_mod2_expected[i] = (unsigned short)(i * 7);
    xor8[i] = xor8_expected[i] = (uint8_t)(i * 11);
  }

  for (int i = 0; i < 67; i++)
  {
    xor16_expected[i] ^= (unsigned short)(i * 5);
    xor16_mod2_expected[i] ^= (unsigned short)((i * 9) % 2);
    xor8_expected[i] ^= (uint8_t)(i * 3);
  }

  {
    unsigned short xor16_src[67], xor16_mod2_src[67];
    uint8_t xor8_src[67];

    for (int i = 0; i < 67; i++)
    {
      xor16_src[i] = (unsigned short)(i * 5);
      xor16_mod2_src[i] = (unsigned short)(i * 9);
      xor8_src[i] = (uint8_t)(i * 3);
    }
    xor_u16(xor16, xor16_src, 67);
    xor_u16_mod2(xor16_mod2, xor16_mod2_src, 67);
    xor_u8(xor8, xor8_src, 67);
  }

  for (int i = 0; i < 67; i++)
  {
    if (xor16[i] != xor16_expected[i] ||
        xor16_mod2[i] != xor16_mod2_expected[i] ||
        xor8[i] != xor8_expected[i])
    {
      fprintf(stderr, "SIMD XOR mismatch at index %d\n", i);
      failed = 1;
      break;
    }
  }

  a.x[0] = 1;
  a.x[1] = 2;
  a.x[2] = 4;
  a.x[3] = 8;
  b.x[0] = 3;
  b.x[1] = 7;
  b.x[3] = 1;

  expected.x[0] = 2;
  expected.x[1] = 5;
  expected.x[2] = 4;
  expected.x[3] = 9;
  expected_product.x[0] = 3;
  expected_product.x[1] = 1;
  expected_product.x[2] = 2;
  expected_product.x[3] = 5;
  expected_product.x[4] = 58;
  expected_product.x[5] = 4;
  expected_product.x[6] = 8;
  expected_vector_product.x[0] = 3;
  expected_vector_product.x[1] = 1;
  expected_vector_product.x[2] = 2;
  expected_vector_product.x[3] = 28;

  failed |= expect_vec("vadd", vadd(a, b), expected.x);
  failed |= expect_vec("polynomial add", poly_add(a, b), expected.x);
  failed |= expect_vec("vmul", vmul(a, b), expected_vector_product.x);
  failed |= expect_vec("fmul", fmul(a, b), expected_product.x);
  failed |= expect_vec("polynomial multiply", poly_mul(a, b), expected_product.x);

  vec modulus = {0}, x = {0}, identity = {0};
  vec quotient = {0}, remainder = {0}, dividend = {0};
  modulus.x[0] = 1;
  modulus.x[1] = 1;
  modulus.x[3] = 1;
  x.x[1] = 1;
  identity.x[0] = 1;
  failed |= expect_vec("polynomial inverse",
                       poly_mod(poly_mul(x, poly_inverse(x, modulus)), modulus),
                       identity.x);

  quotient.x[0] = 2;
  quotient.x[2] = 1;
  remainder.x[0] = 3;
  remainder.x[1] = 1;
  dividend = vadd(poly_mul(quotient, modulus), remainder);
  failed |= expect_vec("polynomial division", poly_div(dividend, modulus), quotient.x);
  failed |= expect_vec("polynomial remainder", poly_mod(dividend, modulus), remainder.x);

  return failed;
}
