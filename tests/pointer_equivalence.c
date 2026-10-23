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
  OP noncanonical = {0};
  unsigned short xor16[67], xor16_expected[67], xor16_mod2[67], xor16_mod2_expected[67];
  uint8_t xor8[67], xor8_expected[67];
  MTX full_rank_matrix = {0};
  MTX singular_matrix = {0}, inverse_matrix = {0};
  MTX product_a = {0}, product_b = {0}, product = {0};
  RNG rng = {123456789, 362436069, 521288629, 88675123};
  OP f, h;
  int failed = 0;

  gen_gf(E, N, 1);
  if (matrix_init(&full_rank_matrix, 8, 8) != 0 ||
      matrix_init(&inverse_matrix, 8, 8) != 0 ||
      matrix_init(&singular_matrix, 2, 2) != 0)
  {
    fprintf(stderr, "matrix allocation failed\n");
    return 1;
  }
  if (random_full_rank_binary_matrix(&full_rank_matrix, &rng) != 0)
  {
    fprintf(stderr, "full-rank matrix generation failed\n");
    return 1;
  }

  if (is_reg(&full_rank_matrix, &inverse_matrix, &rng) != 0)
  {
    fprintf(stderr, "is_reg failed to invert a full-rank matrix\n");
    failed = 1;
  }
  else
  {
    for (int row = 0; row < 8; row++)
    {
      for (int col = 0; col < 8; col++)
      {
        unsigned short product = 0;
        for (int k = 0; k < 8; k++)
          product ^= full_rank_matrix.x[row][k] & inverse_matrix.x[k][col];
        if (product != (unsigned short)(row == col))
        {
          fprintf(stderr, "matrix inverse mismatch at [%d][%d]\n", row, col);
          failed = 1;
        }
      }
    }
  }
  singular_matrix.x[0][0] = 1;
  singular_matrix.x[0][1] = 1;
  singular_matrix.x[1][0] = 1;
  singular_matrix.x[1][1] = 1;
  MTX singular_inverse = {0};
  if (matrix_init(&singular_inverse, 2, 2) != 0)
  {
    fprintf(stderr, "singular inverse allocation failed\n");
    failed = 1;
  }
  else if (is_reg(&singular_matrix, &singular_inverse, &rng) == 0)
  {
    fprintf(stderr, "is_reg accepted a singular matrix\n");
    failed = 1;
  }
  matrix_free(&full_rank_matrix);
  matrix_free(&singular_matrix);
  matrix_free(&singular_inverse);
  matrix_free(&inverse_matrix);

  if (matrix_init(&product_a, 2, 2) != 0 ||
      matrix_init(&product_b, 3, 2) != 0)
  {
    fprintf(stderr, "matrix multiplication allocation failed\n");
    failed = 1;
  }
  else
  {
    product_a.x[0][0] = 1;
    product_a.x[0][1] = 1;
    product_a.x[1][0] = 1;
    product_b.x[0][0] = 1;
    product_b.x[1][1] = 1;
    product_b.x[2][0] = 1;
    if (mulmat(&product_a, &product_b, 1, &product) != 0 ||
        product.rows != 3 || product.cols != 2 ||
        product.x[0][0] != 1 || product.x[0][1] != 1 ||
        product.x[1][0] != 1 || product.x[1][1] != 0 ||
        product.x[2][0] != 1 || product.x[2][1] != 1)
    {
      fprintf(stderr, "binary matrix multiplication mismatch\n");
      failed = 1;
    }
  }
  matrix_free(&product_a);
  matrix_free(&product_b);
  matrix_free(&product);

  MTX target_scrambler = {0}, target_inverse = {0};
  if (matrix_init(&target_scrambler, DEG, DEG) != 0 ||
      matrix_init(&target_inverse, DEG, DEG) != 0)
  {
    fprintf(stderr, "target-size matrix allocation failed\n");
    failed = 1;
  }
  else if (random_full_rank_binary_matrix(&target_scrambler, &rng) != 0 ||
           is_reg(&target_scrambler, &target_inverse, &rng) != 0)
  {
    fprintf(stderr, "target-size binary matrix generation/inversion failed\n");
    failed = 1;
  }
  else
  {
    const size_t sample_indices[] = {0, DEG / 2, DEG - 1};
    for (size_t i = 0; i < SIZE_OF_ARRAY(sample_indices); i++)
      for (size_t j = 0; j < SIZE_OF_ARRAY(sample_indices); j++)
      {
        unsigned short product_value = 0;
        size_t row = sample_indices[i], col = sample_indices[j];
        for (size_t k = 0; k < DEG; k++)
          product_value ^= target_scrambler.x[row][k] & target_inverse.x[k][col];
        if (product_value != (unsigned short)(row == col))
        {
          fprintf(stderr, "target-size matrix inverse mismatch at [%zu][%zu]\n",
                  row, col);
          failed = 1;
        }
      }
  }
  matrix_free(&target_scrambler);
  matrix_free(&target_inverse);

  MTX target_check = {0};
  if (matrix_init(&target_check, N, DEG) != 0)
  {
    fprintf(stderr, "target-size rectangular matrix allocation failed\n");
    failed = 1;
  }
  else
  {
    target_check.x[N - 1][DEG - 1] = 1;
    if (target_check.storage[(size_t)N * DEG - 1] != 1)
    {
      fprintf(stderr, "target-size contiguous matrix indexing mismatch\n");
      failed = 1;
    }
  }
  matrix_free(&target_check);

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

  noncanonical.t[0] = (oterm){2, 5};
  noncanonical.t[1] = (oterm){0, 7};
  noncanonical.t[2] = (oterm){2, 9};
  noncanonical.t[4] = (oterm){3, 1};
  OP other_noncanonical = {0};
  other_noncanonical.t[0] = (oterm){1, 4};
  other_noncanonical.t[2] = (oterm){2, 3};
  other_noncanonical.t[3] = (oterm){0, 6};
  failed |= expect_vec("conv", o2v(conv(noncanonical)), (unsigned short[DEG]){
      [0] = 7, [2] = 9, [3] = 1});
  failed |= expect_vec("conv canonical", o2v(conv(v2o(a))), a.x);
  failed |= expect_vec("oadd noncanonical",
                       o2v(oadd(noncanonical, other_noncanonical)),
                       (unsigned short[DEG]){[0] = 1, [1] = 4, [2] = 10, [3] = 1});

  for (int sample = 0; sample < 32; sample++)
  {
    OP left = {0}, right = {0};
    vec expected_sum = {0};

    for (int i = 0; i < DEG; i++)
    {
      left.t[i].n = (unsigned short)((i * 37 + sample * 11) % (DEG + 8));
      left.t[i].a = (unsigned short)((i * 13 + sample) % M);
      right.t[i].n = (unsigned short)((i * 19 + sample * 23) % (DEG + 8));
      right.t[i].a = (unsigned short)((i * 7 + sample * 3) % M);
    }

    expected_sum = o2v(conv(left));
    vec right_sum = o2v(conv(right));
    for (int i = 0; i < DEG; i++)
      expected_sum.x[i] ^= right_sum.x[i];

    failed |= expect_vec("oadd normalized", o2v(oadd(left, right)), expected_sum.x);
  }

  failed |= expect_vec("vadd", vadd(a, b), expected.x);
  failed |= expect_vec("oadd", o2v(oadd(v2o(a), v2o(b))), expected.x);

  f = v2o(a);
  h = v2o(b);
  failed |= expect_vec("vmul", vmul(a, b), expected_vector_product.x);
  failed |= expect_vec("fmul", fmul(a, b), expected_product.x);
  failed |= expect_vec("vmul vs fmul reference", vmul(a, b), fmul(a, b).x);

  /* Differential tests: vmul must agree with the full-range vector reference. */
  for (int sample = 0; sample < 32; sample++)
  {
    vec left = {0}, right = {0};
    for (int i = 0; i < 20; i++)
    {
      left.x[i] = (unsigned short)((sample * 17 + i * 13 + 1) % M);
      right.x[i] = (unsigned short)((sample * 29 + i * 7 + 1) % M);
    }
    failed |= expect_vec("vmul randomized reference", vmul(left, right),
                         fmul(left, right).x);
  }
  failed |= expect_vec("omul", o2v(omul(f, h)), expected_product.x);

  for (int sample = 0; sample < 32; sample++)
  {
    vec left = {0}, right = {0}, product = {0};

    for (int i = 0; i < 20; i++)
    {
      left.x[i] = (unsigned short)((sample * 17 + i * 13 + 1) % M);
      right.x[i] = (unsigned short)((sample * 29 + i * 7 + 1) % M);
    }
    for (int i = 0; i < 20; i++)
      for (int j = 0; j < 20; j++)
        if (left.x[i] != 0 && right.x[j] != 0)
          product.x[i + j] ^= gf[mlt(fg[left.x[i]], fg[right.x[j]])];

    failed |= expect_vec("omul vector multiplication",
                         o2v(omul(v2o(left), v2o(right))), product.x);
  }

  return failed;
}
