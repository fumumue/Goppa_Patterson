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
  OP f, h;
  int failed = 0;

  gen_gf(E, N, 1);

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
  failed |= expect_vec("omul", o2v(omul(f, h)), expected_product.x);

  return failed;
}
