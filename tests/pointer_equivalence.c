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
  OP f, h;
  int failed = 0;

  gen_gf(E, N, 1);

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
  failed |= expect_vec("oadd", o2v(oadd(v2o(a), v2o(b))), expected.x);

  f = v2o(a);
  h = v2o(b);
  failed |= expect_vec("vmul", vmul(a, b), expected_vector_product.x);
  failed |= expect_vec("fmul", fmul(a, b), expected_product.x);
  failed |= expect_vec("omul", o2v(omul(f, h)), expected_product.x);

  return failed;
}
