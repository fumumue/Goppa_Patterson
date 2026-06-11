// ============================================================
// Patterson_improved.c - PART 2
// Polynomial Division and Power Operations
// ============================================================

// ============================================================
// SECTION 5: POLYNOMIAL ADDITION (OP FORMAT)
// ============================================================

/**
 * Add two polynomials (OP format)
 * @param f First polynomial
 * @param g Second polynomial
 * @return Sum f + g
 */
OP oadd(OP f, OP g)
{
  f = conv(f);
  g = conv(g);

  vec a = {0}, b = {0}, c = {0};
  int i;
  OP h = {0};

  a = o2v(f);
  b = o2v(g);

  for (i = 0; i < DEG; i++) {
    c.x[i] = a.x[i] ^ b.x[i];
  }
  
  h = v2o(c);
  return h;
}

/**
 * Sort polynomial terms in ascending order
 * @param f Polynomial to sort
 * @return Sorted polynomial
 */
OP sort(OP f)
{
  oterm o = {0};
  int i, j, k;

  k = terms(f);
  for (i = 0; i < k + 1; ++i) {
    for (j = i + 1; j < k + 1; ++j) {
      if (f.t[i].n > f.t[j].n) {
        o = f.t[i];
        f.t[i] = f.t[j];
        f.t[j] = o;
      }
    }
  }

  return f;
}

/**
 * Extract leading term from OP polynomial
 * @param f Polynomial
 * @return Leading (highest degree) term
 */
oterm oLT(OP f)
{
  int i, k;
  oterm s = {0};

  k = terms(f);
  if (k == 0) {
    s.a = 0;
    s.n = 0;
    return s;
  }

  s = f.t[0];
  for (i = 0; i < k + 1; i++) {
    if (f.t[i].a > 0 && f.t[i].n > s.n) {
      s.n = f.t[i].n;
      s.a = f.t[i].a;
    }
  }

  return s;
}

// ============================================================
// SECTION 6: TERM MULTIPLICATION
// ============================================================

/**
 * Multiply polynomial by single term (vec format)
 * @param f Polynomial
 * @param t Single term (coeff * x^degree)
 * @return Result f * t
 */
vec vterml(vec f, oterm t)
{
  int i;
  vec h = {0};

  for (i = 0; i < DEG - t.n; i++) {
    if (f.x[i] > 0)
      h.x[i + t.n] = gf[mlt(fg[f.x[i]], fg[t.a])];
  }

  return h;
}

/**
 * Multiply polynomial by single term (OP format)
 * @param f Polynomial
 * @param t Single term
 * @return Result f * t
 */
OP oterml(OP f, oterm t)
{
  f = conv(f);
  
  int i;
  OP h = {0};

  for (i = 0; i < DEG; i++) {
    if (f.t[i].a > 0) {
      h.t[i].n = f.t[i].n + t.n;
      h.t[i].a = gf[mlt(fg[f.t[i].a], fg[t.a])];
    }
  }

  h = conv(h);
  return h;
}

/**
 * Multiply two polynomials (OP format)
 * @param f First polynomial
 * @param g Second polynomial
 * @return Product f * g
 */
OP omul(OP f, OP g)
{
  f = conv(f);
  g = conv(g);
  
  int i, k, l;
  oterm t = {0};
  OP h = {0}, e = {0};

  k = odeg(f);
  l = odeg(g);

  for (i = 0; i <= k; i++) {
    if (g.t[i].a > 0) {
      t = g.t[i];
      e = oterml(f, t);
      h = oadd(h, e);
    }
  }
  
  return h;
}

// ============================================================
// SECTION 7: POLYNOMIAL DIVISION
// ============================================================

/**
 * Extract leading term from vec format
 * @param f Polynomial
 * @return Leading term
 */
oterm lt(vec f)
{
  int i;
  oterm t = {0};

  for (i = DEG - 1; i >= 0; i--) {
    if (f.x[i] > 0) {
      t.n = i;
      t.a = f.x[i];
      break;
    }
  }

  return t;
}

/**
 * Extract leading term from OP format (alternative)
 * @param f Polynomial
 * @return Leading term
 */
oterm LT(OP f)
{
  int i;
  oterm t = {0};

  for (i = DEG - 1; i >= 0; i--) {
    if (f.t[i].a > 0) {
      t.n = f.t[i].n;
      t.a = f.t[i].a;
      break;
    }
  }

  return t;
}

/**
 * Extract first/constant term
 * @param f Polynomial
 * @return First term
 */
oterm LT2(OP f)
{
  oterm t = {0};
  t.n = f.t[0].n;
  t.a = f.t[0].a;
  return t;
}

/**
 * Divide leading term of f by term t
 * @param f Polynomial dividend
 * @param t Term divisor
 * @return Quotient term: LT(f) / t
 */
oterm LTdiv(OP f, oterm t)
{
  oterm tt = {0}, s = {0};

  tt = LT(f);
  
  if (tt.n < t.n) {
    s.n = 0;
    s.a = 0;
  }
  else if (tt.n == t.n) {
    s.n = 0;
    s.a = equ(t.a, tt.a);
  }
  else if (tt.n > t.n) {
    s.n = tt.n - t.n;
    s.a = equ(t.a, tt.a);
  }
  else if (t.n == 0 && t.a > 0) {
    s.a = gf[mlt(fg[tt.a], oinv(t.a))];
    s.n = tt.n;
  }

  return s;
}

/**
 * Make polynomial monic (leading coeff = 1)
 * @param f Polynomial
 * @param d Leading coefficient to divide by
 * @return Monic polynomial
 */
OP coeff(OP f, unsigned short d)
{
  int i, k;

  f = conv(f);
  k = odeg(f) + 1;
  
  for (i = 0; i < k; i++)
    if (f.t[i].a > 0)
      f.t[i].a = gf[mlt(fg[f.t[i].a], oinv(d))];

  return f;
}

/**
 * Print polynomial (vec format)
 * @param a Polynomial to print
 */
void printpol(vec a)
{
  int i, n;

  n = deg(a);

  for (i = n; i >= 0; i--) {
    if (a.x[i] > 0)
      printf("%u*x^%d + ", a.x[i], i);
  }
  printf("\n");
}

/**
 * Polynomial division (quotient only)
 * @param f Dividend polynomial
 * @param g Divisor polynomial
 * @return Quotient q where f = g*q + r
 */
OP odiv(OP f, OP g)
{
  f = conv(f);
  g = conv(g);
  
  int i = 0;
  OP h = {0}, tt = {0};
  oterm b = {0}, c = {0};

  if (LT(g).a == 0) {
    ERROR_PRINT("odiv: division by zero polynomial");
    exit(ERR_DIVIDE_BY_ZERO);
  }
  
  if (LT(g).n == 0 && LT(g).a > 1)
    return coeff(f, LT(g).a);

  b = LT(g);
  
  if (b.a == 1 && b.n == 0)
    return f;
  
  if (odeg(f) < odeg(g))
    return f;

  while (LT(f).n >= LT(g).n && LT(g).a > 0) {
    c = LTdiv(f, b);
    
    tt.t[i] = c;
    i++;

    h = oterml(g, c);
    f = oadd(f, h);
    
    if (odeg(f) == 0 || odeg(g) == 0)
      break;

    if (c.n == 0)
      break;
  }

  // Reverse quotient terms
  OP ret = {0};
  int tt_terms = terms(tt);
  for (i = 0; i < tt_terms; i++) {
    ret.t[i] = tt.t[tt_terms - i - 1];
  }
  ret = conv(ret);
  
  return ret;
}

/**
 * Polynomial modulo (remainder only)
 * @param f Dividend polynomial
 * @param g Divisor polynomial
 * @return Remainder r where f = g*q + r, deg(r) < deg(g)
 */
OP omod(OP f, OP g)
{
  int i = 0;
  OP h = {0};
  oterm b = {0}, c = {0};

  if (LT(f).n < LT(g).n)
    return f;

  b = LT(g);

  while (LT(f).n >= LT(g).n && LT(g).a > 0) {
    c = LTdiv(f, b);
    h = oterml(g, c);
    f = oadd(f, h);
    
    if (odeg(f) == 0 || odeg(g) == 0)
      break;

    if (c.n == 0 || b.n == 0)
      break;
  }

  return f;
}

// ============================================================
// SECTION 8: POWER OPERATIONS
// ============================================================

/**
 * Compute polynomial power (vec format)
 * @param f Polynomial
 * @param n Exponent
 * @return f^n
 */
vec vpow(vec f, int n)
{
  int i;
  vec g = {0};

  g = f;

  for (i = 1; i < n; i++)
    g = vmul(g, f);

  return g;
}

/**
 * Compute polynomial power (OP format)
 * @param f Polynomial
 * @param n Exponent
 * @return f^n
 */
OP opow(OP f, int n)
{
  int i;
  OP g = {0};

  g = f;

  for (i = 1; i < n; i++)
    g = omul(g, f);

  return g;
}

/**
 * Compute f^n mod g (OP format, binary method)
 * @param f Base polynomial
 * @param mod Modulus polynomial
 * @param n Exponent
 * @return f^n mod g
 */
OP opowmod1(OP f, OP mod, int n)
{
  vec v = {0};
  OP ret = {0};

  v.x[0] = 1;
  ret = v2o(v);
  
  while (n > 0) {
    if (n & 1)
      ret = omod(omul(ret, f), mod);
    f = omod(omul(f, f), mod);
    n >>= 1;
  }
  
  return ret;
}

/**
 * Compute polynomial power modulo
 * @param f Base polynomial
 * @param mod Modulus polynomial
 * @param n Exponent
 * @return f^n mod g
 */
OP opowmod(OP f, OP mod, int n)
{
  int i = 1;

  while (i < n + 1) {
    f = omod(omul(f, f), mod);
    i++;
  }

  return f;
}

// ============================================================
// SECTION 9: POLYNOMIAL EVALUATION (TRACE)
// ============================================================

/**
 * Evaluate polynomial at x (vec format)
 * @param f Polynomial
 * @param x Evaluation point
 * @return f(x)
 */
unsigned short vtrace(vec f, unsigned short x)
{
  int i, d;
  unsigned short u = 0;

  d = deg(f);

  for (i = 0; i < d + 1; i++) {
    u ^= gf[mlt(fg[f.x[i]], mltn(i, fg[x]))];
  }

  return u;
}

/**
 * Evaluate polynomial at x (OP format)
 * @param f Polynomial
 * @param x Evaluation point
 * @return f(x)
 */
unsigned short trace(OP f, unsigned short x)
{
  int i, d;
  unsigned short u = 0;

  d = deg(o2v(f));

  for (i = 0; i < d + 1; i++) {
    u ^= gf[mlt(fg[f.t[i].a], mltn(f.t[i].n, fg[x]))];
  }

  return u;
}

// ============================================================
// SECTION 10: UTILITY CONVERSIONS
// ============================================================

/**
 * Convert bit array to integer
 * @param b Bit array
 * @return Integer value
 */
unsigned short b2B(unsigned short b[E])
{
  int i;
  unsigned short a = 0;

  for (i = E - 1; i > -1; i--)
    a ^= (b[E - i - 1] << i);

  return a;
}

/**
 * Convert integer to vector
 * @param n Integer to convert
 * @return Vector representation
 */
vec i2v(unsigned short n)
{
  vec v = {0};
  int i = 0;

  while (n > 0) {
    v.x[i++] = n % 2;
    n = (n >> 1);
  }

  return v;
}

/**
 * Convert vector to integer
 * @param v Vector to convert
 * @return Integer representation
 */
unsigned short v2i(vec v)
{
  unsigned int d = 0, i, e = 0;

  for (i = 0; i < deg(v) + 1; i++) {
    e = v.x[i];
    d ^= (e << i);
  }

  return d;
}

/**
 * Convert 64-bit integer to vector
 * @param n Integer to convert
 * @return Vector representation
 */
vec ull2v(unsigned long long int n)
{
  vec v = {0};
  int i = 0;

  while (n > 0) {
    v.x[i++] = n % 2;
    n = (n >> 1);
  }

  return v;
}

/**
 * Convert vector to 64-bit integer
 * @param v Vector to convert
 * @return 64-bit integer representation
 */
unsigned long long int v2ull(vec v)
{
  unsigned long long int d = 0, i, e = 0;

  for (i = 0; i < deg(v) + 1; i++) {
    e = v.x[i];
    d ^= (e << i);
  }

  return d;
}

/**
 * Set vector from array
 * @param n Number of elements
 * @return Vector set from array
 */
vec Setvec(int n)
{
  int i;
  vec v = {0};

  for (i = 0; i < n; i++) {
    v.x[n - 1 - i] = c[i];
  }

  return v;
}

/**
 * Print vector
 * @param v Vector to print
 */
void printvec(vec v)
{
  int i;

  for (i = 0; i < deg(v) + 1; i++) {
    printf("[%d]=%d\n", i, v.x[i]);
  }
}

/**
 * Compute integer power
 * @param q Base
 * @param u Exponent
 * @return q^u
 */
unsigned int ipow(unsigned int q, unsigned int u)
{
  unsigned int i, m = 1;

  for (i = 0; i < u; i++)
    m *= q;

  DEBUG_PRINT("ipow: %u^%u = %u", q, u, m);

  return m;
}

// ============================================================
// END OF PART 2
// ============================================================
// Part 3 includes: Extended GCD, decoding algorithms
// ============================================================
