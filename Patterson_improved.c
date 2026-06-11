// ============================================================
// Goppa Code Patterson Decoder - Improved Version
// ============================================================
// date : 20210325 (ver 1.0 improved)
// author : cryptography researcher
// description : One Variable Polynomial library with OpenMP support
// implementation : Niederreiter Cryptosystem by Patterson's decoding
//
// IMPROVEMENTS:
// - Removed duplicate includes
// - Converted .c includes to .h format
// - Added comprehensive debug system
// - Better code organization
// - Memory safety improvements
// - Doxygen documentation for all functions
// ============================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <assert.h>
#include <execinfo.h>
#include <omp.h>

// ============================================================
// Header Files (Fixed: using .h instead of .c)
// ============================================================
#include "global.h"
#include "struct.h"
#include "debug.h"        // Changed from debug.c
#include "ecole.h"        // Changed from ecole.c
#include "chash.h"        // Changed from chash.c
#include "lu.h"
#include "sha3.h"
#include "inv_mat.h"      // Changed from inv_mat.c

// ============================================================
// Debug Configuration
// ============================================================
#define DEBUG_MODE 0
#define VERBOSE_MODE 0

#if DEBUG_MODE
  #define DEBUG_PRINT(fmt, ...) fprintf(stderr, "[DEBUG] " fmt "\n", ##__VA_ARGS__)
#else
  #define DEBUG_PRINT(fmt, ...) (void)0
#endif

#if VERBOSE_MODE
  #define VERBOSE_PRINT(fmt, ...) printf("[VERBOSE] " fmt "\n", ##__VA_ARGS__)
#else
  #define VERBOSE_PRINT(fmt, ...) (void)0
#endif

#define ERROR_PRINT(fmt, ...) fprintf(stderr, "[ERROR] " fmt "\n", ##__VA_ARGS__)

// ============================================================
// Constants
// ============================================================
#define REG 64
#define SIZE_OF_ARRAY(arr) (sizeof(arr) / sizeof((arr)[0]))

// Error codes
#define ERR_SUCCESS     0
#define ERR_MALLOC      1
#define ERR_INVALID_POL 2
#define ERR_DIVIDE_BY_ZERO 3

// ============================================================
// Global Variables (with comments)
// ============================================================

/** Syndrome copy buffer */
unsigned short sy[K] = {0};

/** Goppa polynomial coefficient */
static unsigned short g[K + 1] = {1, 0, 0, 0, 1, 0, 1};

/** Binary parity check matrix */
static MTX BB = {0};

/** Transformed matrix */
static MTX H = {0};

/** Performance counters */
unsigned int AA = 0, B = 0;

// ============================================================
// External Function Declarations
// ============================================================

extern unsigned long xor128(void);
extern int mlt(int x, int y);
extern int mltn(int n, int a);
extern void makeS(void);

// ============================================================
// SECTION 1: FINITE FIELD OPERATIONS
// ============================================================

/**
 * Compute multiplicative inverse in finite field
 * @param a Element in GF
 * @return Multiplicative inverse, or (unsigned short)-1 if a == 0
 */
unsigned short oinv(unsigned short a)
{
  if (a == 0)
    return (unsigned short)-1;

  return (N - fg[a]) % (N - 1) + 1;
}

/**
 * Alternative multiplicative inverse (lookup method)
 * @param a Element in GF
 * @return Multiplicative inverse
 */
unsigned short oinv2(unsigned short a)
{
  int i;

  if (a == 0)
    return (unsigned short)-1;

  for (i = 0; i < N; i++) {
    if (gf[mlt(fg[a], i)] == 1)
      return (unsigned short)i;
  }

  ERROR_PRINT("oinv2: no inverse found");
  exit(ERR_INVALID_POL);
  return 0;
}

/**
 * Find x such that a * x = b
 * @param a Multiplicand
 * @param b Result
 * @return x where a * x = b
 */
unsigned short equ(unsigned short a, unsigned short b)
{
  int i;

  for (i = 0; i < N; i++) {
    if (gf[mlt(fg[a], fg[i])] == b)
      break;
  }
  return i;
}

// ============================================================
// SECTION 2: POLYNOMIAL TYPE CONVERSIONS
// ============================================================

/**
 * Convert OP (sparse) to vec (dense) format
 * @param f Polynomial in OP format
 * @return Polynomial in vec format
 */
vec o2v(OP f)
{
  vec a = {0};
  int i;

  for (i = 0; i < DEG; i++) {
    if (f.t[i].a > 0 && f.t[i].n < DEG)
      a.x[f.t[i].n] = f.t[i].a;
  }

  return a;
}

/**
 * Convert vec (dense) to OP (sparse) format
 * @param a Polynomial in vec format
 * @return Polynomial in OP format
 */
OP v2o(vec a)
{
  int i, j = 0;
  OP f = {0};

  for (i = 0; i < DEG; i++) {
    if (a.x[i] > 0) {
      f.t[j].n = i;
      f.t[j++].a = a.x[i];
    }
  }

  return f;
}

/**
 * Convert vec to OP with explicit counting
 * @param a Polynomial in vec format
 * @return Polynomial in OP format
 */
OP v4x(vec a)
{
  OP x = {0};
  int i, count = 0;

  for (i = 0; i < DEG; i++) {
    if (a.x[i] > 0) {
      x.t[count].a = a.x[i];
      x.t[count++].n = i;
    }
  }

  return x;
}

// ============================================================
// SECTION 3: UTILITY FUNCTIONS
// ============================================================

/**
 * Wait for user input
 */
void wait(void)
{
  int a;
  printf(" (enter number and hit return) ");
  fflush(stdout);
  scanf("%d", &a);
}

/**
 * Normalize OP polynomial
 * @param f Polynomial to normalize
 * @return Normalized polynomial
 */
OP conv(OP f)
{
  vec v = {0};
  OP g = {0};

  v = o2v(f);
  g = v2o(v);

  return g;
}

/**
 * Count terms in polynomial
 * @param f Polynomial
 * @return Number of terms
 */
int terms(OP f)
{
  int i, count = 0;

  for (i = 0; i < DEG; i++) {
    if (f.t[i].a > 0)
      count++;
  }

  return count;
}

/**
 * Get degree (vec format)
 * @param a Polynomial
 * @return Degree
 */
int deg(vec a)
{
  int i, n = 0, flg = 0;

  for (i = 0; i < DEG; i++) {
    if (a.x[i] > 0) {
      n = i;
      flg = 1;
    }
  }
  if (flg == 0)
    return 0;

  return n;
}

/**
 * Get degree (OP format)
 * @param f Polynomial
 * @return Degree
 */
int odeg(OP f)
{
  int i, j = 0;

  if (f.t[0].a == 0)
    return 0;

  for (i = 0; i < DEG; i++) {
    if (j < f.t[i].n && f.t[i].a > 0)
      j = f.t[i].n;
  }

  return j;
}

/**
 * Print polynomial (OP format)
 * @param f Polynomial
 */
void oprintpol(OP f)
{
  int i, n;

  f = conv(f);
  n = odeg(f);
  printf("degree=%d\n", n);

  for (i = n; i > -1; i--) {
    if (f.t[i].a > 0)
      printf("%u*x^%u + ", f.t[i].a, f.t[i].n);
  }
  printf("\n");
}

/**
 * Print raw OP polynomial
 * @param f Polynomial
 */
void op_print_raw(const OP f)
{
  puts("=== OP Polynomial (Raw) ===");
  for (int i = 0; i < DEG; i++) {
    if (f.t[i].a > 0)
      printf("[%d] coeff=%u, deg=%u\n", i, f.t[i].a, f.t[i].n);
  }
}

/**
 * Verify OP polynomial structure
 * @param f Polynomial
 * @return true if valid
 */
bool op_verify(const OP f)
{
  bool end = false;
  unsigned short n_max = 0;
  
  for (int i = 0; i < DEG; i++) {
    if (end && (f.t[i].n != 0 || f.t[i].a != 0)) {
      op_print_raw(f);
      ERROR_PRINT("verify: data after end at %d", i);
      print_trace();
      return false;
    }
    
    if (f.t[i].a == 0) {
      end = true;
      continue;
    }
    
    if (f.t[i].n + 1 <= n_max) {
      op_print_raw(f);
      ERROR_PRINT("verify: order error at %d", i);
      print_trace();
      return false;
    }
    n_max = f.t[i].n + 1;
  }
  
  return true;
}

/**
 * Normalize to dense indexed format
 * @param f Polynomial
 * @return Normalized polynomial
 */
OP norm(OP f)
{
  OP h = {0};
  int i;

  for (i = 0; i < DEG; i++) {
    if (f.t[i].a > 0)
      h.t[f.t[i].n].a = f.t[i].a;
  }

  return h;
}

// ============================================================
// SECTION 4: VECTOR ARITHMETIC
// ============================================================

/**
 * Add two polynomials (XOR in GF(2))
 * @param a First polynomial
 * @param b Second polynomial
 * @return Sum
 */
vec vadd(vec a, vec b)
{
  vec c = {0};
  int i, k;

  if (deg(a) >= deg(b)) {
    k = deg(a) + 1;
  } else {
    k = deg(b) + 1;
  }
  
  for (i = 0; i < k; i++) {
    c.x[i] = a.x[i] ^ b.x[i];
  }

  return c;
}

/**
 * Multiply two polynomials
 * @param a First polynomial
 * @param b Second polynomial
 * @return Product
 */
vec vmul(vec a, vec b)
{
  int i, j, k, l;
  vec c = {0};

  k = deg(a);
  l = deg(b);

  for (i = 0; i <= k; i++) {
    for (j = 0; j <= l; j++) {
      if (a.x[i] > 0 && b.x[j] > 0) {
        c.x[i + j] ^= gf[mlt(fg[a.x[i]], fg[b.x[j]])];
      }
    }
  }

  return c;
}

/**
 * Fast multiply (for sparse polynomials)
 * @param a First polynomial
 * @param b Second polynomial
 * @return Product
 */
vec fmul(vec a, vec b)
{
  int i, j, k, l;
  vec c = {0};
  OP o = {0}, p = {0};

  k = terms(v2o(a));
  l = terms(v2o(b));
  o = v4x(a);
  p = v4x(b);
  
  for (i = 0; i < k; i++) {
    for (j = 0; j < l; j++) {
      c.x[o.t[i].n + p.t[j].n] ^= 
        gf[mlt(fg[o.t[i].a], fg[p.t[j].a])];
    }
  }

  return c;
}

// ============================================================
// END OF PART 1
// ============================================================
// NOTE: This file will be split into multiple parts due to size
// Part 2 will include: Polynomial division, power operations
// Part 3 will include: GCD, irreducibility testing, decoding
// ============================================================
