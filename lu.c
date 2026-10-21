#include "common.h"

static unsigned char packed_bit(const uint64_t *row, size_t bit)
{
  return (unsigned char)((row[bit / 64] >> (bit % 64)) & 1U);
}

int matrix_init(MTX *matrix, size_t rows, size_t cols)
{
  if (matrix == NULL || rows == 0 || cols == 0 ||
      rows > SIZE_MAX / cols ||
      rows * cols > SIZE_MAX / sizeof(*matrix->storage))
    return -1;

  MTX allocated = {0};
  allocated.x = calloc(rows, sizeof(*allocated.x));
  allocated.storage = calloc(rows * cols, sizeof(*allocated.storage));
  if (allocated.x == NULL || allocated.storage == NULL)
  {
    free(allocated.x);
    free(allocated.storage);
    return -1;
  }

  allocated.rows = rows;
  allocated.cols = cols;
  for (size_t row = 0; row < rows; row++)
    allocated.x[row] = allocated.storage + row * cols;

  matrix_free(matrix);
  *matrix = allocated;
  return 0;
}

void matrix_free(MTX *matrix)
{
  if (matrix == NULL)
    return;
  free(matrix->x);
  free(matrix->storage);
  memset(matrix, 0, sizeof(*matrix));
}

int matrix_copy(MTX *destination, const MTX *source)
{
  if (destination == NULL || source == NULL ||
      source->x == NULL || source->storage == NULL)
    return -1;
  if (destination == source)
    return 0;
  if (matrix_init(destination, source->rows, source->cols) != 0)
    return -1;
  memcpy(destination->storage, source->storage,
         source->rows * source->cols * sizeof(*source->storage));
  return 0;
}

int random_full_rank_binary_matrix(MTX *matrix, RNG *rng)
{
  /* Keep one reduced basis row per pivot, packed to avoid cubic bit loops. */
  if (matrix == NULL || matrix->x == NULL || matrix->storage == NULL ||
      matrix->rows == 0 || matrix->rows != matrix->cols || rng == NULL)
  {
    fprintf(stderr, "invalid or unavailable storage for binary scrambler matrix\n");
    return -1;
  }
  size_t n = matrix->rows;
  size_t words = (n + 63) / 64;
  if (n > SIZE_MAX / words || n * words > SIZE_MAX / sizeof(uint64_t))
    return -1;
  uint64_t *basis = calloc(n * words, sizeof(*basis));
  unsigned char *used = calloc(n, sizeof(*used));
  uint64_t *candidate = malloc(words * sizeof(*candidate));
  uint64_t *reduced = malloc(words * sizeof(*reduced));
  if (basis == NULL || used == NULL || candidate == NULL || reduced == NULL)
  {
    free(basis);
    free(used);
    free(candidate);
    free(reduced);
    fprintf(stderr, "failed to allocate binary scrambler generation workspace\n");
    return -1;
  }
  memset(matrix->storage, 0, matrix->rows * matrix->cols * sizeof(*matrix->storage));

  for (size_t row = 0; row < n; row++)
  {
    bool accepted = false;
    for (int attempt = 0; attempt < 4096 && !accepted; attempt++)
    {
      memset(candidate, 0, words * sizeof(*candidate));
      for (size_t col = 0; col < n; col++)
        candidate[col / 64] |= (uint64_t)(xor128(rng) & 1UL) << (col % 64);
      memcpy(reduced, candidate, words * sizeof(*reduced));

      for (size_t pivot = 0; pivot < n; pivot++)
      {
        if (!used[pivot])
          continue;
        if (packed_bit(reduced, pivot))
        {
          const uint64_t *basis_row = basis + pivot * words;
          for (size_t word = pivot / 64; word < words; word++)
            reduced[word] ^= basis_row[word];
        }
      }

      size_t pivot = n;
      for (size_t word = 0; word < words; word++)
      {
        uint64_t bits = reduced[word];
        if (bits != 0)
        {
          pivot = word * 64 + (size_t)__builtin_ctzll(bits);
          break;
        }
      }

      if (pivot == n)
        continue;

      memcpy(basis + pivot * words, reduced, words * sizeof(*reduced));
      used[pivot] = true;

      for (size_t col = 0; col < n; col++)
        matrix->x[row][col] = packed_bit(reduced, col);

      accepted = true;
    }

    if (!accepted)
    {
      memset(matrix->storage, 0,
             matrix->rows * matrix->cols * sizeof(*matrix->storage));
      free(basis);
      free(used);
      free(candidate);
      free(reduced);
      return -1;
    }
  }

  free(basis);
  free(used);
  free(candidate);
  free(reduced);
  return 0;
}

int is_reg(const MTX *cc_input, MTX *R, RNG *rng)
{
  (void)rng;
  if (cc_input == NULL || R == NULL || cc_input->x == NULL ||
      cc_input->rows == 0 || cc_input->rows != cc_input->cols ||
      R->x == NULL || R->rows != cc_input->rows ||
      R->cols != cc_input->cols)
  {
    fprintf(stderr, "invalid matrix dimensions for binary matrix inversion\n");
    return -1;
  }
  size_t n = cc_input->rows;
  size_t words = (n + 63) / 64;
  if (n > SIZE_MAX / words || n * words > SIZE_MAX / sizeof(uint64_t))
    return -1;
  uint64_t *reduced = calloc(n * words, sizeof(*reduced));
  uint64_t *inverse = calloc(n * words, sizeof(*inverse));
  if (reduced == NULL || inverse == NULL)
  {
    free(reduced);
    free(inverse);
    fprintf(stderr, "failed to allocate binary matrix inversion workspace\n");
    return -1;
  }
  memset(R->storage, 0, R->rows * R->cols * sizeof(*R->storage));

  for (size_t row = 0; row < n; row++)
  {
    for (size_t col = 0; col < n; col++)
    {
      unsigned short value = cc_input->x[row][col];
      if (value > 1)
      {
        free(reduced);
        free(inverse);
        return -1;
      }
      if (value != 0)
        reduced[row * words + col / 64] |= UINT64_C(1) << (col % 64);
    }
    inverse[row * words + row / 64] |= UINT64_C(1) << (row % 64);
  }

  for (size_t col = 0; col < n; col++)
  {
    size_t pivot = col;
    while (pivot < n && !packed_bit(reduced + pivot * words, col))
      pivot++;
    if (pivot == n)
    {
      free(reduced);
      free(inverse);
      return -1;
    }

    if (pivot != col)
    {
      for (size_t word = 0; word < words; word++)
      {
        uint64_t value = reduced[col * words + word];
        reduced[col * words + word] = reduced[pivot * words + word];
        reduced[pivot * words + word] = value;
        value = inverse[col * words + word];
        inverse[col * words + word] = inverse[pivot * words + word];
        inverse[pivot * words + word] = value;
      }
    }

    for (size_t row = 0; row < n; row++)
    {
      if (row == col || !packed_bit(reduced + row * words, col))
        continue;

      for (size_t word = 0; word < words; word++)
      {
        reduced[row * words + word] ^= reduced[col * words + word];
        inverse[row * words + word] ^= inverse[col * words + word];
      }
    }
  }

  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      R->x[row][col] = packed_bit(inverse + row * words, col);

  free(reduced);
  free(inverse);
  return 0;
}
