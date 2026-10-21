#include <stdio.h>
#include <stdlib.h>
#include "common.h"

#define F (K * E)

static unsigned short parity_product(const uint64_t *left,
                                     const uint64_t *right,
                                     size_t words)
{
    unsigned int parity = 0;
    for (size_t word = 0; word < words; ++word)
        parity ^= (unsigned int)__builtin_parityll(left[word] & right[word]);
    return (unsigned short)parity;
}

void rp(unsigned short *a)
{
    unsigned short *pool = malloc((size_t)N * sizeof(*pool));
    if (pool == NULL) {
        perror("rp: malloc");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < N; ++i)
        pool[i] = (unsigned short)i;
    for (int i = 0; i < N; ++i) {
        unsigned int j = (unsigned int)rand() % (unsigned int)(N - i);
        a[i] = pool[j];
        pool[j] = pool[N - i - 1];
    }
    free(pool);
}

int Inv(unsigned short b)
{
    return (b == 0) ? 0 : (M - fg[b]);
}

static int matrix_inverse(const MTX *input, int n, MTX *inverse)
{
    if (input == NULL || inverse == NULL || n <= 0 ||
        input->rows < (size_t)n || input->cols < (size_t)n)
        return -1;

    size_t dimension = (size_t)n;
    if (dimension > SIZE_MAX / dimension ||
        dimension * dimension > SIZE_MAX / sizeof(unsigned short))
        return -1;
    size_t cells = dimension * dimension;
    unsigned short *left = malloc(cells * sizeof(*left));
    unsigned short *right = calloc(cells, sizeof(*right));
    if (left == NULL || right == NULL) {
        free(left);
        free(right);
        return -1;
    }
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c)
            left[(size_t)r * n + c] = input->x[r][c];
        right[(size_t)r * n + r] = 1;
    }

    for (int col = 0; col < n; ++col) {
        int pivot = col;
        while (pivot < n && left[(size_t)pivot * n + col] == 0)
            ++pivot;
        if (pivot == n) {
            free(left);
            free(right);
            return 0;
        }
        if (pivot != col) {
            for (int c = 0; c < n; ++c) {
                unsigned short tmp = left[(size_t)col * n + c];
                left[(size_t)col * n + c] = left[(size_t)pivot * n + c];
                left[(size_t)pivot * n + c] = tmp;
                tmp = right[(size_t)col * n + c];
                right[(size_t)col * n + c] = right[(size_t)pivot * n + c];
                right[(size_t)pivot * n + c] = tmp;
            }
        }

        unsigned short scale = gf[Inv(fg[left[(size_t)col * n + col]])];
        for (int c = 0; c < n; ++c) {
            left[(size_t)col * n + c] =
                (unsigned short)mlt(left[(size_t)col * n + c], scale);
            right[(size_t)col * n + c] =
                (unsigned short)mlt(right[(size_t)col * n + c], scale);
        }

        for (int r = 0; r < n; ++r) {
            if (r == col)
                continue;
            unsigned short factor = left[(size_t)r * n + col];
            if (factor == 0)
                continue;
            for (int c = 0; c < n; ++c) {
                left[(size_t)r * n + c] ^=
                    (unsigned short)mlt(factor, left[(size_t)col * n + c]);
                right[(size_t)r * n + c] ^=
                    (unsigned short)mlt(factor, right[(size_t)col * n + c]);
            }
        }
    }

    if (matrix_init(inverse, (size_t)n, (size_t)n) != 0) {
        free(left);
        free(right);
        return -1;
    }
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c)
            inverse->x[r][c] = right[(size_t)r * n + c];
    free(left);
    free(right);
    return 1;
}

int gauss(const MTX *input, MTX *inverse)
{
    return matrix_inverse(input, F, inverse);
}

int matinv(const MTX *input, int n, MTX *inverse)
{
    return matrix_inverse(input, n, inverse);
}

int mulmat(const MTX *A, const MTX *B, int mode, MTX *result)
{
    if (A == NULL || B == NULL || result == NULL ||
        result == A || result == B ||
        A->x == NULL || B->x == NULL)
        return -1;

    size_t rows, cols, inner;
    if (mode == 1) {
        if (A->rows != A->cols || A->cols != B->cols)
            return -1;
        rows = B->rows;
        cols = A->rows;
        inner = A->cols;
    } else if (mode == 2) {
        size_t half = (size_t)E * (K / 2 + 1);
        if (A->rows < half || A->cols < half ||
            B->rows < N || B->cols < half)
            return -1;
        rows = N;
        cols = half;
        inner = half;
    } else if (mode == 3) {
        if (A->cols != B->rows)
            return -1;
        rows = A->rows;
        cols = B->cols;
        inner = A->cols;
    } else {
        return -1;
    }

    if (matrix_init(result, rows, cols) != 0)
        return -1;
    if (mode == 1 || mode == 2) {
        size_t words = (inner + 63) / 64;
        if (cols > SIZE_MAX / words || rows > SIZE_MAX / words ||
            cols * words > SIZE_MAX / sizeof(uint64_t) ||
            rows * words > SIZE_MAX / sizeof(uint64_t))
            return -1;
        uint64_t *packed_a = calloc(cols * words, sizeof(*packed_a));
        uint64_t *packed_b = calloc(rows * words, sizeof(*packed_b));
        if (packed_a == NULL || packed_b == NULL) {
            free(packed_a);
            free(packed_b);
            return -1;
        }
        for (size_t r = 0; r < cols; ++r) {
            for (size_t c = 0; c < inner; ++c) {
                if (A->x[r][c] > 1) {
                    free(packed_a);
                    free(packed_b);
                    return -1;
                }
                packed_a[r * words + c / 64] |=
                    (uint64_t)A->x[r][c] << (c % 64);
            }
        }
        for (size_t r = 0; r < rows; ++r) {
            for (size_t c = 0; c < inner; ++c) {
                if (B->x[r][c] > 1) {
                    free(packed_a);
                    free(packed_b);
                    return -1;
                }
                packed_b[r * words + c / 64] |=
                    (uint64_t)B->x[r][c] << (c % 64);
            }
        }
        for (size_t r = 0; r < rows; ++r) {
            for (size_t c = 0; c < cols; ++c) {
                const uint64_t *a_row = packed_a + c * words;
                const uint64_t *b_row = packed_b + r * words;
                result->x[r][c] = parity_product(a_row, b_row, words);
            }
        }
        free(packed_a);
        free(packed_b);
        return 0;
    }

    for (size_t r = 0; r < rows; ++r) {
        for (size_t c = 0; c < cols; ++c) {
            unsigned short value = 0;
            for (size_t k = 0; k < inner; ++k)
                value ^= (unsigned short)mlt(A->x[r][k], B->x[k][c]);
            result->x[r][c] = value;
        }
    }
    return 0;
}

void mmul(const MTX *A, const MTX *B)
{
    if (A == NULL || B == NULL || A->cols != B->rows) {
        fprintf(stderr, "mmul: incompatible matrix dimensions\n");
        return;
    }
    for (size_t r = 0; r < A->rows; ++r) {
        for (size_t c = 0; c < B->cols; ++c) {
            unsigned short value = 0;
            for (size_t k = 0; k < A->cols; ++k)
                value ^= (unsigned short)mlt(A->x[r][k], B->x[k][c]);
            printf("%u%c", value, c + 1 == B->cols ? '\n' : ' ');
        }
    }
}
