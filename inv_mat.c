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


// mulmat_gf2 より前に置く
static inline void pack_elems(const unsigned short *elem, size_t kcnt,
                              uint64_t *brow)   // brow は 0 クリア済み
{
    for (size_t j = 0; j < kcnt; ++j) {
        size_t pos = j * E;
        size_t w = pos / 64, s = pos % 64;
        uint64_t v = elem[j] & ((1u << E) - 1);
        brow[w] |= v << s;
        if (s + E > 64)
            brow[w + 1] |= v >> (64 - s);
    }
}


#define CH 8

static int mulmat_gf2(const MTX *A, const MTX *B,
                      const unsigned short (*elems)[K],   // NULL なら B から読む
                      size_t rows, size_t cols, size_t inner, MTX *result)
{
    const size_t win  = (inner + 63) / 64;
    const size_t wout = (cols  + 63) / 64;
    const size_t nch  = (inner + CH - 1) / CH;
    const size_t TS   = (size_t)1 << CH;

    uint64_t *at  = calloc(inner * wout, sizeof *at);
    uint64_t *tbl = calloc(nch * TS * wout, sizeof *tbl);
    if (!at || !tbl) { free(at); free(tbl); return -1; }

    for (size_t c = 0; c < cols; ++c)
        for (size_t k = 0; k < inner; ++k) {
            if (A->x[c][k] > 1) { free(at); free(tbl); return -1; }
            at[k * wout + c / 64] |= (uint64_t)A->x[c][k] << (c % 64);
        }

    for (size_t ch = 0; ch < nch; ++ch) {
        uint64_t *t = tbl + ch * TS * wout;
        for (unsigned v = 1; v < TS; ++v) {
            unsigned low = v & (0u - v);
            size_t k = ch * CH + (size_t)__builtin_ctz(low);
            const uint64_t *prev = t + (size_t)(v ^ low) * wout;
            for (size_t w = 0; w < wout; ++w)
                t[(size_t)v * wout + w] =
                    prev[w] ^ (k < inner ? at[k * wout + w] : 0);
        }
    }

    int bad = 0;
    #pragma omp parallel for schedule(static) reduction(|:bad)
    for (size_t r = 0; r < rows; ++r) {
        uint64_t brow[win], acc[wout];
        memset(brow, 0, sizeof brow);
        memset(acc, 0, sizeof acc);

        if (elems) {
            pack_elems(elems[r], inner / E, brow);
        } else {
            for (size_t k = 0; k < inner; ++k) {
                unsigned short b = B->x[r][k];
                if (b > 1) { bad = 1; break; }
                brow[k / 64] |= (uint64_t)b << (k % 64);
            }
            if (bad) continue;
        }

        for (size_t ch = 0; ch < nch; ++ch) {
            size_t bit = ch * CH;
            unsigned v = (unsigned)(brow[bit / 64] >> (bit % 64)) & (unsigned)(TS - 1);
            const uint64_t *t = tbl + (ch * TS + v) * wout;
            for (size_t w = 0; w < wout; ++w) acc[w] ^= t[w];
        }
        for (size_t c = 0; c < cols; ++c)
            result->x[r][c] = (acc[c / 64] >> (c % 64)) & 1;
    }
    free(at); free(tbl);
    return bad ? -1 : 0;
}

int pubkey_pack(const MTX *scrambler,
                unsigned short field_matrix[N][K],
                const unsigned short permutation[N],
                unsigned short syndrome_matrix[N][K])
{
    const size_t nbits = (size_t)K * E;            // 3328
    const size_t words = (nbits + 63) / 64;        // 52
    const size_t nch   = (nbits + CH - 1) / CH;
    const size_t TS    = (size_t)1 << CH;

    if (scrambler == NULL || scrambler->x == NULL ||
        scrambler->rows != nbits || scrambler->cols != nbits)
        return -1;

    uint64_t *at  = calloc(nbits * words, sizeof *at);
    uint64_t *tbl = calloc(nch * TS * words, sizeof *tbl);
    uint64_t *out = calloc((size_t)N * words, sizeof *out);
    if (!at || !tbl || !out) { free(at); free(tbl); free(out); return -1; }

    // S^T を詰める: at[k] の第 c ビット = S[c][k]
    for (size_t c = 0; c < nbits; ++c)
        for (size_t k = 0; k < nbits; ++k) {
            if (scrambler->x[c][k] > 1) { free(at); free(tbl); free(out); return -1; }
            at[k * words + c / 64] |= (uint64_t)scrambler->x[c][k] << (c % 64);
        }

    // 表: tbl[ch][v] = その塊で v のビットが立つ行の XOR
    for (size_t ch = 0; ch < nch; ++ch) {
        uint64_t *t = tbl + ch * TS * words;
        for (unsigned v = 1; v < TS; ++v) {
            unsigned low = v & (0u - v);
            size_t k = ch * CH + (size_t)__builtin_ctz(low);
            const uint64_t *prev = t + (size_t)(v ^ low) * words;
            for (size_t w = 0; w < words; ++w)
                t[(size_t)v * words + w] =
                    prev[w] ^ (k < nbits ? at[k * words + w] : 0);
        }
    }

    // 各行: 要素 → 詰めたビット → 表引きで XOR
    #pragma omp parallel for schedule(static)
    for (size_t r = 0; r < N; ++r) {
        uint64_t brow[words];
        uint64_t *acc = out + r * words;
        memset(brow, 0, sizeof brow);
        pack_elems(field_matrix[r], K, brow);
        for (size_t ch = 0; ch < nch; ++ch) {
            size_t bit = ch * CH;
            unsigned v = (unsigned)(brow[bit / 64] >> (bit % 64)) & (unsigned)(TS - 1);
            const uint64_t *t = tbl + (ch * TS + v) * words;
            for (size_t w = 0; w < words; ++w) acc[w] ^= t[w];
        }
    }

    // permutation で並べ替えながら 13 ビット要素に戻す(toByte の代わり)
    int bad = 0;
    for (size_t i = 0; i < N; ++i) {
        if (permutation[i] >= N) { bad = 1; break; }
        const uint64_t *o = out + (size_t)permutation[i] * words;
        for (size_t j = 0; j < K; ++j) {
            size_t pos = j * E, w = pos / 64, s = pos % 64;
            uint64_t v = o[w] >> s;
            if (s + E > 64) v |= o[w + 1] << (64 - s);
            syndrome_matrix[i][j] = (unsigned short)(v & (((uint64_t)1 << E) - 1));
        }
    }
    free(at); free(tbl); free(out);
    return bad ? -1 : 0;
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

    if (mode == 1 || mode == 2)
        return mulmat_gf2(A, B, NULL, rows, cols, inner, result); 
    
    // mode 3 は元のまま
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


int mulmat2(const MTX *A, const MTX *B, int mode, MTX *result)
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
    
    //#pragma omp parallel for schedule(static)
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
