#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"


// 公開鍵を [I|T] に変換し、T (NB行 × TW語) を T_out に返す。
// 左 NB×NB ブロックが非正則なら -1(失敗した列は *fail_col)。
int pubkey_systematic(const unsigned short syndrome_matrix[N][K],
                      uint64_t *T_out, size_t *fail_col)
{
    uint64_t *A = calloc(NB * RW, sizeof *A);
    if (!A) return -1;

    // 列 i のビット c を A の行 c, ビット i に並べ替える(転置)
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < K; ++j) {
            unsigned v = syndrome_matrix[i][j];
            for (size_t b = 0; b < E; ++b)
                if ((v >> b) & 1)
                    A[(j * E + b) * RW + i / 64] |= (uint64_t)1 << (i % 64);
        }

    // 左 NB 列を単位行列にするガウス・ジョルダン消去
    for (size_t col = 0; col < NB; ++col) {
        size_t w = col / 64;
        uint64_t m = (uint64_t)1 << (col % 64);

        size_t piv = col;
        while (piv < NB && !(A[piv * RW + w] & m)) ++piv;
        if (piv == NB) { if (fail_col) *fail_col = col; free(A); return -1; }

        if (piv != col)
            for (size_t k = 0; k < RW; ++k) {
                uint64_t t = A[piv * RW + k];
                A[piv * RW + k] = A[col * RW + k];
                A[col * RW + k] = t;
            }

        const uint64_t *prow = A + col * RW;
        #pragma omp parallel for schedule(static)
        for (size_t r = 0; r < NB; ++r) {
            if (r == col) continue;
            uint64_t *row = A + r * RW;
            if (row[w] & m)
                for (size_t k = w; k < RW; ++k) row[k] ^= prow[k];
        }
    }

    // I を捨て、右 N-NB 列(語 52 から)だけを取り出す
    for (size_t r = 0; r < NB; ++r)
        memcpy(T_out + r * TW, A + r * RW + NB / 64, TW * sizeof *T_out);

    free(A);
    return 0;
}

// T の概要と先頭の行を表示する(全体は約2MBなので一部だけ)
void pubkey_systematic_print(const uint64_t *tt, size_t show_rows, size_t show_words)
{
    printf("T: %zu rows x %zu bits = %zu bytes\n",
           NB, TW * 64, NB * TW * sizeof *tt);
    for (size_t r = 0; r < show_rows && r < NB; ++r) {
        printf("row %4zu:", r);
        for (size_t k = 0; k < show_words && k < TW; ++k)
            printf(" %016llx", (unsigned long long)tt[r * TW + k]);
        printf(" ...\n");
    }
}

// 全体をバイナリで書き出す(1行 = 608 バイト、行ごとに連続)
int pubkey_systematic_write(const char *path, const uint64_t *tt)
{
    FILE *f = fopen(path, "wb");
    if (!f) return -1;
    size_t n = fwrite(tt, sizeof *tt, NB * TW, f);
    fclose(f);
    return n == NB * TW ? 0 : -1;
}