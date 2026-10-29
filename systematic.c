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

// systematic.c
#ifndef NB
#define NB ((size_t)K * E)            // 3328 : シンドローム長(行数)
#endif
#ifndef RW
#define RW ((size_t)N / 64)           // 128  : 長さ N の行の語数
#endif
#ifndef TW
#define TW ((size_t)(N - NB) / 64)    // 76   : T の1行の語数
#endif
_Static_assert(((size_t)K * E) % 64 == 0, "K*E must be a multiple of 64");
_Static_assert(N % 64 == 0, "N must be a multiple of 64");

// H(field_matrix と permutation)から、S を使わずに直接 [I|T] を作り、
// I を除いた T を tmat に返す。
//   field_matrix[r] : サポート r の列(K 個の 13 ビット要素)
//   permutation[i]  : 位置 i が使うサポート番号
//   tmat            : NB 行 × TW 語(呼び出し側で NB*TW*8 バイト確保)
//   colperm[j]      : 並べ替え後の位置 j → 元の位置 i
//                     前半 NB 個がピボット列、後半 N-NB 個が T の列
//   rank_out        : 到達した階数(NULL 可)
// 戻り値: 0 成功 / -1 失敗(階数が NB 未満、または不正な入力、メモリ不足)
int pubkey_systematic_pivot_H(const unsigned short field_matrix[N][K],
                              const unsigned short permutation[N],
                              uint64_t *tmat, uint16_t colperm[N],
                              size_t *rank_out)
{
    if (rank_out) *rank_out = 0;

    for (size_t i = 0; i < N; ++i)
        if (permutation[i] >= N) return -1;

    uint64_t *hm = calloc(NB * RW, sizeof *hm);
    uint16_t *pivcol = malloc(NB * sizeof *pivcol);
    unsigned char *is_piv = calloc(N, 1);
    if (!hm || !pivcol || !is_piv) {
        free(hm); free(pivcol); free(is_piv);
        return -1;
    }

    // 転置: 64 列ずつのブロックで並列化(ブロックごとに書き込む語が別)
    #pragma omp parallel for schedule(static)
    for (size_t blk = 0; blk < RW; ++blk) {
        for (size_t t = 0; t < 64; ++t) {
            const size_t i = blk * 64 + t;
            const unsigned short *colv = field_matrix[permutation[i]];
            for (size_t j = 0; j < K; ++j) {
                unsigned v = colv[j];
                for (size_t b = 0; b < E; ++b)
                    if ((v >> b) & 1)
                        hm[(j * E + b) * RW + blk] |= (uint64_t)1 << t;
            }
        }
    }

    // 消去: 並列領域は 1 回だけ作る
    size_t rank = 0, col = 0, piv_row = 0;
    int piv_found = 0;

    #pragma omp parallel 
    {
        while(1) {
            // ピボット探索と行の入れ替えは 1 スレッドで行う
            #pragma omp single
            {
                piv_found = 0;
                while (col < N && rank < NB) {
                    const size_t w = col / 64;
                    const uint64_t m = (uint64_t)1 << (col % 64);
                    size_t p = rank;
                    while (p < NB && !(hm[p * RW + w] & m)) ++p;
                    if (p < NB) {
                        if (p != rank)
                            for (size_t k = w; k < RW; ++k) {
                                uint64_t t = hm[p * RW + k];
                                hm[p * RW + k] = hm[rank * RW + k];
                                hm[rank * RW + k] = t;
                            }
                        piv_row = rank;
                        piv_found = 1;
                        break;                 // この col がピボット列
                    }
                    ++col;                     // 独立でない列は飛ばす
                }
            }   // single の終わりで全スレッドが待ち合わせる

            if (!piv_found) break;             // 全スレッドが同じ値を見て抜ける

            const size_t w = col / 64;
            const uint64_t m = (uint64_t)1 << (col % 64);
            const uint64_t *prow = hm + piv_row * RW;

            #pragma omp for schedule(static)
            for (size_t r = 0; r < NB; ++r) {
                if (r == piv_row) continue;
                uint64_t *row = hm + r * RW;
                if (row[w] & m)
                    for (size_t k = w; k < RW; ++k) row[k] ^= prow[k];
            }   // for の終わりでも待ち合わせる

            #pragma omp single
            {
                pivcol[rank++] = (uint16_t)col;
                is_piv[col] = 1;
                ++col;
            }
        }
    }

    if (rank_out) *rank_out = rank;
    if (rank < NB) {                           // 本当の階数不足
        free(hm); free(pivcol); free(is_piv);
        return -1;
    }

    // 並べ替え: ピボット列 → 残りの列(昇順)
    size_t pos = 0;
    for (size_t r = 0; r < NB; ++r) colperm[pos++] = pivcol[r];
    for (size_t c = 0; c < N; ++c)
        if (!is_piv[c]) colperm[pos++] = (uint16_t)c;

    // T を詰めて取り出す(I は捨てる)
    memset(tmat, 0, NB * TW * sizeof *tmat);
    #pragma omp parallel for schedule(static)
    for (size_t r = 0; r < NB; ++r) {
        for (size_t k = 0; k < N - NB; ++k) {
            size_t orig = colperm[NB + k];
            uint64_t bit = (hm[r * RW + orig / 64] >> (orig % 64)) & 1;
            tmat[r * TW + k / 64] |= bit << (k % 64);
        }
    }

    free(hm); free(pivcol); free(is_piv);
    return 0;
}


