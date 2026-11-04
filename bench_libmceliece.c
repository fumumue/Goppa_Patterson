/*
 * bench_libmceliece.c
 * libmceliece (Classic McEliece) の mceliece8192128 を、
 * 自作実装の tests/bench.c と同じ形式で計測するベンチマーク。
 *
 * ビルド(通常版 / f版):
 *   gcc -O2 -o bench_lm   bench_libmceliece.c -lmceliece
 *   gcc -O2 -DUSE_F -o bench_lmf bench_libmceliece.c -lmceliece
 *   (リンクエラーが出たら -lrandombytes や -lcpucycles を追加)
 *
 * 実行(1コアに固定して比較するのがおすすめ):
 *   taskset -c 0 ./bench_lm 5
 *
 * 計測内容: keygen / encrypt(enc) / decode(dec) と、その合計。
 * 時間は CLOCK_MONOTONIC の実時間(ms)。
 */
#define _POSIX_C_SOURCE 200809L
#include <mceliece.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef USE_F
#define KEYPAIR mceliece8192128f_keypair
#define ENC     mceliece8192128f_enc
#define DEC     mceliece8192128f_dec
#define NAME    "mceliece8192128f"
#else
#define KEYPAIR mceliece8192128_keypair
#define ENC     mceliece8192128_enc
#define DEC     mceliece8192128_dec
#define NAME    "mceliece8192128"
#endif

/* f版と非f版は鍵・暗号文のサイズが同じなので、サイズ定数は共通 */
#define PKBYTES mceliece8192128_PUBLICKEYBYTES
#define SKBYTES mceliece8192128_SECRETKEYBYTES
#define CTBYTES mceliece8192128_CIPHERTEXTBYTES
#define KBYTES  mceliece8192128_BYTES

#define STAGES 4 /* keygen, encrypt, decode, total */

static double now_ms(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1e3 + ts.tv_nsec / 1e6;
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

int main(int argc, char **argv)
{
    int runs = (argc > 1) ? atoi(argv[1]) : 5;
    if (runs < 1) runs = 1;

    /* 公開鍵は約1.3MB。スタックではなくヒープに置く */
    unsigned char *pk = malloc(PKBYTES);
    unsigned char *sk = malloc(SKBYTES);
    unsigned char ct[CTBYTES];
    unsigned char k_enc[KBYTES], k_dec[KBYTES];
    double (*t)[STAGES] = malloc(sizeof(double[STAGES]) * (size_t)runs);
    if (!pk || !sk || !t) {
        fprintf(stderr, "malloc failed\n");
        return 1;
    }

    /* ウォームアップ(初回のページフォルトや遅延初期化を除外) */
    KEYPAIR(pk, sk);
    ENC(ct, k_enc, pk);
    DEC(k_dec, ct, sk);

    int pass = 0, fail = 0;
    for (int i = 0; i < runs; i++) {
        double t0 = now_ms();
        KEYPAIR(pk, sk);
        double t1 = now_ms();
        ENC(ct, k_enc, pk);
        double t2 = now_ms();
        DEC(k_dec, ct, sk);
        double t3 = now_ms();

        t[i][0] = t1 - t0;
        t[i][1] = t2 - t1;
        t[i][2] = t3 - t2;
        t[i][3] = t3 - t0;

        if (memcmp(k_enc, k_dec, KBYTES) == 0) pass++; else fail++;
    }

    static const char *names[STAGES] = {"keygen", "encrypt", "decode", "total"};
    printf("%s  (runs=%d)\n", NAME, runs);
    printf("%-8s %10s %10s %10s %10s\n",
           "stage", "min(ms)", "median(ms)", "mean(ms)", "max(ms)");
    for (int s = 0; s < STAGES; s++) {
        double *v = malloc(sizeof(double) * (size_t)runs);
        double sum = 0;
        for (int i = 0; i < runs; i++) { v[i] = t[i][s]; sum += v[i]; }
        qsort(v, (size_t)runs, sizeof(double), cmp_double);
        double med = (runs % 2) ? v[runs / 2]
                                : (v[runs / 2 - 1] + v[runs / 2]) / 2.0;
        printf("%-8s %10.1f %10.1f %10.1f %10.1f\n",
               names[s], v[0], med, sum / runs, v[runs - 1]);
        free(v);
    }
    printf("\nverify: %s %d / FAIL %d\n", fail ? "FAIL" : "PASS", pass, fail);

    free(pk); free(sk); free(t);
    return fail ? 1 : 0;
}