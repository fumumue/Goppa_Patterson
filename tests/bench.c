/*
 * tests/bench.c
 * 鍵生成・暗号化(シンドローム計算)・復号(Patterson)の各工程を計測するベンチマーク。
 *
 * 使い方:
 *   make bench                       (既定: 5回, シード1)
 *   ./build/bench 10                 (10回)
 *   ./build/bench 10 42              (10回, 乱数シード42)
 *   ./build/bench 10 42 bench.csv    (1回ごとの結果をCSVにも書き出す)
 *
 * 毎回、鍵生成 -> 暗号化 -> 復号 を往復し、誤り位置が一致するかを確認する。
 * 1回でも不一致があれば、終了コードを失敗(1)にする。
 */
#include "../common.h"
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* main.c の verify_decode_result と同じ判定(static のためここに複製) */
static int verify_decode_result(vec decoded,
                                const unsigned short permutation[N],
                                const unsigned short expected[N],
                                int *decoded_count,
                                int *expected_count,
                                int *matching_count)
{
    unsigned char decoded_positions[N] = {0};
    int valid = 1;

    *decoded_count = 0;
    *expected_count = 0;
    *matching_count = 0;

    for (int i = 0; i < T * 2; i++)
    {
        unsigned short position = decoded.x[i];
        if (i > 0 && position == 0)
            break;
        if (position >= N || decoded_positions[position])
        {
            valid = 0;
            continue;
        }
        decoded_positions[position] = 1;
        (*decoded_count)++;
    }

    for (int i = 0; i < N; i++)
    {
        int was_expected = expected[i] != 0;
        int was_decoded = decoded_positions[permutation[i]] != 0;
        *expected_count += was_expected;
        *matching_count += was_expected && was_decoded;
        if (was_expected != was_decoded)
            valid = 0;
    }
    return valid;
}

static void report(const char *name, const double *v, int n)
{
    double *s = malloc((size_t)n * sizeof(*s));
    double sum = 0.0;
    if (s == NULL)
        return;
    memcpy(s, v, (size_t)n * sizeof(*s));
    qsort(s, (size_t)n, sizeof(*s), cmp_double);
    for (int i = 0; i < n; i++)
        sum += v[i];
    double median = (n % 2) ? s[n / 2] : (s[n / 2 - 1] + s[n / 2]) / 2.0;
    printf("%-10s %10.1f %10.1f %10.1f %10.1f\n",
           name, s[0] * 1e3, median * 1e3, sum / n * 1e3, s[n - 1] * 1e3);
    free(s);
}

int main(int argc, char **argv)
{
    int iters = (argc > 1) ? atoi(argv[1]) : 5;
    unsigned int seed = (argc > 2) ? (unsigned int)strtoul(argv[2], NULL, 10) : 1u;
    const char *csv_path = (argc > 3) ? argv[3] : NULL;
    int pass = 0, fail = 0;

    if (iters < 1)
        iters = 1;

    unsigned short *z1 = calloc(N, sizeof(*z1));
    unsigned short *permutation = malloc((size_t)N * sizeof(*permutation));
    unsigned short (*syndrome_matrix)[K] = calloc(N, sizeof(*syndrome_matrix));
    double *t_key = calloc((size_t)iters, sizeof(double));
    double *t_enc = calloc((size_t)iters, sizeof(double));
    double *t_dec = calloc((size_t)iters, sizeof(double));
    double *t_all = calloc((size_t)iters, sizeof(double));
    RNG rng = {123456789, 362436069, 521288629, 88675123};
    FILE *csv = NULL;

    if (!z1 || !permutation || !syndrome_matrix || !t_key || !t_enc || !t_dec || !t_all)
    {
        fprintf(stderr, "bench: failed to allocate buffers\n");
        return EXIT_FAILURE;
    }
    if (K > N)
    {
        fprintf(stderr, "configuration error: K is too big\n");
        return EXIT_FAILURE;
    }

    srand(seed);
    gen_gf(E, N, 1);

    if (csv_path != NULL)
    {
        csv = fopen(csv_path, "w");
        if (csv != NULL)
            fprintf(csv, "iter,keygen_s,encrypt_s,decode_s,total_s,verify\n");
    }

    printf("N=%d K=%d E=%d T=%d iterations=%d seed=%u\n", N, K, E, T, iters, seed);
#ifdef _OPENMP
    printf("OpenMP max threads: %d\n", omp_get_max_threads());
#endif

    for (int it = 0; it < iters; it++)
    {
        MTX inverse_scrambler = {0};
        vec w, f, v;
        double t0, t1, t2, t3;

        memset(z1, 0, (size_t)N * sizeof(*z1));

        /* 鍵生成 */
        t0 = now_sec();
        w = pubkeygen(permutation, &inverse_scrambler, &rng, syndrome_matrix);
        t1 = now_sec();

        /* 暗号化(誤りベクトル生成 + シンドローム計算) */
        if (mkerr(z1, K, &rng) != 0)
        {
            fprintf(stderr, "bench: mkerr failed (iter %d)\n", it);
            matrix_free(&inverse_scrambler);
            fail++;
            break;
        }
        f = zin(&inverse_scrambler, syndrome_matrix, z1);
        t2 = now_sec();

        /* 復号 */
        v = patterson(w, f);
        t3 = now_sec();

        /* 正しさの確認(時間には含めない) */
        int dc = 0, ec = 0, mc = 0;
        unsigned short ya[N] = {0};
        int ok = verify_decode_result(v, permutation, z1, &dc, &ec, &mc);
        int reported = ero2(v, permutation, ya);
        ok = ok && reported == dc;

        t_key[it] = t1 - t0;
        t_enc[it] = t2 - t1;
        t_dec[it] = t3 - t2;
        t_all[it] = t3 - t0;
        if (ok)
            pass++;
        else
            fail++;

        printf("iter %2d: keygen %.3fs  encrypt %.3fs  decode %.3fs  [%s]\n",
               it, t_key[it], t_enc[it], t_dec[it], ok ? "PASS" : "FAIL");
        if (csv != NULL)
            fprintf(csv, "%d,%.6f,%.6f,%.6f,%.6f,%s\n",
                    it, t_key[it], t_enc[it], t_dec[it], t_all[it], ok ? "PASS" : "FAIL");

        matrix_free(&inverse_scrambler);
    }

    int done = pass + fail;
    if (done > 0)
    {
        printf("\n%-10s %10s %10s %10s %10s\n", "stage", "min(ms)", "median(ms)", "mean(ms)", "max(ms)");
        report("keygen", t_key, done);
        report("encrypt", t_enc, done);
        report("decode", t_dec, done);
        report("total", t_all, done);
    }
    printf("\nverify: PASS %d / FAIL %d\n", pass, fail);

    if (csv != NULL)
        fclose(csv);
    free(z1);
    free(permutation);
    free(syndrome_matrix);
    free(t_key);
    free(t_enc);
    free(t_dec);
    free(t_all);
    return fail == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
