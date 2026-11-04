// bench_lm.c : libmceliece (mceliece8192128 / mceliece8192128f) の鍵生成・暗号化・復号の時間を測る。
//
// ビルド(libmceliece をインストール済みの場合):
//   gcc -O3 -march=native -o bench_lm  bench_lm.c -lmceliece
//   gcc -O3 -march=native -DUSE_F -o bench_lmf bench_lm.c -lmceliece
//
// libmceliece を /usr/local などに入れた場合は、必要に応じて次も付ける:
//   -I/usr/local/include -L/usr/local/lib -Wl,-rpath,/usr/local/lib
//
// 使い方:
//   ./bench_lm            (測定回数は環境変数 RUNS_LM、既定 20)
//   RUNS_LM=5 taskset -c 0 ./bench_lm
//
// 注意: 関数名と定数名は libmceliece のヘッダ(mceliece.h)にあるものを使っている。
//       コンパイルエラーが出たら、ヘッダで実際の名前を確認すること。
//         grep -n "8192128" $(find / -name mceliece.h 2>/dev/null | head -1) | head -20

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mceliece.h>

#ifdef USE_F
#define NAME(x) mceliece8192128f_##x
#define LABEL "mceliece8192128f"
#else
#define NAME(x) mceliece8192128_##x
#define LABEL "mceliece8192128"
#endif

static double now_ms(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static void report(const char *name, double *v, int n)
{
    qsort(v, n, sizeof(double), cmp_double);
    double median = (n % 2) ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2.0;
    printf("%-8s 中央値 %9.3f ms (最小 %9.3f, 最大 %9.3f)\n", name, median, v[0], v[n - 1]);
}

int main(void)
{
    int runs = 20;
    const char *env = getenv("RUNS_LM");
    if (env && atoi(env) > 0)
        runs = atoi(env);

    unsigned char *pk = malloc(NAME(PUBLICKEYBYTES));
    unsigned char *sk = malloc(NAME(SECRETKEYBYTES));
    unsigned char *ct = malloc(NAME(CIPHERTEXTBYTES));
    unsigned char key1[NAME(BYTES)], key2[NAME(BYTES)];
    double *tk = malloc(sizeof(double) * runs);
    double *te = malloc(sizeof(double) * runs);
    double *td = malloc(sizeof(double) * runs);
    double *tt = malloc(sizeof(double) * runs);
    if (!pk || !sk || !ct || !tk || !te || !td || !tt)
    {
        fprintf(stderr, "malloc failed\n");
        return 1;
    }

    int pass = 0, fail = 0;
    for (int r = 0; r < runs; r++)
    {
        double t0 = now_ms();
        // 戻り値の型(void か int か)がバージョンで違っても通るよう、戻り値は使わない。
        // 成否は、暗号化側と復号側の鍵が一致するかどうかで判定する。
        NAME(keypair)(pk, sk);
        double t1 = now_ms();
        NAME(enc)(ct, key1, pk);
        double t2 = now_ms();
        NAME(dec)(key2, ct, sk);
        double t3 = now_ms();

        tk[r] = t1 - t0;
        te[r] = t2 - t1;
        td[r] = t3 - t2;
        tt[r] = t3 - t0;
        if (memcmp(key1, key2, sizeof key1) == 0)
            pass++;
        else
            fail++;
    }

    printf("%s, %d 回\n", LABEL, runs);
    report("keygen", tk, runs);
    report("encrypt", te, runs);
    report("decode", td, runs);
    report("total", tt, runs);
    printf("verify PASS %d / FAIL %d\n", pass, fail);
    return fail ? 1 : 0;
}