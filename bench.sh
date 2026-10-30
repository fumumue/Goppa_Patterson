#!/usr/bin/env bash
# compare_bench.sh
# 自作ベンチ(tests/bench.c)と libmceliece のベンチを、条件をそろえて順に実行し、
# 結果をファイルに保存する。ベンチのプログラム自体は書き直さなくてよい。
#
# 使い方(Goppa_Patterson のフォルダで):
#   chmod +x compare_bench.sh
#   ./compare_bench.sh
#
# 設定は環境変数で変えられる:
#   BENCH_CMD  自作ベンチを動かすコマンド   (既定: "make bench")
#   BENCH_LM   libmceliece 非f版のベンチ    (既定: ./bench_lm)
#   BENCH_LMF  libmceliece f版のベンチ      (既定: ./bench_lmf)
#   RUNS_LM    libmceliece側の測定回数      (既定: 20)
#   CPU1       1コア実行で使うCPU番号       (既定: 0)
#   CPUS2      2コア実行で使うCPU番号       (既定: 0,1)
#   ROUNDS     全体を何周するか             (既定: 1)
#
# 例: BENCH_CMD="./tests/bench" ROUNDS=2 ./compare_bench.sh
#
# 注意: コード内で omp_set_num_threads() などを呼んでいる場合、
#       環境変数 OMP_NUM_THREADS は無視される。そのため、コア数の制御は
#       taskset(CPUの固定)で行う。1コアに固定すれば、スレッドが何本あっても
#       実際に使えるのは1コアだけになる。

set -u

BENCH_CMD="${BENCH_CMD:-make bench}"
BENCH_LM="${BENCH_LM:-./bench_lm}"
BENCH_LMF="${BENCH_LMF:-./bench_lmf}"
RUNS_LM="${RUNS_LM:-20}"
CPU1="${CPU1:-0}"
CPUS2="${CPUS2:-0,1}"
ROUNDS="${ROUNDS:-1}"
OUT="results_$(date +%Y%m%d_%H%M%S).txt"

if ! command -v taskset >/dev/null 2>&1; then
    echo "taskset が見つかりません(util-linux が必要です)" >&2
    exit 1
fi

run() {
    local title="$1"; shift
    {
        echo "=== $title ==="
        echo "\$ $*"
    } | tee -a "$OUT"
    "$@" 2>&1 | tee -a "$OUT"
    echo | tee -a "$OUT"
}

{
    echo "# compare_bench  $(date '+%F %T')"
    echo "# CPU: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ //')"
    echo "# CPU番号とコア番号の対応(同じCOREの行はSMTの兄弟):"
    lscpu -e=CPU,CORE 2>/dev/null | sed 's/^/#   /'
    if [ -r /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor ]; then
        echo "# governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor)"
    fi
    echo "# 1コア用CPU=$CPU1  2コア用CPU=$CPUS2  ROUNDS=$ROUNDS"
    echo
} | tee -a "$OUT"

for r in $(seq 1 "$ROUNDS"); do
    echo "########## round $r / $ROUNDS ##########" | tee -a "$OUT"

    # 自作: 1コアに固定
    run "自作 / 1コア固定 (CPU $CPU1)" \
        env OMP_NUM_THREADS=1 taskset -c "$CPU1" bash -c "$BENCH_CMD"

    # 自作: 2コアに固定
    run "自作 / 2コア固定 (CPU $CPUS2)" \
        env OMP_NUM_THREADS=2 taskset -c "$CPUS2" bash -c "$BENCH_CMD"

    # 自作: 固定なし(普段の実行)
    run "自作 / 固定なし(普段どおり)" \
        bash -c "$BENCH_CMD"

    # libmceliece は常に1コア
    if [ -x "$BENCH_LMF" ]; then
        run "libmceliece 8192128f / 1コア固定 (CPU $CPU1)" \
            taskset -c "$CPU1" "$BENCH_LMF" "$RUNS_LM"
    else
        echo "($BENCH_LMF がないのでスキップ)" | tee -a "$OUT"
    fi
    if [ -x "$BENCH_LM" ]; then
        run "libmceliece 8192128 / 1コア固定 (CPU $CPU1)" \
            taskset -c "$CPU1" "$BENCH_LM" "$RUNS_LM"
    else
        echo "($BENCH_LM がないのでスキップ)" | tee -a "$OUT"
    fi
done

echo "結果を保存しました: $OUT"
