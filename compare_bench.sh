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
#
# 実行順(1周あたり):
#   1. 自作       1コア  (taskset -c $CPU1)
#   2. 自作       2コア  (taskset -c $CPUS2)
#   3. libmceliece 非f版  1コア
#   4. libmceliece f版    1コア
# 周回数が2以上のときは、この順を繰り返す(測定条件のずれを均すため)。
# 結果は results/compare_YYYYmmdd_HHMMSS.txt に保存し、画面にも表示する。

set -u

BENCH_CMD="${BENCH_CMD:-make bench}"
BENCH_LM="${BENCH_LM:-./bench_lm}"
BENCH_LMF="${BENCH_LMF:-./bench_lmf}"
RUNS_LM="${RUNS_LM:-20}"
CPU1="${CPU1:-0}"
CPUS2="${CPUS2:-0,1}"
ROUNDS="${ROUNDS:-1}"

if ! command -v taskset >/dev/null 2>&1; then
    echo "taskset が見つかりません(util-linux が必要です)。" >&2
    exit 1
fi

mkdir -p results
OUT="results/compare_$(date +%Y%m%d_%H%M%S).txt"

# 実行ファイルの存在確認(自作は BENCH_CMD の先頭の語が make 以外のときだけ確認)
check_exec() {
    local label="$1" path="$2"
    if [ ! -x "$path" ]; then
        echo "警告: $label が見つからない、または実行できません: $path" >&2
        return 1
    fi
    return 0
}

run_block() {
    local title="$1"
    shift
    {
        echo
        echo "=================================================="
        echo "## $title"
        echo "## コマンド: $*"
        echo "=================================================="
    } | tee -a "$OUT"
    "$@" 2>&1 | tee -a "$OUT"
}

{
    echo "compare_bench.sh  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "host: $(uname -srm)"
    if command -v lscpu >/dev/null 2>&1; then
        lscpu | grep -E 'Model name|^CPU\(s\)|Thread|Core' | sed 's/^/  /'
    fi
    echo "BENCH_CMD=$BENCH_CMD"
    echo "BENCH_LM=$BENCH_LM  BENCH_LMF=$BENCH_LMF  RUNS_LM=$RUNS_LM"
    echo "CPU1=$CPU1  CPUS2=$CPUS2  ROUNDS=$ROUNDS"
} | tee "$OUT"

for round in $(seq 1 "$ROUNDS"); do
    echo | tee -a "$OUT"
    echo "######## 周回 $round / $ROUNDS ########" | tee -a "$OUT"

    # 1. 自作 1コア
    run_block "自作 1コア (CPU $CPU1) 周回$round" \
        taskset -c "$CPU1" bash -c "$BENCH_CMD"

    # 2. 自作 2コア
    run_block "自作 2コア (CPU $CPUS2) 周回$round" \
        taskset -c "$CPUS2" bash -c "$BENCH_CMD"

    # 3. libmceliece 非f版 1コア
    if check_exec "libmceliece 非f版 (BENCH_LM)" "$BENCH_LM"; then
        run_block "libmceliece 8192128 (非f) 1コア (CPU $CPU1) 周回$round" \
            env RUNS_LM="$RUNS_LM" taskset -c "$CPU1" "$BENCH_LM"
    fi

    # 4. libmceliece f版 1コア
    if check_exec "libmceliece f版 (BENCH_LMF)" "$BENCH_LMF"; then
        run_block "libmceliece 8192128f 1コア (CPU $CPU1) 周回$round" \
            env RUNS_LM="$RUNS_LM" taskset -c "$CPU1" "$BENCH_LMF"
    fi
done

echo | tee -a "$OUT"
echo "結果を保存しました: $OUT" | tee -a "$OUT"