#!/usr/bin/env bash
# Every job of tests/benchmark_scale.py, two at a time, niced, each under an
# address-space cap so a runaway dies instead of swapping the production vhosts.
# Rows append to tests/_test_benchmark_scale.jsonl. Run from the repo root, in tmux:
#   tmux new -d -s bench 'bash tests/run_benchmark_scale.sh'
set -u
cd "$(dirname "$0")/.."
PY="${PY:-$(conda info --base 2>/dev/null)/envs/simplicity/bin/python}"
ROOT="${ROOT:-/tmp/test_benchmark_scale}"
OUT="${OUT:-tests/_test_benchmark_scale.jsonl}"
PROCS="${PROCS:-2}"
CAP_KB="${CAP_KB:-3600000}"
mkdir -p "$ROOT"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PY" -I tests/benchmark_scale.py --list | grep -v '^$' |
  xargs -P "$PROCS" -I{} bash -c \
    'ulimit -v '"$CAP_KB"'; nice -n 19 "$0" -I tests/benchmark_scale.py --one {} --root "$1" --out "$2" > "$1/{}.log" 2>&1; echo "$(date +%T) done {} (exit $?)"' \
    "$PY" "$ROOT" "$OUT"
echo "$(date +%T) all done"
rm -rf "$ROOT"
