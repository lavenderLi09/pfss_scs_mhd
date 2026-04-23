#!/usr/bin/env bash
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$CASE_DIR"

mkdir -p logs data figs analysis

PAR_FILE="${1:-amrvac.par}"
TAG="${2:-run_$(date +%Y%m%d_%H%M%S)}"

OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" ./amrvac -i "$PAR_FILE" \
  > "logs/${TAG}.out" \
  2> "logs/${TAG}.err"
