#!/bin/bash
set -euo pipefail

case_dir="/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched"
cd "$case_dir"

stages=(pc_init pc_smoke pc_smoke_offpole pc_medium pc_long)

if [[ $# -gt 1 ]]; then
  echo "Usage: $0 [all|pc_init|pc_smoke|pc_smoke_offpole|pc_medium|pc_long]"
  exit 2
fi

target="${1:-all}"
if [[ "$target" != "all" ]]; then
  stages=("$target")
fi

run_stage() {
  local stage="$1"
  local par_file="${stage}.par"
  local stdout_file="data/${stage}.stdout"
  local log_file="data/${stage}.log"

  echo "== running ${stage} =="
  rm -f "data/${stage}"*.dat "data/${stage}"*.vtu "data/${stage}.log" "data/${stage}.stdout"

  OMP_NUM_THREADS=1 ./amrvac -i "$par_file" > "$stdout_file" 2>&1

  if [[ ! -f "$log_file" ]]; then
    echo "Missing log file for ${stage}: ${log_file}"
    return 1
  fi

  if grep -Eqi '(^|[^[:alpha:]])(nan|inf)([^[:alpha:]]|$)|negative|dtmin|abort|mpistop|fatal|segmentation fault' "$stdout_file" "$log_file"; then
    echo "Detected a hard-failure marker in ${stage}; inspect ${stdout_file} and ${log_file}."
    return 1
  fi

  echo "Last log line for ${stage}:"
  tail -n 1 "$log_file"
  echo
}

for stage in "${stages[@]}"; do
  run_stage "$stage"
done
