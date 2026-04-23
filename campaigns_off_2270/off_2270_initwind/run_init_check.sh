#!/bin/bash
set -euo pipefail

case_dir="/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_initwind"
cd "$case_dir"

export AMRVAC_DIR="/Users/zhaoyan/Documents/codes/amrvac3.2"

check_no_active_conversion_calls() {
  if rg -n "^[[:space:]]*call[[:space:]]+mhd_to_(primitive|conserved)\\(" mod_usr.t; then
    echo "Found active mhd_to_primitive/mhd_to_conserved calls in mod_usr.t; aborting."
    return 1
  fi
}

check_no_active_conversion_calls

echo "== build: make clean -> setup.pl -> make =="
if [[ -f makefile ]]; then
  make clean
else
  echo "makefile not found (first build); skip make clean this time"
fi
"$AMRVAC_DIR/setup.pl" -d=3 -arch=debug
make

mkdir -p data analysis

echo "== run init_check.par =="
OMP_NUM_THREADS=1 ./amrvac -i init_check.par > data/init_check.stdout 2>&1

if [[ ! -f data/init_check.log ]]; then
  echo "Missing data/init_check.log"
  exit 1
fi

if rg -n -i "(^|[^[:alpha:]])(nan|inf)([^[:alpha:]]|$)|negative|dtmin|abort|mpistop|fatal|segmentation fault" data/init_check.stdout data/init_check.log; then
  echo "Detected a hard-failure marker in init_check run"
  exit 1
fi

check_no_active_conversion_calls

echo "== done =="
ls -1 data/init_check* | sed 's#^#  #'
