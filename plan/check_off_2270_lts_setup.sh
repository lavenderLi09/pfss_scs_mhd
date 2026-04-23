#!/bin/zsh
set -euo pipefail

root="/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic"
plan_file="$root/plan/2026-04-02-off_2270-lts-first-amrvac32.md"
case_dir="$root/off_2270_lts_32_linde"

fail() {
  printf 'FAIL: %s\n' "$1" >&2
  exit 1
}

[[ -f "$plan_file" ]] || fail "missing plan file: $plan_file"
grep -Fq "This case does not use a theta-phi AMR window; refinement is not window-restricted." "$plan_file" || \
  fail "plan file missing no-theta-phi-window clarification"

[[ -d "$case_dir" ]] || fail "missing experiment case directory: $case_dir"
[[ -f "$case_dir/amrvac.par" ]] || fail "missing experiment amrvac.par"
[[ -f "$case_dir/mod_usr.t" ]] || fail "missing experiment mod_usr.t"
[[ -f "$case_dir/makefile" ]] || fail "missing experiment makefile"

grep -Fq "time_stepper='threestep'" "$case_dir/amrvac.par" || fail "experiment must keep threestep"
grep -Fq "local_timestep=.true." "$case_dir/amrvac.par" || fail "experiment must enable local_timestep"
grep -Fq "typecourant='maxsum'" "$case_dir/amrvac.par" || fail "experiment must keep typecourant maxsum"
grep -Fq "courantpar=0.5d0" "$case_dir/amrvac.par" || fail "experiment must keep courantpar 0.5d0"
grep -Fq "mhd_gravity=.true." "$case_dir/amrvac.par" || fail "experiment must keep gravity enabled"
grep -Fq "typedivbfix='linde'" "$case_dir/amrvac.par" || fail "experiment must switch divB fix to linde"
if grep -Fq "typedivbfix='glm'" "$case_dir/amrvac.par"; then
  fail "experiment amrvac.par still contains glm divB fix"
fi

if grep -n "inwin" "$case_dir/mod_usr.t" >/dev/null; then
  fail "experiment mod_usr.t should not keep inwin/theta-phi window logic"
fi

grep -Fq "(absj/absb)> joverb_cr" "$case_dir/mod_usr.t" || \
  fail "experiment mod_usr.t must keep the primary global J/B AMR threshold"
grep -Fq "(absj/absb)>(joverb_cr-0.3/dxmin)" "$case_dir/mod_usr.t" || \
  fail "experiment mod_usr.t must keep the secondary global J/B AMR threshold"

printf 'PASS: plan file and experiment case are configured as requested\n'
