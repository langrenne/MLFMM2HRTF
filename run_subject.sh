#!/bin/bash
# Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
# SPDX-License-Identifier: MIT
# run_subject.sh : HRTFs of one subject, both ears, from the case files to the SOFA file.
# Usage: ./run_subject.sh [options] CASE_L CASE_R      (./run_subject.sh -h for the options)
set -u

usage() {
cat <<'EOF'
Usage: run_subject.sh [options] CASE_L CASE_R

Computes the HRTFs of one subject and writes the SOFA file. CASE_L and CASE_R
are the case files of the left and right ears (config.md), in the same
directory; everything runs in that directory.

  1. solver, left ear then right ear: main_HRTF_Ray if the source (value 7) is
     a triangle number, main_HRTF_Diff if it is a point x y z. Skipped for an
     ear whose result file <stem>.txt already exists (--force to recompute);
  2. if the frequency line is a list (mixed grid): interpolation to 512
     frequencies (python/interpolate_hrtf.py);
  3. HRIRs with diffuse-field equalisation (python/hrtf2hrir.py);
  4. SOFA file <subject>.sofa (python/npz2sofa.py; metadata from
     sofa_meta.txt in the run directory if present).

Options:
  -s, --subject NAME  subject name (default: result stem of CASE_L without a
                      final _L)
  -w, --weights FILE  quadrature weights of the directions (default:
                      weights_<grid>_<N>.txt if the external points file is
                      pext_<grid>_<N>.txt (Lebedev or Fibonacci) and that file
                      exists; else uniform)
      --wav           also write one WAV file per direction and ear, with an
                      index, in wav_<NAME>/
  -t, --threads N     OpenMP threads of the solver (default: OMP_NUM_THREADS
                      if set, else the number of cores)
      --post-only     skip step 1 (the result files must exist)
      --force         recompute an ear even if its result file exists
      --wait          before each solver run, wait until no other
                      main_HRTF_Ray / main_HRTF_Diff runs on the machine
      --python PY     Python interpreter (default: $PYTHON, else python3)
  -h, --help          this help

Progress: run_<NAME>.log ; solver screen output: <stem>.out
Example: ./run_subject.sh --wav config_S01_L.txt config_S01_R.txt
EOF
}

REPO=$(cd "$(dirname "$0")" && pwd)
SUBJECT= ; WEIGHTS= ; WAV=0 ; THREADS=${OMP_NUM_THREADS:-$(nproc)}
POST_ONLY=0 ; FORCE=0 ; WAIT=0 ; PY=${PYTHON:-python3}
CASES=()
while [ $# -gt 0 ]; do
  case $1 in
    -s|--subject) SUBJECT=$2; shift 2 ;;
    -w|--weights) WEIGHTS=$2; shift 2 ;;
    --wav)        WAV=1; shift ;;
    -t|--threads) THREADS=$2; shift 2 ;;
    --post-only)  POST_ONLY=1; shift ;;
    --force)      FORCE=1; shift ;;
    --wait)       WAIT=1; shift ;;
    --python)     PY=$2; shift 2 ;;
    -h|--help)    usage; exit 0 ;;
    -*)           echo "run_subject.sh: unknown option $1 (-h for help)" >&2; exit 2 ;;
    *)            CASES+=("$1"); shift ;;
  esac
done
[ ${#CASES[@]} -eq 2 ] || { usage >&2; exit 2; }

die() { echo "run_subject.sh: $*" >&2; [ -n "${LOG:-}" ] && echo "$(date '+%F %T')  ERROR $*" >> "$LOG"; exit 1; }

for f in "${CASES[@]}"; do [ -f "$f" ] || die "case file not found: $f"; done
DIR_L=$(cd "$(dirname "${CASES[0]}")" && pwd); DIR_R=$(cd "$(dirname "${CASES[1]}")" && pwd)
[ "$DIR_L" = "$DIR_R" ] || die "both case files must be in the same directory"
CASE_L=$(basename "${CASES[0]}"); CASE_R=$(basename "${CASES[1]}")
cd "$DIR_L" || exit 1

# value n of a case file (value lines = non-empty lines not starting with '!')
value() { tr -d '\r' < "$1" | awk '!/^[[:space:]]*!/ && NF' | sed -n "${2}p" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//'; }

STEM_L=$(value "$CASE_L" 6); STEM_R=$(value "$CASE_R" 6)
[ -n "$STEM_L" ] && [ -n "$STEM_R" ] || die "cannot read the result stems (value 6) of the case files"
[ "$STEM_L" != "$STEM_R" ] || die "both ears write the same result file $STEM_L.txt"
[ -n "$SUBJECT" ] || SUBJECT=${STEM_L%_L}
LOG=run_$SUBJECT.log
log() { echo "$(date '+%F %T')  $*" | tee -a "$LOG"; }

solver_of() {   # Ray: source = one integer ; Diff: source = x y z
  local n; n=$(value "$1" 7 | wc -w)
  if [ "$n" -eq 1 ]; then echo main_HRTF_Ray; elif [ "$n" -ge 3 ]; then echo main_HRTF_Diff
  else die "$1: cannot read the source (value 7)"; fi
}

export OPENBLAS_NUM_THREADS=1      # parallelism is on the frequencies (OpenMP)
export OMP_NUM_THREADS=$THREADS

log "subject $SUBJECT : $CASE_L -> $STEM_L, $CASE_R -> $STEM_R"

# ---- 1. solver, one ear after the other
if [ $POST_ONLY -eq 0 ]; then
  for c in "$CASE_L" "$CASE_R"; do
    stem=$(value "$c" 6); prog=$(solver_of "$c")
    [ -x "$REPO/$prog" ] || die "$REPO/$prog not found: build it first (make)"
    if [ -s "$stem.txt" ] && [ $FORCE -eq 0 ]; then
      log "$stem.txt exists, solver skipped (--force to recompute)"; continue
    fi
    if [ $WAIT -eq 1 ]; then
      while pgrep -x 'main_HRTF_Ray|main_HRTF_Diff' > /dev/null; do sleep 60; done
    fi
    log "$prog $c ($OMP_NUM_THREADS threads)"
    "$REPO/$prog" "$c" > "$stem.out" 2>&1 || die "$prog $c failed (see $stem.out and the log file of the case)"
    [ -s "$stem.txt" ] || die "$prog $c: no result file $stem.txt (see $stem.out)"
    log "done: $stem.txt"
  done
fi
for s in "$STEM_L" "$STEM_R"; do [ -s "$s.txt" ] || die "result file $s.txt missing"; done

# ---- 2. mixed grid: interpolation to 512 frequencies
POST_L=$CASE_L; RES_R=$STEM_R
if [ "$(value "$CASE_L" 3 | tr ',' ' ' | wc -w)" -gt 1 ]; then
  for c in "$CASE_L" "$CASE_R"; do
    "$PY" "$REPO/python/interpolate_hrtf.py" "$c" >> "$LOG" 2>&1 || die "interpolation of $c failed (see $LOG)"
  done
  POST_L=config_${STEM_L}_512.txt; RES_R=${STEM_R}_512
  log "interpolated to 512 frequencies: ${STEM_L}_512.txt, ${STEM_R}_512.txt"
fi

# ---- 3. HRIRs
if [ -z "$WEIGHTS" ]; then
  pext=$(value "$CASE_L" 5)
  case $pext in
    *pext_*.txt) b=$(basename "$pext")                 # pext_lebedev_N.txt -> weights_lebedev_N.txt
                 w=$(dirname "$pext")/weights_${b#pext_}
                 [ -f "$w" ] && WEIGHTS=$w ;;
  esac
fi
ARGS=(--result-R "$RES_R" -o "hrir_$SUBJECT.npz")
[ -n "$WEIGHTS" ] && ARGS+=(--weights "$WEIGHTS")
[ $WAV -eq 1 ] && ARGS+=(--wav-dir "wav_$SUBJECT")
log "HRIR: hrtf2hrir.py $POST_L ${ARGS[*]}"
"$PY" "$REPO/python/hrtf2hrir.py" "$POST_L" "${ARGS[@]}" >> "$LOG" 2>&1 || die "hrtf2hrir.py failed (see $LOG)"

# ---- 4. SOFA
"$PY" "$REPO/python/npz2sofa.py" "hrir_$SUBJECT.npz" "$SUBJECT.sofa" --subject "$SUBJECT" >> "$LOG" 2>&1 \
  || die "npz2sofa.py failed (see $LOG)"
log "written: $SUBJECT.sofa$([ $WAV -eq 1 ] && echo ", wav_$SUBJECT/")"
