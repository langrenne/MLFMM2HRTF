#!/bin/bash
# Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
# SPDX-License-Identifier: MIT
# Regression tests of the HRTF solvers.
#
#   ./run_test.sh            run the 2 cases and compare them to the references
#   ./run_test.sh --genref   write the references (ref_*.log) from the current code
#
# 2 cases on a sphere of 17 cm (sphere_test.msh, 4010 triangles), default (fast) settings,
# each at 2 kHz (octree levels in BF and HF) and 24 kHz (HF only: reduced CHI depth,
# third band, 33-point near-field quadrature, high orders):
#   RAY   main_HRTF_Ray,  Kreuzer source (default), rigid
#   DIFF  main_HRTF_Diff, monopole source, skin impedance 0.99 / 0.9
# Only the CHECK lines of the logs (verbose 2) are compared: they must be identical.

EXEC_DIR=..
TEST_DIR=$(cd "$(dirname "$0")" && pwd)
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
GENREF=0; [ "$1" = "--genref" ] && GENREF=1
CASES="RAY:Ray DIFF:Diff"
NFAIL=0

echo "Running the solvers..."
for c in $CASES; do
    (cd "$TEST_DIR" && "$EXEC_DIR/main_HRTF_${c##*:}" "config_${c%%:*}.txt" > /dev/null 2>&1
     echo $? > "$TMP/status_${c%%:*}") &
done
wait

for c in $CASES; do
    name=${c%%:*}
    new="$TEST_DIR/new_${name}.log"; ref="$TEST_DIR/ref_${name}.log"
    printf '%-5s ' "$name"
    if [ "$(cat "$TMP/status_$name")" != "0" ] || [ ! -f "$new" ]; then
        echo "ERROR: the solver failed (see $new)"; NFAIL=$((NFAIL+1)); continue
    fi
    if [ $GENREF = 1 ]; then cp "$new" "$ref"; echo "reference written: $ref"; continue; fi
    if [ ! -f "$ref" ]; then echo "no reference: run ./run_test.sh --genref first"; NFAIL=$((NFAIL+1)); continue; fi
    if diff <(grep '^CHECK' "$ref") <(grep '^CHECK' "$new") > "$TMP/diff_$name"; then
        echo "OK ($(grep -c '^CHECK' "$new") CHECK lines)"
    else
        echo "FAIL, differences (< reference, > new):"; cat "$TMP/diff_$name"; NFAIL=$((NFAIL+1))
    fi
done

echo
if [ $NFAIL = 0 ]; then echo "=== all tests OK ==="; else echo "=== $NFAIL test(s) FAILED ==="; exit 1; fi
