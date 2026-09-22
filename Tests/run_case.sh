#!/bin/sh
# Runs a ninjaMesher case via mpirun --bind-to none -np "${NINJA_NP:-1}" and asserts exit 0.
# Usage: run_case.sh <ninjaMesher-exe> <caseDir>
set -eu

EXE="$1"
CASE_DIR="$2"

mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" "$CASE_DIR"
