#!/bin/bash
# Registered C4 runs for one designated background (case label[/mode]); two lanes can run in parallel. Resumable (checkpoints). Usage: run_c4.sh MANIFEST OUTDIR CASE COUPLED(1/0)
export OMP_NUM_THREADS=1
python3 t3fund_proca_c4.py --manifest $1 --out $2 --case $3 --coupled $4 --N 16384 --dt 0.025
python3 t3fund_proca_c4.py --manifest $1 --out $2 --case $3 --coupled $4 --N 32768 --dt 0.0125
