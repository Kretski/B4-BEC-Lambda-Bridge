#!/bin/bash
export OMP_NUM_THREADS=1
for bg in K3 K4; do
 python3 t3fund_hs2_c4.py --bg $bg --coupled $1 --N 16384 --dt 0.025 --out .
 python3 t3fund_hs2_c4.py --bg $bg --coupled $1 --N 32768 --dt 0.0125 --out .
done
echo LANE_DONE
