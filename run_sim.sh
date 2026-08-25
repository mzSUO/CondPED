#!/bin/bash
TAG=$1; SCRIPT=$2
OUT=/data2/smz/CondPED/inst/formal-simu/output/$TAG
mkdir -p $OUT
# --user scope: runs as the calling user (no sudo, correct PATH/ownership);
# --wait is incompatible with --scope on this systemd; Nice not supported.
systemd-run --user --scope --unit=condped-$TAG \
  -p CPUWeight=10000 -p IOWeight=10000 \
  -p MemoryHigh=70% -p MemoryMax=80% \
  --same-dir --collect \
  bash -lc "cd /data2/smz/CondPED && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/smz/miniconda3/bin/Rscript $SCRIPT 2>&1 | tee $OUT/run.log"
