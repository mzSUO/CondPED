#!/bin/bash
TAG=$1; SCRIPT=$2
mkdir -p output/$TAG
sudo systemd-run --scope --unit=condped-$TAG \
  -p CPUWeight=10000 -p IOWeight=10000 -p Nice=-10 \
  -p MemoryHigh=70% -p MemoryMax=80% \
  --same-dir --wait --collect \
  bash -lc "cd /data2/smz/CondPED && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript $SCRIPT 2>&1 | tee output/$TAG/run.log"
