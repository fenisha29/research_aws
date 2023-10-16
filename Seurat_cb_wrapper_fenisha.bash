#!/bin/bash


declare -a samples
samples=(sample1 sample2 sample3)

for i in ${samples[@]}; do
  nohup Seurat_cb_doubletfinder.R $i > cb_${i}.out 2>&1 &
done
