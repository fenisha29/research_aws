#!/bin/bash


declare -a samples
samples=(RU1066_PRI_LU RU1080C_MET_KI RU1108a_REC_LU_b RU1108a_REC_LU_bf RU1108a_REC_LU_rpmi RU1124A_MET_LN RU1144_MET_LN RU1144_REC_LU RU1145_PRI_LU RU1152_MET_LN)

for i in ${samples[@]}; do
  nohup ./Seurat_cb_doubletfinder.R $i > cb_${i}.out 2>&1 &
done
