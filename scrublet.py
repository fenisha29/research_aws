# -*- coding: utf-8 -*-
"""Untitled1.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1ojurFCpG9Z65aiyO11mXZmZvK-sWxTW_
"""

#!/usr/bin/env python3

import sys
sys.path.append('/usr/local/lib/python3.6/dist-packages/')

import scrublet as scr
import scipy.io
import sys

#patient=sys.argv[1]
patient = "RU1065C_MET_LI"
#doublet_rate=float(sys.argv[2])
doublet_rate= 0.0644

print('Patient:', patient, 'Doublet_rate:', doublet_rate)
print(sys.argv)

counts_matrix = scipy.io.mmread(''.join(['data/', patient, '/matrix_', patient, '_raw.mtx'])).T.tocsc()
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=doublet_rate)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=1, min_cells=1, min_gene_variability_pctl=50, n_prin_comps=20)
print(doublet_scores, predicted_doublets)

outFI = open(''.join(['data/', patient, '/doublets_', patient, '_raw.txt']),'w')
outFI.write('predicted_doublets\tdoublet_scores\n')
for i in range(len(predicted_doublets)):
  outFI.write(str(predicted_doublets[i]).upper()+'\t'+str(doublet_scores[i])+'\n')

outFI.close()
