#!/usr/bin/env python
import sys
import numpy as np
import os
import re

# Read the file into a dictionary, labelled by event->timeseries
sep = ';'
with open(sys.argv[1],"r") as fnin:
    line = fnin.readline().strip()
    colnames = re.sub(r'[^\x00-\x7F]','', line).split(sep)
    ensembles = { colname:[] for colname in colnames}
    while line:
        line = fnin.readline()
        linesplit = line.strip().split(sep)
        for i in range(len(linesplit)):
            try:
                ensembles[colnames[i]].append(float(linesplit[i].replace(',','.')))
            except:
                ensembles[colnames[i]].append(linesplit[i])

first_perturbed = 0
last_perturbed = 100
standard = 113

label2nr = {}
rownames = ensembles['Scenario']
for i in range(len(rownames)):
    label2nr[rownames[i]] = i
selected = ['Sum%d'%dd for dd in range (first_perturbed,last_perturbed+1)]    # selected labels for realizations
selnr = [label2nr[rowname] for rowname in selected]

irow_standard = label2nr['Sum%d'%standard]
for icol in range(1,len(colnames)): 
    standard=ensembles[colnames[icol]][irow_standard]
    sample=np.array([ensembles[colnames[icol]][i] for i in selnr])
    print ("\n\n%s : %f\n"%(colnames[icol],standard))
    print (sample)
    # call beta-fit procedure on this sample 
    # print ens


sys.exit()


