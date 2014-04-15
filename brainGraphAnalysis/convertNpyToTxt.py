#!/usr/bin/python

"""
Converts .npy data file with lcc to a txt file with index followed by 
xyz location

1: raw data file directory
2: processed data file directory
3: file prefix
"""

import sys
import struct
import numpy as np
import scipy.sparse as sp
import scipy.io
import zindex

rawDataDir = sys.argv[1]
procDataDir = sys.argv[2]
filePre = sys.argv[3]

inputfile = ''.join((rawDataDir, filePre, '_big_lcc.npy'))
outputfile = ''.join((procDataDir, filePre, '_big_lcc.txt'))

#load connected components from file
mat = np.load(inputfile).item().toarray()

#keep only the largest connected component
mat[mat != 1] = 0

#creat sparse matrix and pull out rows
mat = sp.lil_matrix(mat)
mat_rows = np.array(mat.rows.item())
mat_rows.size

#convert rows into xyz coordinates and store them in one matrix
mat_loc = np.zeros((mat_rows.shape[0], 4), dtype=np.int64)
index = 0
for i in mat_rows:
	mat_loc[index, ] = np.append(i, zindex.MortonXYZ(i))
	index += 1;

#output data 
np.savetxt(outputfile, mat_loc, fmt="%9d", delimiter=" ")
