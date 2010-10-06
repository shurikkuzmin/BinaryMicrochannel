#!/usr/bin/python

import sys
import numpy as np

a = np.load(sys.argv[1])

bnd = np.argwhere(a['phi'] < 0)
midx = (np.min(bnd[:,1]) + np.max(bnd[:,1])) / 2
phi = a['phi'][:,midx]

for i in np.argwhere(phi[1:] * phi[:-1] < 0):
    i = i[0]
    w = abs(phi[i]) + abs(phi[i+1])
    pos = (i * abs(phi[i+1]) + (i+1) * abs(phi[i])) / w 
    print pos
