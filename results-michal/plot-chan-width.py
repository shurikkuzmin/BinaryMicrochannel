#!/usr/bin/python

import matplotlib
matplotlib.use('cairo')

import glob
import numpy as np

from pylab import *

for x in sorted(glob.glob('ca0.05_cw???_1000000.npz')):
    print x
    a = np.load(x)
    width = a['phi'].shape[0] - 2
    mid = width/2

    x1 = np.min(np.where(a['phi'][mid,:] < 0.0))
    x2 = np.max(np.where(a['phi'][mid,:] < 0.0))

    if (x2 - x1 > a['phi'].shape[1] / 2):
        x = x2
    else:
        x = (x2 + x1) / 2

    xdat = np.linspace(0, 1.0, width+1)
    xmax = np.sum(xdat < 0.15)

    plot(xdat[:xmax], a['phi'][:xmax,x], '-', label='%d' % width)

legend()
savefig('chan-width.pdf')


