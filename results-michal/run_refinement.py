#!/usr/bin/python

Ca = 0.05
gamma = 0.03771

for y in range(100, 250, 25):
    ny = y + 2
    nx = y * 15

    iwidth = int(0.06 * y)
    force = Ca * 8.0 * gamma / y**2 

    print "./bmc.py --lat_nx=%d --lat_ny=%d --iwidth=%d --force=%e --output=ca0.05_cw%d_ --every=100000 --max_iters=1000000" % (nx, ny, iwidth, force, y)
