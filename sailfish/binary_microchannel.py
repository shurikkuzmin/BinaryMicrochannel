#!/usr/bin/python

import numpy as np
from sailfish import lbm
from sailfish import geo

import optparse
from optparse import OptionGroup, OptionParser, OptionValueError

# TODO: [(736, 51), (1471, 100), (2206, 149), (2941, 198)]

class GeoFE(geo.LBMGeo2D):

    def define_nodes(self):
        hy, hx = np.mgrid[0:self.lat_ny, 0:self.lat_nx]
        wall_map = np.logical_or(hy == 0, hy == self.lat_ny-1)
        self.set_geo(wall_map, self.NODE_WALL)
 
    def init_fields(self):
        hy, hx = np.mgrid[0:self.lat_ny, 0:self.lat_nx]
        width = 10

        self.sim.phi[:] = 1.0
        self.sim.phi[np.logical_and(
                       np.logical_and(hx >= (self.lat_nx-1)/4.0,
                                      hx <= 3*(self.lat_nx-1)/4.0),
                       np.logical_and(hy >= width, hy <= self.lat_ny-width-1))] = -1.0
        self.sim.rho[:] = 1.0

class FESim(lbm.BinaryFluidFreeEnergy):
    filename = 'fe_seperation_2d'

    def __init__(self, geo_class, defaults={}):
        settings = {'verbose': True, 'lat_nx': 736,
                    'lat_ny': 51, 'grid': 'D2Q9',
                    'kappa': 0.04, 'Gamma': 1.0, 'A': 0.04,
                    'scr_scale': 2,
                    'tau_a': 2.5, 'tau_b': 0.7, 'tau_phi': 1.0,
                    'periodic_x': True, 'periodic_y': True}
        settings.update(defaults)

        lbm.BinaryFluidFreeEnergy.__init__(self, geo_class, options=[], defaults=settings)

        self.add_body_force((6.0e-6, 0.0), grid=0, accel=False)

        # Use the fluid velocity in the relaxation of the order parameter field,
        # and the molecular velocity in the relaxation of the density field.
        self.use_force_for_eq(None, 0)
        self.use_force_for_eq(0, 1)


if __name__ == '__main__':
    sim = FESim(GeoFE)
    sim.run()
