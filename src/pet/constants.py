#!/usr/bin/env python

import netCDF4 as nc
import numpy as np
from collections import namedtuple

################################################################################
################################################################################
#
#
#
# ELR 29-07-2020
#
################################################################################
################################################################################

# specific heat of air
cp = 1010.0

# gas constant of dry air (J kg-1 K-1)
r = 287.05

# latent heat (J kg-1)
l = 2.5e6

# gamma: psychometric constant for specific humidity calculations (K-1)
gamma = 0.0004

# MORECS stomatal resistance for grass (m s-1)
rsc_day = [80.,80.,60.,50.,40.,60.,60.,70.,70.,70.,80.,80.]
rsc_night = 2500.

# MORECS wet bare soil surface resistance m s-1
rss = 100.0

# MORECS grass LAI
lai = [2.,2.,3.,4.,5.,5.,5.,5.,4.,3.,2.5,2.0]

# albedos
# Dry bare soil (high-medium-low AWC) albedo MORECS
albedo_s_dry = [0.1, 0.2, 0.3]
# Wet bare soil (high-medium-low-low AWC) albedo MORECS
albedo_s_wet = [0.05, 0.1, 0.15]
# Full crop albedo MORECS grass
albedo_c = 0.25

# MORECS emissivity
emiss = 0.95

# Use MORECS canopy height
canht = 0.15

# ratio mol weight water vapour/dry air
epsilon = 0.622

# steam point temperature
Ts = 373.15

# steam point pressure
Ps = 101325.0

# Stefan-Boltzmann constant (W m-2 K-4)
sigma = 5.67e-8

# coeffs for fit of saturated specific humidity (need citation)
a_coeffs = [13.3185, -1.9760, -0.6445, -0.1299]

# ground heat storage (W hr m-2) From MORECS v2 (Hydrol Mem 45)
P = [-137., -75.,  30.,  167.,  236.,  252.,
      213.,  69., -85., -206., -256., -206.]

# MORECS interception enhancement factor
enhance = [1.0, 1.0, 1.2, 1.4, 1.6, 2.0, 2.0, 2.0, 1.8, 1.4, 1.2, 1.0]

# time constants
daylen_sec = 86400.0
daylen_hour = 24.0
hourlen_sec = 3600.0

# Solar constants

# Longitude of perihelion
# https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
lonperi = 102.94719 * np.pi / 180.

