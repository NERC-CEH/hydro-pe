#!/usr/bin/env python

import netCDF4 as nc
import numpy as np

from . import constants
from . import utils

################################################################################
################################################################################
#
#
#
# ELR 03-08-2020
#
################################################################################
################################################################################

def sun_times(daynumber, year, lat, lon):

    delta = solar_declination(daynumber, year)

    # convert lat and lon to radians
    phi = lat * np.pi / 180.0
    lambdal = lon * np.pi /180.0

    tanp_tand = np.tan(phi) * np.tan(delta)

    # Perpetual day if tantan<-1.0, so sun is apparent from midnight to midnight
    # Perpetual night if tantan>1.0, so sun 'rises' and 'sets' at exactly
    # midday, resulting in a day length of zero
    omega_s = np.where(tanp_tand<-1.0, 0.0, 
                       np.where(tanp_tand>1, np.pi,
                                np.arccos(-1.0*tanp_tand)))

    time_up = (1.0 - (omega_s + lambdal)/np.pi) * constants.daylen_hour / 2.0
    time_down = (1.0 + (omega_s - lambdal)/np.pi) * constants.daylen_hour / 2.0

    return time_up, time_down

def time_max_temperature(time_up, time_down, tmax_offset):
    # From HCTN96 (Williams and Clark, 2014), but ultimately from IMOGEN code
    return 0.5 * (time_up+time_down) + tmax_offset * (time_down - time_up)

def solar_declination(daynumber, year):
    # From metprod. Equivalent to MORECS
    dy_in_yr = utils.days_in_year(year)

    return 0.409 * np.sin((2.0*np.pi*daynumber/dy_in_yr) - 1.39)


def solar_declination_2(daynumber, year):
    # From the inclined plane calculations in metprod

    dy_in_yr = utils.days_in_year(year)

    gammad = 2.0 * np.pi * (daynumber-1.0) / dy_in_yr 

    dec = 0.006919 - 0.39912 * np.cos(gammad) + 0.070257 * np.sin(gammad) \
         -0.006758 * np.cos(2.0 * gammad) + 0.000907 * np.sin(2.0 * gammad) \
         -0.002697 * np.cos(3 * gammad) + 0.00146 * np.sin(3.0 * gammad)


    return 0.409 * np.sin((2.0*np.pi*daynumber/dy_in_yr) - 1.39)

