#!/usr/bin/env python
###############################################################################
#
# A CEH version of interception-corrected potential evaporation
#
###############################################################################

import numpy as np
import sys
import string
import argparse
import netCDF4 as nc
import datetime as dt
import calendar as cal

from . import constants as const
from . import solar

###############################################################################
###############################################################################
#
# Define things
#
###############################################################################
###############################################################################


###############################################################################
# Get day/night LW down
###############################################################################
def get_daynight_lwdown(lwdown, dtr, tair, time_down, time_up, time_max):

    lwdown_day = lwdown_subdaily(lwdown, time_up, time_down, time_max,
                                 tair, dtr)
    lwdown_night = lwdown_subdaily(lwdown, time_down-const.daylen_hour,
                                   time_up, time_max, tair, dtr)

    return lwdown_day, lwdown_night


def lwdown_subdaily(lwdown_mean, t1, t2, time_max, tair, dtr):
    D = const.daylen_hour
    lwdown = lwdown_mean * (1.0 + ((np.sin(2.0*np.pi*(t2 - time_max)/D) -
                                    np.sin(2.0*np.pi*(t1 - time_max)/D)) *
                            dtr * D / (tair*np.pi*(t2-t1))))
    return lwdown


###############################################################################
# Get day/night humidities
###############################################################################
def get_daynight_qair(qmean, tair_mean, psurf, tair_day, tair_night):

    relhum = get_relhum_from_qair(qmean, tair_mean, psurf)

    qair_day = get_qair_from_relhum(relhum, tair_day, psurf)
    qair_night = get_qair_from_relhum(relhum, tair_night, psurf)

    return qair_day, qair_night


def get_qair_from_relhum(relhum, tair, psurf):

    # Get reference temperature (K)
    tr = tref(tair)

    # Get qsat as function of tref (kg kg-1)
    qs = qsat(tr, psurf)

    return (relhum/100.0) * qs * psurf / (psurf + (qs * ((relhum/100.0)-1.0)))


def get_relhum_from_qair(qair, tair, psurf):

    # Get reference temperature (K)
    tr = tref(tair)

    # Get qsat as function of tref (kg kg-1)
    qs = qsat(tr, psurf)

    return 100.0*qair*(1.0-qs)/(qs*(1.0-qair))


###############################################################################
# Get day/night temperatures
###############################################################################
def get_daynight_tair(tair_mean, dtr, time_up, time_down, time_max):

    tair_day = tair_subdaily(tair_mean, time_up, time_down, time_max, dtr)

    tair_night = tair_subdaily(tair_mean, time_down-const.daylen_hour,
                               time_up, time_max, dtr)

    return tair_day, tair_night


def tair_subdaily(tair_mean, t1, t2, time_max, dtr):
    D = const.daylen_hour
    return tair_mean + (dtr/2.0) * (D / (2.0*np.pi*(t2 - t1))) * \
        (np.sin(2.0*np.pi*(t2 - time_max)/D) -
         np.sin(2.0*np.pi*(t1 - time_max)/D))


###############################################################################
# Reference temperature Tr
###############################################################################
def tref(tair):
    # Tr = 1 - Ts/Ta
    tref = 1.0-(const.Ts/tair)
    return tref


###############################################################################
# Vapour pressure at saturation qsat
###############################################################################
def esat(tref):
    sumat = sum([const.a_coeffs[i]*(tref**(i+1)) for i in range(4)])
    esat = const.Ps*np.exp(sumat)
    return esat


###############################################################################
# Specific humidity at saturation qsat
###############################################################################
def qsat(tref, psurf):
    es = esat(tref)
    qsat = const.epsilon * es/(psurf-((1.0-const.epsilon)*es))
    return qsat


###############################################################################
# Gradient of qsat curve
###############################################################################
def Del(tair, psurf):
    tr = tref(tair)
    es = esat(tr)
    qs = qsat(tr, psurf)
    sumat_deriv = sum([const.a_coeffs[i]*(i+1)*(tr**i) for i in range(4)])
    Del = psurf * const.Ts * sumat_deriv * qs / \
        (tair * tair * (psurf - ((1.0 - const.epsilon) * es)))
    return Del

###############################################################################
# Apply threshold
###############################################################################
def apply_threshold(data, thresh, type):
    if type=='max':
        if np.ma.is_masked(data):
            data[np.logical_and(~data.mask,data>thresh)] = thresh
        else:
            data[data>thresh] = thresh
    elif type=='min':
        if np.ma.is_masked(data):
            data[np.logical_and(~data.mask,data<thresh)] = thresh
        else:
            data[data<thresh] = thresh
    else:
        sys.exit("Error: unknown type <%s>"%type)

    return

###############################################################################
# Aerodynamic resistance ra
# From MORECS for short crops (eq 4.36)
###############################################################################
def raero(wind10, h):
    z0 = 0.1 * h
    ra = 6.25 * np.log(10./z0) * np.log(6./z0) / wind10
    return ra


###############################################################################
# canopy resistance (as MORECS)
###############################################################################
def canres_day(lai, rsc, rss):
    f = 0.7
    A = f**lai
    rs = 1./(((1.0-A)/rsc)+(A/rss))
    return rs


def canres_night(lai, rsc, rss):
    rs = 1.0 / ((lai/2500.)+(1.0/rss))
    return rs


###############################################################################
# Albedo
###############################################################################
def get_albedo(lai, a_c, a_s_d, a_s_w, precip,
               awc=None, awc_low=None, awc_high=None):

    alb_soil_wet = np.ones_like(precip)*a_s_w[1]
    alb_soil_dry = np.ones_like(precip)*a_s_d[1]

    if awc is not None:
        alb_soil_wet[np.logical_and(~alb_soil_wet.mask, awc <= awc_low)] \
            = a_s_w[2]
        alb_soil_wet[np.logical_and(~alb_soil_wet.mask, awc >= awc_high)] \
            = a_s_w[0]
        alb_soil_dry[np.logical_and(~alb_soil_dry.mask, awc <= awc_low)] \
            = a_s_d[2]
        alb_soil_dry[np.logical_and(~alb_soil_dry.mask, awc >= awc_high)] \
            = a_s_d[0]

    albedo = np.ma.where(lai > 4,
                         a_c,
                         np.ma.where(precip > 0,
                                     alb_soil_wet+0.25*(a_c-alb_soil_wet)*lai,
                                     alb_soil_dry+0.25*(a_c-alb_soil_dry)*lai))

    return albedo


###############################################################################
# Net radiation
# Rnet = SWnet + LWnet
#   = SWdown-(albedo*SWdown) + LWdown - (emiss*sigma*T^4 + (1-emiss)LWdown)
#   = SWdown*(1-albedo) + emiss*(LWdown-sigma*T^4)
###############################################################################
def net_SW(swdown, albedo):
    SWnet = (swdown*(1.0-albedo))
    return SWnet

def net_LW(lwdown, tair, emiss):
    LWnet = (emiss*(lwdown-(const.sigma*tair*tair*tair*tair)))
    return LWnet

def net_radiation(swdown, lwdown, tair, albedo, emiss):
    A = net_SW(swdown, albedo) + net_LW(lwdown, tair, emiss)
    return A


###############################################################################
# Ground heat flux
###############################################################################
def ground_heat_flux_day(Rnet):
    Gd = 0.2 * Rnet
    return Gd


def ground_heat_flux_night(Gd, t1, t2, P):
    Gn = (P - ((t2-t1)*Gd))/(24-(t2-t1))
    return Gn


###############################################################################
# Air density
###############################################################################
def rhoair(psurf, tair):
    rhoair = psurf/(const.r*tair)
    return rhoair


###############################################################################
# Potential evapotranspiration
###############################################################################
def get_pet(tair, qair, psurf, A, ra, rs, emiss, qdefnotneg=False,
            docorr=True):

    # Get reference temperature
    tr = tref(tair)

    # Get qsat as function of tref
    qs = qsat(tr, psurf)

    # Get Del (as function of tair)
    D = Del(tair, psurf)

    # Get air density
    rhoa = rhoair(psurf, tair)
    cpra = const.cp*rhoa

    # humidity defecit
    qdef = qs-qair
    if qdefnotneg:
        qdef[np.where(qdef < 0.0)] = 0.0

    # Correction for using air temperature to calculate net radiation
    if docorr:
        corr = 4.0 * emiss * const.sigma * tair**3

        # Calculate Penman-Monteith potential evaporation
        # PE =  Del*A + cp*rhoair*(qsat-qair)*(1+b*ra/(rhoair*cp))/ra
        #      ----------------------------------
        #          Del + gamma*(1+rs/ra)*(1+b*ra/(rhoair*cp))
        pet = ((D*A) + (cpra*qdef*(1.0+(corr*ra/cpra))/ra)) / \
            (D + (const.gamma*(1.0+(corr*ra/cpra))*(1.0+(rs/ra))))
    else:
        pet = ((D*A) + (cpra*qdef/ra)) / \
            (D + (const.gamma*(1.0+(rs/ra))))

    return pet


###############################################################################
# Intercepted proportion of rain as for MORECS
###############################################################################
def interc_rf(ppt, lai, enhance):

    p = 1.0 - 0.5**lai
    interc = p*ppt

    if np.ma.isMA(interc):
        interc[np.logical_and(interc > (0.2*lai),~interc.mask)] = 0.2*lai
    else:
        interc[interc > (0.2*lai)] = 0.2*lai

    interc *= enhance

    if np.ma.isMA(interc):
        interc[np.logical_and(interc > ppt, ~interc.mask)] = ppt[np.logical_and(interc > ppt, ~interc.mask)]
    else:
        interc[interc > ppt] = ppt[interc > ppt]

    return interc


###############################################################################
# Interception correction to PE
###############################################################################
def interc_corr(pet_t, pet_i, interc, D):

    # Time to evaporate all intercepted water (hours)
    t_to_dry = interc/pet_i * D

    pet_i_corr = interc
    pet_t_corr = (D - t_to_dry)*pet_t/D

    pet_i_corr[t_to_dry >= D] = pet_i[t_to_dry >= D]
    pet_t_corr[t_to_dry >= D] = 0.0

    # mm
    return pet_t_corr, pet_i_corr


###############################################################################
# CO2 function from RC.
# function to calculate surface resistance using intitial MORECS values and
# atmospheric CO2 concentration: varying through year
###############################################################################
def rs_fn_Rudd(co2, month, co2_baseline, lai, rsc):
    co2_delta = np.where(co2-co2_baseline < 1000., co2 - co2_baseline, 1000.0)

    # decrease in stomatal and crop conductance with increased CO2
    # (Kruijt, 2008. Journal of Hydrology)
    g_sc_dec = -0.093 * co2_delta

    # r_sc multiplier
    r_sc_mult = 1./(1.+(g_sc_dec/100.))

    # future r_sc
    F_rsc = r_sc_mult * rsc

    return F_rsc
