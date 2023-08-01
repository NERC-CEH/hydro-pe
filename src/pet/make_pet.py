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

from .pet import *
from . import constants as const
from . import utils

###############################################################################
# Parse the input
###############################################################################
def parse_input():
    parser = argparse.ArgumentParser(description="Calculate PETI")

    # optional
    parser.add_argument("-t", "--tairvar",
                        help="Air temperature variable",
                        required=False, default="tas")
    parser.add_argument("-d", "--dtrvar",
                        help="Daily air temperature range variable",
                        required=False, default="dtr")
    parser.add_argument("-q", "--qairvar",
                        help="Specific humidity variable",
                        required=False, default="huss")
    parser.add_argument("-w", "--windvar",
                        help="Wind speed (10m) variable",
                        required=False, default="sfcWind")
    parser.add_argument("-W", "--windthresh",
                        help="Minimum threshold for wind speed",
                        required=False, default=None, type=float)
    parser.add_argument("-s", "--swvar",
                        help="SWdown variable",
                        required=False, default="rsds")
    parser.add_argument("-l", "--lwvar",
                        help="LWdown variable",
                        required=False, default="rlds")
    parser.add_argument("-p", "--psurfvar",
                        help="Surface air pressure variable",
                        required=False, default="psurf")
    parser.add_argument("-P", "--precipvar",
                        help="Precipitation variable",
                        required=False, default="precip")
    parser.add_argument("--latvar",
                        help="Latitude variable",
                        required=False, default="lat")
    parser.add_argument("--lonvar",
                        help="Longitude variable",
                        required=False, default="lon")
    parser.add_argument("-T", "--toffset",
                        help="Temperature offset",
                        required=False, type=float, default=0.0)
    parser.add_argument("-S", "--pscale",
                        help="Air pressure scale",
                        required=False, type=float, default=1.0)
    parser.add_argument("-L", "--precipscale",
                        help="Precipitation scale",
                        required=False, type=float, default=1.0)
    parser.add_argument("-Z", "--zlib",
                        help="Create compressed netCDF files",
                        required=False, action="store_true")
    parser.add_argument("-z", "--zeroneg",
                        help="Negative values set to zero",
                        required=False, action="store_true")
    parser.add_argument("-v", "--verbose",
                        help="Verbose",
                        required=False, action="store_true")
    parser.add_argument("-V", "--version",
                        help="Version for the files",
                        required=False, default="v1.0")
    parser.add_argument("-Q", "--qdefnotneg",
                        help="Don't allow humidity defecit to be negative",
                        required=False, action="store_true")
    parser.add_argument("-f", "--fullout",
                        help="Output day and night (and times) "
                             "as well as total(s)",
                        required=False, action="store_true")
    parser.add_argument("-G", "--precipscalegridfilevarmn",
                        help="Precip scale grid file, variable and month variable", nargs=3,
                        default=[None, None, None])
    parser.add_argument("-C", "--co2filevaryr",
                        help="CO2 file, variable and baseline year", nargs=3,
                        default=[None, None, None])
    parser.add_argument("-N", "--netradiation",
                        help="Input net radiation", action="store_true")
    parser.add_argument("-c", "--docorr",
                        help="Apply correction for net radiation calculated "
                             "with air temperature instead of surface "
                             "temperature in LWup", action="store_true")
    parser.add_argument("--timevar",
                        help="Time variable",
                        required=False, default="time")
    parser.add_argument("--groundheatflux",
                        help="How to calculate ground heat flux. "
                             "Current options: morecs or zero.",
                        required=False, default="morecs")
    parser.add_argument("--gridmapvar",
                        help="Grid mapping variable",
                        required=False, default=None)
    parser.add_argument("-u", "--user",
                        help="Name of user who created this data",
                        required=False, default=None)
    parser.add_argument("-e", "--email",
                        help="Email of user who created this data",
                        required=False, default=None)
    parser.add_argument("-i", "--interception",
                        help="Do interception correction",
                        required=False, action="store_true")
    parser.add_argument("-D", "--daynight",
                        help="Do separate day and night calculations",
                        required=False, action="store_true")

    # positional
    parser.add_argument("datafiletmplt",
                        help="Template of file containing "
                             "input ($varn)")
    parser.add_argument("outfilename",
                        help="Output file name")

    args = parser.parse_args()

    if args.netradiation:
        varns = {
            "tair": args.tairvar,
            "qair": args.qairvar,
            "wind10": args.windvar,
            "swnet": args.swvar,
            "lwnet": args.lwvar,
            "psurf": args.psurfvar
            }

        # Check if the correction factor is required
        docorr = args.docorr
    else:

        varns = {
            "tair": args.tairvar,
            "qair": args.qairvar,
            "wind10": args.windvar,
            "swdown": args.swvar,
            "lwdown": args.lwvar,
            "psurf": args.psurfvar
            }

        # Always do the correction if we're calculating the net radiation
        # within the code
        docorr = True

    if args.interception:
        varns["precip"] = args.precipvar

    if args.daynight:
        varns["dtr"] = args.dtrvar

    if args.co2filevaryr[2] is None:
        co2yr = None
    else:
        co2yr = int(args.co2filevaryr[2])

    if args.groundheatflux not in ["morecs", "zero"]:
        sys.exit("Error: unknown groundheatflux type: <%s>"%args.groundheatflux)

    if args.daynight and args.netradiation:
        sys.exit("Error: Can't do daily disaggregation with net radiation")

    return string.Template(args.datafiletmplt), \
        args.outfilename, \
        varns, \
        args.interception, \
        args.daynight, \
        args.zlib, \
        args.zeroneg, \
        args.version, \
        args.verbose, \
        args.qdefnotneg, \
        args.fullout, \
        args.co2filevaryr[0], \
        args.co2filevaryr[1], \
        co2yr, \
        args.precipscalegridfilevarmn[0], \
        args.precipscalegridfilevarmn[1], \
        args.precipscalegridfilevarmn[2], \
        args.netradiation, \
        args.docorr, \
        args.groundheatflux, \
        args.toffset, \
        args.pscale, \
        args.precipscale, \
        args.windthresh, \
        args.timevar, \
        args.latvar, \
        args.lonvar, \
        args.gridmapvar, \
        args.user, \
        args.email


###############################################################################
###############################################################################
#
# Start the main routine
#
###############################################################################
###############################################################################
def main():

    ########################################################################
    # Parse the input
    datafiletmplt, outfilename, varns, interception, daynight, \
        zlib, zneg, version, verbose, qdefnotneg, fullout, \
        co2file, co2var, co2yr, precipscalegridfile, precipscalegridvar, \
        precipscalegridmn, netradiation, docorr, groundheatflux, \
        toffset, pscale, precipscale, windthresh, timevar, latvar, lonvar, \
        gridmapvar, user, email = parse_input()

    ########################################################
    # Read CO2
    if co2file is None:
        co2 = None
    else:
        if verbose:
            print("Read CO2 file")
        co2, co2yrs, co2_baseline = utils.read_co2_file(co2file, co2var,
                                                        co2yr,
                                                        timevar=timevar)

    ########################################################
    # Read precip scale grid
    if precipscalegridfile is None:
        precipscalegrid = None
    else:
        if verbose:
            print("Read precip scale grid file")
        precipscalegrid, precipmonths = utils.read_precipscalegrid_file(
                            precipscalegridfile, precipscalegridvar,
                            monthvar=precipscalegridmn)

    ########################################################
    # Loop over required variables
    if verbose:
        print("Read variables")
    data, dims, dimvardata = utils.read_input_data(datafiletmplt, varns,
                                                   verbose=verbose,
                                                   timevar=timevar)

    # get shape
    ne, nt, ny, nx = utils.get_shape(data["tair"])

    # Rearrange data if necessary
    if ne > 0:
        if verbose:
            print("Reorder data to cope with ensemble members")
        for varn in data:
            if len(data[varn].shape) == 4:
                data[varn] = data[varn].swapaxes(0, 1)

    # get attributes
    if interception:
        outvarn = "peti"
    else:
        outvarn = "pet"
    outvarns = {
        outvarn: utils.get_var_attrs(datafiletmplt.substitute(varn=varns["tair"]),
        varns["tair"])
    }

    # Check fill value
    if '_FillValue' in outvarns[outvarn]:
        mv = outvarns[outvarn]['_FillValue']
    else:
        mv = -1e20

    # Get the co2 index if necessary
    if co2file is not None:
        co2indx = utils.get_co2_indx(co2yrs, dimvardata["year"])

    # Apply wind threshold
    # threshold is 0.1 on 2m wind, for consistency with chess
    if windthresh is not None:
        if verbose:
            print("Applying minimum threshold on wind speed %f" % windthresh)
        apply_threshold(data["wind10"], windthresh, 'min')

    # Offset temperature
    if toffset != 0:
        if verbose:
            print("Offsetting temperature by %f" % toffset)
        data["tair"] += toffset

    # Scale air pressure
    if pscale != 1:
        if verbose:
            print("Scaling air pressure by %f" % pscale)
        data["psurf"] *= pscale

    if precipscale != 1:
        if verbose:
            print("Scaling precipitation by %f" % precipscale)
        data["precip"] *= precipscale

    if precipscalegrid is not None:
        if verbose:
            print("Scaling precipitation by %s in %s" % (precipscalegridvar,
                                                         precipscalegridfile))
        data["precip"] = utils.scale_precip_grid(data["precip"],
                                                 dimvardata["month"],
                                                 precipscalegrid,
                                                 precipmonths)

    # Make the new array
    if ne == 0:
        pet = np.ma.masked_equal(np.ones([nt, ny, nx])*mv, mv)
    else:
        pet = np.ma.masked_equal(np.ones([ne, nt, ny, nx])*mv, mv)

    # Loop over timesteps
    for t in range(nt):
        if verbose:
            print("Timestep ", t)

        ########################################################
        # Calculate length of day
        if daynight:
            if verbose:
                print("Calculate day length")
            tmax_offset = 0.15
            time_up, time_down = solar.sun_times(dimvardata["dayno"][t],
                                                 dimvardata["year"][t],
                                                 dimvardata[latvar],
                                                 dimvardata[lonvar])

            time_max = solar.time_max_temperature(time_up, time_down,
                                                  tmax_offset)

            day_hours = time_down - time_up
            night_hours = const.daylen_hour - day_hours

        else:
            time_up = None
            time_down = None
            time_max = None
            day_hours = None
            night_hours = None

        ########################################################
        # MET CALCULATIONS
        lai = const.lai[dimvardata["month"][t]-1]
        if daynight:
            # Get day/night air temperature
            tair = get_daynight_tair(data["tair"][t, :],
                                     data["dtr"][t, :],
                                     time_up, time_down, time_max)
            ########################################################
            # Get day/night humidity, assuming relative humidity
            # is constant
            qair = get_daynight_qair(data["qair"][t, :],
                                     data["tair"][t, :],
                                     data["psurf"][t, :],
                                     tair[0], tair[1])

            if netradiation:
                sys.exit("Error: can't do separate day/night calculations "
                         "with net radiation")
            else:
                ######################################################
                # Get day/night SW down. All of the SW is in the day
                # and none at night (ignoring twilight)
                swdown = [data["swdown"][t, :] * const.daylen_hour / day_hours,
                          np.zeros_like(tair[1])]


                ########################################################
                # Get day/night LW down. Uses the assumed sinusoid as
                # in JULES
                lwdown = get_daynight_lwdown(data["lwdown"][t, :],
                                             data["dtr"][t, :],
                                             data["tair"][t, :],
                                             time_down, time_up, time_max)

            # Canopy resistance

            if verbose:
                print("    canopy resistance")
            if co2 is None:
                rsc_day = const.rsc_day[dimvardata["month"][t]-1]
            else:
                rsc_day = rs_fn_Rudd(co2[co2indx[t]],
                                     dimvardata["month"][t]-1,
                                     co2_baseline, lai,
                                     const.rsc_day[dimvardata["month"][t]-1])

            rs_day = canres_day(lai, rsc_day, const.rss)

            rs_night = canres_night(lai, const.rsc_night, const.rss)

            rs = [rs_day, rs_night]

            # Time period length
            Ds = [const.hourlen_sec*day_hours,
                  const.hourlen_sec*night_hours]

        else:
            tair = [data["tair"][t, :],]
            qair = [data["qair"][t, :],]
            if netradiation:
                swnet = [data["swnet"][t, :],]
                lwnet = [data["lwnet"][t, :],]
            else:
                swdown = [data["swdown"][t, :],]
                lwdown = [data["lwdown"][t, :],]

            ########################################################
            # canopy resistance
            if verbose:
                print("    canopy resistance")
            if co2 is None:
                rsc = const.rsc_day[dimvardata["month"][t]-1]
            else:
                rsc = rs_fn_Rudd(co2[co2indx[t]], dimvardata["month"][t]-1,
                                 co2_baseline, lai,
                                 const.rsc_day[dimvardata["month"][t]-1])

            rs = [canres_day(lai, rsc, const.rss),]

            # Time period length
            Ds = [const.daylen_sec,]


        if verbose:
            print("    available energy")
        if netradiation:
            # Radiation is the sum of the net components
            Rn = [lw+sw for lw, sw in zip(lwnet, swnet)]

        else:
            ########################################
            # Calculate albedo
            if interception:
                albedo = get_albedo(lai, const.albedo_c, const.albedo_s_dry,
                                    const.albedo_s_wet, data["precip"][t, :])
            else:
                albedo = get_albedo(lai, const.albedo_c, const.albedo_s_dry,
                                    const.albedo_s_wet,
                                    np.zeros_like(data["tair"][t, :]))

            ########################################
            # Calculate net radiation
            Rn = [net_radiation(sw, lw, ta, albedo, const.emiss)
                  for sw, lw, ta in zip(swdown, lwdown, tair)]

        if daynight:
            # Ground heat flux
            if groundheatflux == "morecs":
                G_day = ground_heat_flux_day(Rn[0])
                P = const.P[dimvardata["month"][t]-1]
                G_night = ground_heat_flux_night(G_day, time_up, time_down, P)
                Ghf = [G_day, G_night]
            elif groundheatflux == "zero":
                Ghf = [0.0, 0.0]
            else:
                sys.exit("Error: unknown groundheatflux <%s>"%groundheatflux)
        else:
            # ground heat flux
            if groundheatflux == "morecs":
                Ghf = [const.P[dimvardata["month"][t]-1]/24.0,]
            elif groundheatflux == "zero":
                Ghf = [0.0, ]
            else:
                sys.exit("Error: unknown groundheatflux <%s>"%groundheatflux)

        # Available energy is net radiation - ground heat flux
        AE = [R - G for R, G in zip(Rn, Ghf)]

        ########################################
        # Aerodynamic resistance from wind speed
        # wind measured at 10m
        if verbose:
            print("    aerodynamic resistance")
        ra = raero(data["wind10"][t, :], const.canht)

        ########################################
        # Get lambdaE (W m-2)
        # TRANSPIRATION
        if verbose:
            print("    transpiration")
        pet_t_p = [get_pet(ta, qa,
                          data["psurf"][t, :],
                          A, ra, rss, const.emiss, qdefnotneg,
                          docorr=docorr)
                          for ta, qa, A, rss in zip(tair, qair, AE, rs)]

        ########################################
        # Convert to mm/time period
        pet_t_p = [ptp*D/const.l for D, ptp in zip(Ds, pet_t_p)]

        if interception:
            if verbose:
                print("    interception")

            ########################################
            # INTERCEPTION
            # Intercepted proportion of rain
            interc = interc_rf(data["precip"][t, :], lai,
                               const.enhance[dimvardata["month"][t]-1])

            pet_i_p = [get_pet(ta, qa, data["psurf"][t, :],
                              A, ra, 0.0, const.emiss, qdefnotneg,
                              docorr=docorr)
                              for ta, qa, A in zip(tair, qair, AE)]

            ########################################
            # Convert to mm/time period
            pet_i_p = [pti*D/const.l for pti, D in zip(pet_i_p, Ds)]

            ########################################
            # Interception correction
            # if there has been no rain then
            # pet_t = pet_t_pot and pet_i = 0,
            # otherwise the correction is applied
            if verbose:
                print("    interception correction")
            if daynight:

                # Apply correction to daytime PET, assuming all rain
                # is in the canopy at the start of the day
                pet_t_day, pet_i_day = interc_corr(pet_t_p[0], pet_i_p[0],
                                                   interc, day_hours)

                # Reduce the canopy storage by the amount evaporated
                # during the day (can't be negative)
                interc_start_night = np.where(interc - pet_i_day > 0,
                                              interc - pet_i_day,
                                              0.0)

                # Apply correction to nighttime PET, using remaining
                # canopy storage
                pet_t_night, pet_i_night = interc_corr(pet_t_p[1],
                                                       pet_i_p[1],
                                                       interc_start_night,
                                                       night_hours)

                # Add day/night components
                pet_t = pet_t_day + pet_t_night
                pet_i = pet_i_day + pet_i_night

            else:
                # Apply canopy storage to total daily PET
                pet_t, pet_i = interc_corr(pet_t_p[0], pet_i_p[0],
                                           interc, const.daylen_hour)

            ########################################
            # add interception and transpiration
            pet[t, :] = pet_t + pet_i

        else:
            if daynight:
                # Add day and night components
                pet[t, :] = pet_t_p[0] + pet_t_p[1]
            else:
                # Use the total daily PET
                pet[t, :] = pet_t_p[0]

    # Rarrange the dimensions if required
    if ne > 0:
        if verbose:
            print("Re-swap axes")
        pet = pet.swapaxes(0, 1)

    ################################################
    # If asked, we set negative values to
    # zero
    if zneg:
        if verbose:
            print("    set negative values to zero")
        pet[np.where(np.logical_and(pet < 0, ~pet.mask))] = 0.0

    ########################################################
    # Create an output file
    if verbose:
        print("    creat output file")
        print("    ", outfilename)

    inf = nc.Dataset(datafiletmplt.substitute(varn=varns["tair"]), 'r')

    of = utils.copy_nc_file(inf, outfilename, exclude=[varns["tair"]],
                            doglobattrs=False)
    inf.close()

    ########################################################
    # PETI metadata
    outvarattrs = {
        outvarn: {
            "units": "mm/day",
            "standard_name": "water_potential_evaporation_amount",
            "cell_methods": "time: mean"
            }
        }
    if interception:
        outvarattrs[outvarn]["long_name"] = "Potential evapotranspiration with interception correction"
    else:
        outvarattrs[outvarn]["long_name"] = "Potential evapotranspiration"
    if gridmapvar is not None:
        outvarattrs[outvarn]["grid_mapping"] = gridmapvar

    utils.create_vars(of, {outvarn: outvarns[outvarn]}, {outvarn: outvarattrs[outvarn]})

    ########################################################
    # Global metadata
    globattrs = {
        "institution": "CEH Wallingford - NERC",
        "date_created": dt.datetime.now().strftime("%Y-%m-%d"),
        "version": version
    }
    if interception:
        globattrs["title"] = "Potential evapotranspiration with interception " \
                   "correction."
        globattrs["description"] = "Daily potential evapotranspiration with " \
                         "interception correction (mm)"
    else:
        globattrs["title"] = "Potential evapotranspiration"
        globattrs["description"] = "Daily potential evapotranspiration"

    if user is not None:
        globattrs["creator_name"] = user
    if email is not None:
        globattrs["creator_email"] = email

    for attn, attr in globattrs.items():
        of.setncattr(attn, attr)


    ########################################################
    # create PETi variable
    if verbose:
        print("    writing peti")
    if "_FillValue" not in outvarns[outvarn]:
        if np.ma.is_masked(pet):
            sys.exit("Error: not expecting masked data")
        else:
            of.variables[outvarn][:] = pet.data[:]
    else:
        of.variables[outvarn][:] = pet[:]

    of.close()
    if verbose:
        print("Done")

if __name__ == "__main__":
    main()
