#!/usr/bin/env python
# Need everything except time to match
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
    parser.add_argument("-P", "--precipvar",
                        help="Precipitation variable",
                        required=False, default="precip")
    parser.add_argument("--latvar",
                        help="Latitude variable",
                        required=False, default="lat")
    parser.add_argument("--lonvar",
                        help="Longitude variable",
                        required=False, default="lon")
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
    parser.add_argument("-G", "--precipscalegridfilevarmn",
                        help="Precip scale grid file, variable and month variable", nargs=3,
                        default=[None, None, None])
    parser.add_argument("--timevar",
                        help="Time variable",
                        required=False, default="time")
    parser.add_argument("--gridmapvar",
                        help="Grid mapping variable",
                        required=False, default=None)
    parser.add_argument("-u", "--user",
                        help="Name of user who created this data",
                        required=False, default=None)
    parser.add_argument("-e", "--email",
                        help="Email of user who created this data",
                        required=False, default=None)

    # positional
    parser.add_argument("componentfiletmplt",
                        help="Template of files containing "
                             "monthly components ($varn)")
    parser.add_argument("precipfilename",
                        help="File containing daily precipitation")
    parser.add_argument("outfilename",
                        help="Output file name")

    args = parser.parse_args()

    precipvarn = {"precip": args.precipvar}

    return string.Template(args.componentfiletmplt), \
        args.precipfilename, \
        args.outfilename, \
        precipvarn, \
        args.zlib, \
        args.zeroneg, \
        args.version, \
        args.verbose, \
        args.precipscalegridfilevarmn[0], \
        args.precipscalegridfilevarmn[1], \
        args.precipscalegridfilevarmn[2], \
        args.precipscale, \
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
    componentfiletmplt, precipfilename, outfilename, precipvarn, \
        zlib, zneg, version, verbose, \
        precipscalegridfile, precipscalegridvar, \
        precipscalegridmn, precipscale, timevar, latvar, lonvar, \
        gridmapvar, user, email = parse_input()

    varns = {
            'pet': 'pet',
            'pei': 'pei'
            }

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
    # Read pet and pei
    if verbose:
        print("Read components")
    data, dims, dimvardata = utils.read_input_data(componentfiletmplt, varns,
                                                   verbose=verbose,
                                                   timevar=timevar)

    # get shape
    ne, nt, ny, nx = utils.get_shape(data["pet"])

    # Rearrange data if necessary
    if ne > 0:
        if verbose:
            print("Reorder data to cope with ensemble members")
        for varn in data:
            if len(data[varn].shape) == 4:
                data[varn] = data[varn].swapaxes(0, 1)

    ########################################################
    # Read precip
    if verbose:
        print("Read precip")
    pdata, pdims, pdimvardata = utils.read_input_data(string.Template(precipfilename), precipvarn,
                                                   verbose=verbose,
                                                   timevar=timevar)

    # get shape
    nep, ntp, nyp, nxp = utils.get_shape(pdata["precip"])
    if nep!=ne or nyp!=ny or nxp!=nx:
        # Need everything except time to match
        sys.exit("Error: shapes don't match")

    # Rearrange data if necessary
    if nep > 0:
        if verbose:
            print("Reorder precip to cope with ensemble members")
        for varn in pdata:
            if len(pdata[varn].shape) == 4:
                pdata[varn] = pdata[varn].swapaxes(0, 1)

    # get attributes
    outvarn = "peti"
    outvarns = {
        outvarn: utils.get_var_attrs(precipfilename,
        precipvarn["precip"])
    }

    # Check fill value
    if '_FillValue' in outvarns["peti"]:
        mv = outvarns["peti"]['_FillValue']
    else:
        mv = -1e20

    # scale precip
    if precipscale != 1:
        if verbose:
            print("Scaling precipitation by %f" % precipscale)
        pdata["precip"] *= precipscale

    if precipscalegrid is not None:
        if verbose:
            print("Scaling precipitation by %s in %s" % (precipscalegridvar,
                                                         precipscalegridfile))
        pdata["precip"] = utils.scale_precip_grid(pdata["precip"],
                                                 pdimvardata["month"],
                                                 precipscalegrid,
                                                 precipmonths)

    # Make the new array (need daily, so copy the precip shape)
    if nep == 0:
        peti = np.ma.masked_equal(np.ones([ntp, nyp, nxp])*mv, mv)
    else:
        peti = np.ma.masked_equal(np.ones([ntp, nep, nyp, nxp])*mv, mv)

    # Loop over timesteps (daily, so precip)
    for t in range(ntp):

        if verbose:
            print("Timestep ", t)

        # find which monthly pet/pei to use
        mn = np.where(np.logical_and(dimvardata["month"] == pdimvardata["month"][t],dimvardata["year"] == pdimvardata["year"][t]))

        if len(mn) > 1:
            sys.exit("Error: Too many indices")
        elif len(mn) == 0 or len(mn[0]) == 0:
            sys.exit("Error: month not found")
        elif len(mn[0]) > 1:
            sys.exit("Error: Found too many possible months")
        else:
            mn = mn[0][0]

        # get lai
        lai = const.lai[pdimvardata["month"][t]-1]

        ########################################
        # INTERCEPTION
        # Intercepted proportion of rain
        interc = interc_rf(pdata["precip"][t, :], lai,
                           const.enhance[pdimvardata["month"][t]-1])

        ########################################
        # Interception correction
        # if there has been no rain then
        # pet_t = pet_t_pot and pet_i = 0,
        # otherwise the correction is applied
        if verbose:
            print("    interception correction")

        # Apply canopy storage to total daily PET
        pet_t, pet_i = interc_corr(data["pet"][mn], data["pei"][mn],
                                   interc, const.daylen_hour)

        ########################################
        # add interception and transpiration
        peti[t, :] = pet_t + pet_i

    # Rarrange the dimensions if required
    if ne > 0:
        if verbose:
            print("Re-swap axes")
            peti = peti.swapaxes(0, 1)

    ################################################
    # If asked, we set negative values to
    # zero (only for total peti even if outputting
    # components)
    if zneg:
        if verbose:
            print("    set negative values to zero")
        peti[np.where(np.logical_and(peti < 0, ~peti.mask))] = 0.0

    ########################################################
    # Create an output file
    if verbose:
        print("    creat output file")
        print("    ", outfilename)

    inf = nc.Dataset(precipfilename, 'r')

    of = utils.copy_nc_file(inf, outfilename, exclude=[precipvarn["precip"]],
                            doglobattrs=False)
    inf.close()

    ########################################################
    # PETI metadata
    outvarattrs = {
        outvarn: {
            "units": "mm/day",
            "standard_name": "water_potential_evaporation_amount",
            "cell_methods": "time: mean",
            "long_name": "Potential evapotranspiration with interception correction"
            }
        }

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
    globattrs["title"] = "Potential evapotranspiration with interception " \
               "correction."
    globattrs["description"] = "Daily potential evapotranspiration with " \
                     "interception correction (mm)"

    if user is not None:
        globattrs["creator_name"] = user
    if email is not None:
        globattrs["creator_email"] = email

    for attn, attr in globattrs.items():
        of.setncattr(attn, attr)


    ########################################################
    # create PETI variable
    if verbose:
        print("    writing peti")
    if "_FillValue" not in outvarns[outvarn]:
        if np.ma.is_masked(peti):
            sys.exit("Error: not expecting masked data")
        else:
            of.variables[outvarn][:] = peti.data[:]
    else:
        of.variables[outvarn][:] = peti[:]



    of.close()
    if verbose:
        print("Done")

if __name__ == "__main__":
    main()
