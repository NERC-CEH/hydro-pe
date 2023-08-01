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

from .pet import *
from . import constants as const
from . import utils

###############################################################################
# Parse the input
###############################################################################
def parse_input():
    parser = argparse.ArgumentParser(description='Mask file and copy points as necessary')

    # optional
    parser.add_argument('--pointsfile', help='Point map file')

    # positional
    parser.add_argument('maskfile', help='Mask file')
    parser.add_argument('maskvarn', help='Mask variable')
    parser.add_argument('invarn', help='Input variable')
    parser.add_argument('infiles', help='Input files', nargs='+')

    args = parser.parse_args()

    return \
            args.maskfile, \
            args.maskvarn, \
            args.invarn, \
            args.infiles, \
            args.pointsfile

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
    maskfile, maskvarn, invarn, infiles, pointsfile = parse_input()

    mf = nc.Dataset(maskfile, 'r')
    mask = mf.variables[maskvarn][:].mask

    if pointsfile is not None:
        ny, nx = mask.shape
        filly, fillx, donory, donorx = utils.read_points_file(pointsfile)
        filly = [ny - y - 1 for y in filly]
        donory = [ny - y - 1 for y in donory]

    mask = mask[np.newaxis, ...]

    for infile in infiles:
        inf = nc.Dataset(infile, 'r')
        if infile[-3:] != '.nc':
            sys.exit('Error unrecognised file name <%s>'%infile)
        outfname = infile[:-3] + '_masked.nc'

        outf = utils.copy_nc_file(inf, outfname, nodat=[invarn,])

        outvarn = np.ma.masked_array(inf.variables[invarn][:])

        outvarn.mask = mask

        if pointsfile is not None:
            utils.fill_points(outvarn, filly, fillx, donory, donorx)

        outf.variables[invarn][:] = outvarn

        outf.close()






###############################################################################
if __name__ == '__main__':
    main()
