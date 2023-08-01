#!/usr/bin/env python

import netCDF4 as nc
import numpy as np
import calendar as _cal
import sys as _sys

################################################################################
################################################################################
#
#
#
# ELR 29-07-2020
#
################################################################################
################################################################################


################################################################################
# Check number of days in year
################################################################################
def days_in_year(year, calendar='gregorian'):
    if calendar == 'gregorian':
        if _cal.isleap(year):
            return 366
        else:
            return 365
    elif calendar == '360_day':
        return 360
    else:
        sys.exit("Error: unknown calendar")


def days_in_year_arr(year, calendar='gregorian'):
    return np.array([days_in_year(yr, calendar=calendar) for yr in year])


################################################################################
def create_nc_file(outfname, dims, varns, varattrs, globattrs):

    outf = nc.Dataset(outfname, 'w')

    create_dims(outf, dims)

    create_vars(outf, varns, varattrs)

    for attn, attr in globattrs.items():
        outf.setncattr(attn, attr)

    return outf

def create_dims(outf, dims):

    for dim, dlen in dims.items():
        outf.createDimension(dim,dlen)

    return

def create_vars(outf, varns, varattrs):

    for varn, varstuff in varns.items():
        if '_FillValue' in varstuff:
            outf.createVariable(varn, varstuff['dtype'], varstuff['dimensions'],
                    fill_value = varstuff['_FillValue'])
        else:
            outf.createVariable(varn, varstuff['dtype'], varstuff['dimensions'])

    for varn, varattr in varattrs.items():
        for attn, attr in varattr.items():
            outf.variables[varn].setncattr(attn, attr)

    return

def copy_nc_file_structure(inf, outfname, exclude=[]):

    dims = {dim: len(inf.dimensions[dim]) for dim in inf.dimensions}

    varns = {varn:
        {'dtype': inf.variables[varn].dtype,
                    'dimensions': inf.variables[varn].dimensions,
                    '_FillValue': inf.variables[varn]._FillValue}
        if '_FillValue' in inf.variables[varn].ncattrs() else
        {'dtype': inf.variables[varn].dtype,
                    'dimensions': inf.variables[varn].dimensions}
                    for varn in inf.variables
                if varn not in exclude}

    varattrs = {varn: {attn: inf.variables[varn].getncattr(attn)
                       for attn in inf.variables[varn].ncattrs()
                       if attn not in ['_FillValue',]}
                for varn in inf.variables
                if varn not in exclude}

    outf = create_nc_file(outfname, dims, varns, varattrs, {})

    return outf

def copy_nc_file(inf, outfname, nodat=[], exclude=[], doglobattrs=True):


    outf = copy_nc_file_structure(inf, outfname, exclude=exclude)

    for varn in inf.variables:
        if varn not in nodat and varn not in exclude:
            outf.variables[varn][:] = inf.variables[varn][:]

    if doglobattrs:
        globattrs = {attn: inf.getncattr(attn) for attn in inf.ncattrs()}
        for attn, attr in globattrs.items():
            outf.setncattr(attn, attr)

    return outf

def read_points_file(infname):
    inf = open(infname,'r')

    hdr = inf.readline()

    dat = [[int(el) for el in lin.strip().split(',')[:4]] for lin in inf.readlines()]
    fillx, filly, donorx, donory = [l for l in zip(*dat)]

    return filly, fillx, donory, donorx

def fill_points(dat, filly, fillx, donory, donorx):

    warnpts = []
    fillpts = 0
    for fy, fx, dy, dx in zip(filly, fillx, donory, donorx):
        fillpts+=1
        if np.ma.is_masked(dat[..., dy, dx]):
            _sys.exit("Error donor point y=%d, x=%d has some masking"%(dy, dx))
        if not np.ma.is_masked(dat[..., fy, fx]):
            print("Warning fill point y=%d, x=%d already has data in it"%(fy, fx))
            warnpts.append([fy, fx])

        dat[..., fy, fx] = dat[..., dy, dx]

    print("Filled %d points"%fillpts)
    if len(warnpts)>0:
        print("Warned %d points"%len(warnpts))
        print(warnpts)
    return

def read_co2_file(infname, varn, co2yr, timevar='time'):

    cf = nc.Dataset(infname, "r")
    co2 = cf.variables[varn][:]

    co2yrs = [t.year for t in nc.num2date(cf.variables[timevar][:],
                                          cf.variables[timevar].units,
                                          cf.variables[timevar].calendar)]
    co2yrs = np.array(co2yrs)

    iyr = np.where(co2yrs == co2yr)
    if len(iyr[0]) < 1:
        sys.exit("Error: year <%d> not found in CO2 data" % co2yr)
    elif len(iyr[0]) > 1:
        sys.exit("Error: year <%d> found too many "
                 "(%d) times in CO2 data!" % (co2yr, len(iyr[0])))

    iyr = iyr[0][0]

    # Make the years before the baseline year equal to the baseline
    # so that the r_s doesn't change
    co2_baseline = co2[iyr]
    co2[:iyr] = co2_baseline
    cf.close()

    return co2, co2yrs, co2_baseline

def read_precipscalegrid_file(infname, varn, monthvar='Months'):

    pf = nc.Dataset(infname, "r")
    precipscalegrid = pf.variables[varn][:]

    months = pf.variables[monthvar][:]

    if np.any(months<1) or np.any(months>12):
        sys.exit("Error: months must be an integer in the range [1,12]")
    pf.close()

    return precipscalegrid, months

def scale_precip_grid(indata, months, scalegrid, scalemonths):

    if np.ma.isMA(indata) or np.ma.isMA(scalegrid):
        outdata = np.ma.array([ind*scalegrid[scalemonths==mn][0] for ind, mn in zip(indata, months)])
    else:
        outdata = np.array([ind*scalegrid[scalemonths==mn][0] for ind, mn in zip(indata, months)])

    return outdata

def read_input_data(datafiletmplt, varns, verbose=False,
                    timevar="time", exclude=[]):
    data = {}
    dims = {}
    dimvardata = {}
    read_dimvarns = True
    for varn, invarn in varns.items():
        if verbose:
            print("    read ", varn)

        ################################################
        # Read the file
        fname = datafiletmplt.substitute(varn=invarn)
        f = nc.Dataset(fname, "r")

        ################################################
        # Get the relevant data
        data[varn] = f.variables[invarn][:]

        ########################################
        # Get the shape of the data
        if len(data[varn].shape) == 4:
            if data[varn].shape[0] == 1:
                data[varn] = data[varn][0, :]
            else:
                sys.exit("Can't currently handle this sort of data")

        ################################################
        # Get dimension variables from one of the files
        if read_dimvarns:
            for dvarn in f.variables:
                if dvarn not in varns:
                    if verbose:
                        print("    ", dvarn)

                    dimvardata[dvarn] = f.variables[dvarn][:]

                    if dvarn == timevar:
                        nt = len(f.dimensions[dvarn])
                        dimvardata["datetime"] = \
                            nc.num2date(dimvardata[dvarn][:],
                                        f.variables[dvarn].units,
                                        f.variables[dvarn].calendar)
                        dimvardata["dayno"] = np.array([t.dayofyr
                                                  for t in dimvardata["datetime"]])
                        dimvardata["month"] = np.array([t.month
                                                  for t in dimvardata["datetime"]])
                        dimvardata["year"] = np.array([t.year
                                                 for t in dimvardata["datetime"]])
                        if "calendar" in f.variables[dvarn].ncattrs():
                            calendar = f.variables[dvarn].calendar
                        else:
                            calendar = "gregorian"

                        tunits = f.variables[dvarn].units


            for dim in f.variables[invarn].dimensions:
                if dim not in dims:
                    dims[dim] = len(f.dimensions[dim])

            # Check that fill values are consistent. Some versions of
            # netCDF4 apply fill values to variables that don't have a fill
            # value in the file, and this can cause problems with arithmetic
            # at a later point if it's not consistent with the defined fill
            # value.
            for dvarn in dimvardata:
                if np.ma.isMA(dimvardata[dvarn]) and np.ma.isMA(data[varn]) \
                        and dimvardata[dvarn].shape != ():
                    if dimvardata[dvarn].fill_value != data[varn].fill_value:
                        dimvardata[dvarn] = np.ma.array(dimvardata[dvarn].data,
                                            mask=dimvardata[dvarn].mask,
                                            fill_value=data[varn].fill_value)

            read_dimvarns = False

        ################################################
        # Close file
        f.close()
    return data, dims, dimvardata

def get_shape(data):

    if len(data.shape) == 4:
        ne, nt, ny, nx = data.shape
    else:
        nt, ny, nx = data.shape
        ne = 0

    return ne, nt, ny, nx

def get_var_attrs(fname, varn):
    f = nc.Dataset(fname, 'r')
    attrs = {
        'dtype': f.variables[varn].dtype,
        'dimensions': f.variables[varn].dimensions
        }
    if "_FillValue" in f.variables[varn].ncattrs():
        attrs['_FillValue'] = f.variables[varn]._FillValue

    f.close()


    return attrs

def get_co2_indx(co2yrs, years):
    return np.array([np.where(co2yrs == year)[0][0]
                        for year in years])
