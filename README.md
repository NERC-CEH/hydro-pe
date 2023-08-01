# hydro-pe
A package to calculate Hydro-PE potential evapotranspiration from netCDF inputs.

# Package contents

This provides three scripts that can be called on the command line:

- `make_pet` calculates potential evapotranspiration (PET) from meteorological data, optionally with an interception correction (PETI)
- `mask_file` applies a land-sea mask to a file and fills in points that have been incorrectly categorised as sea with values from 'donor' land points
- `combine_peti_components` combines monthly potential evapotranspiration (PET) and potential interception (PEI) rates to calculate daily PETI by scaling them in proportion to the amount of daily precipitation

It also provides several libraries:

- `pet` contains functions used to calculate PET and PETI from the input meteorological variables
- `pet.constants` contains useful constant values
- `pet.solar` functions to calculate solar variables from latitude and longitude
- `pet.utils` contains functions for input and output file handling

# Installation

    python ./setup.py build
    python ./setup.py install

# Hydro-PE calculations

Details of the calculations are available in the associated data paper in Earth System Science Data [[1]](#references).

The output units of PET and PETI are mm d⁻¹.

Hydro-PE calculates Penman-Monteith potential evapotranspiration (PET) parameterised for a short grass surface, with an optional correction for interception (PETI). It can be applied to observation based meteorological data sets or to the output of weather and climate models. Much of the methodology and parameterisations are based on the equations and parameters in MORECS 2.0 [[2]](#references), which was developed by the UK Met Office. This is motivated by the historical use of MORECS as a driver for hydrological models and allows for the calculation of equivalent PETI with other meteorological inputs.

In addition, a tool is provided to combine monthly potential evapotranspiration and potential interception rates with daily precipitation to calculate daily PETI. This is useful in cases where precipitation is available daily, but other variables are only available monthly, for example for some climate model outputs.

## Useage: `make_pet`

The basic useage is

`make_pet DATAFILETEMPLATE OUTFILENAME`

- `DATAFILETEMPLATE` is either the name of the single input file containing all the variables, or a template for the individual files containing the variables, with the variable name substituted with `${varn}`
- `OUTFILENAME` is the name of the output file

## Options: `make_pet`

`--interception` / `-i` The default behaviour is to calculate PET. This option enables the interception correction.

`--docorr` / `-c` If net longwave radiation is not known, then the upward component can be approximated by calculating it using air temperature as a proxy for surface temperature. In this case, a correction to the Penman-Monteith equation is applied [[2]](#references). If the radiation inputs are net radiation variables that have been calculated using this approximation, then this flag enables the Penman-Monteith correction, otherwise net radiation inputs will be assumed to be calculated with surface temperature. If the radiation inputs are downward radiation variables, then this flag is overridden, as the upward longwave is calculated using the air temperature approximation in the code, so the Penman-Monteith correction is required.

`--netsw` The input shortwave radiation is net radiation. Otherwise, the input is downward shortwave radiation.

`--netlw` The input longwave radiation is net radiation. Otherwise, the input is downward longwave radiation.

`--netradiation` Both the input short- and longwave radiation are net radiation. Otherwise they may be specified individually.

`--qdefnotneg` / `-Q` The Penman-Monteith calculation uses specific humidity, following Stewart (1989) [[3]](#references). The specific humidity at saturation is calculated using an empirical fit to vapour pressure at saturation [[4]](#references). In cases where the input specific humidity is greater than the calculated specific humidity at saturation, then this option can be used to set the humidity defecit to zero.

`--co2filevaryr FILENAME VARNAME YEAR` / `-C FILENAME VARNAME YEAR` In order to account for the response of plant stomatal conductance to rising CO₂ levels, a fertilisation effect may be applied [[5]](#references) by using this option to specify 
- `FILENAME` a netCDF file containing annual CO₂ concentrations in PPMV
- `VARNAME` the variable name in the file
- `YEAR` the year to use as the baseline (before which the stomatal conductance will not be adjusted).

`--daynight` / `-D` The PET and PETI may be calcuated directly from daily mean values, or it may optionally be separated into daytime and nighttime components by applying a diurnal cycle [[6]](#references). This option `is not available if net radiation inputs are used. 

`--zeroneg` / `-z` If the output PET or PETI is negative, this option will set all negative values to zero.

`--groundheatflux OPTION` This option allows selection of the method to calculate ground heat flux. The available options are 
- `morecs` calculate ground heat flux following MORECS 2.0 (default)
- `zero` does not calculate ground heat flux.

`--windthresh THRESHOLD` / `-W THRESHOLD` Applies a minimum threshold to the wind speed variable. All wind speed values less than `THRESHOLD` are set equal to the threshold value.

`--toffset OFFSET` / `-T OFFSET` Applies an offset to the air temperature variable. If the input air temperature is in °C then an offset of 273.15 can be used to convert to the expected input units of K.

`--pscale SCALE` / `-S SCALE` Applies a scale factor to the surface air pressure variable. If the input surface air pressur is in hPa then a scale of 100 can be used to convert to the expected input units of Pa.

`--precipscale SCALE` / `-L SCALE` Applies a scale factor to the precipitation variable. If the input precipitation is in kg m⁻² s⁻¹ then a scale of 86400 can be used to convert to the expected input units of mm d⁻¹.

`--precipscalegridfilevarmn FILENAME VARNAME MNVARNAME` / `-G FILENAME VARNAME MNVARNAME` Applies a grid of monthly scale factors to the precipitation variable. The scale factors are defined in the file FILENAME, in the variable VARNAME. The variable containing the month is MNVARNAME.

## Input variables: `make_pet`

Inputs to `make_pet` must be in netCDF format. Variables may all be in the same file, or may be in individual files, with file names tempated on the variable name `${varn}`. The code requires each variale to have the expected units. Some variables may be input in other units, and then a specified offset or scale factor applied. If the variable names are not specified, then a default value is used.

The radiation variables may either be input as net or downward. In the latter case the upward components are calculated in the code. If net radiation is input, then the split into day- and nighttime calculations cannot be used.

| Variable                     | Units     | Long option   | Short option | Default name | Scale/offset           | Notes                                                |
|:-----------------------------|:----------|:--------------|:-------------|:-------------|:-----------------------|:-----------------------------------------------------|
| Daily mean air temperature   | K         | `--tairvar`   | `-t`         | `tas`        | `--toffset` / `-T`     |                                                      |
| Daily mean specific humidity | kg kg⁻¹   | `--qairvar`   | `-q`         | `huss`       |                        |                                                      |
| Daily mean wind speed at 10m | m s⁻¹     | `--windvar`   | `-w`         | `sfcWind`    |                        | Minimum threshold specified by `--windthresh` / `-W` |
| Downward shortwave radiation | W m⁻²     | `--swvar`     | `-s`         | `rsds`       |                        | Only used if `--netradiation` is not specified       |
| Net shortwave radiation      | W m⁻²     | `--swvar`     | `-s`         | `rsds`       |                        | Only used if `--netradiation` is specified           |
| Downward longwave radiation  | W m⁻²     | `--lwvar`     | `-l`         | `rlds`       |                        | Only used if `--netradiation` is not specified       |
| Net longwave radiation       | W m⁻²     | `--lwvar`     | `-l`         | `rlds`       |                        | Only used if `--netradiation` is specified           |
| Surface air pressure         | Pa        | `--psurfvar`  | `-p`         | `psurf`      | `--pscale` / `-S`      |                                                      |
| Daily precipitation          | mm d⁻¹    | `--precipvar` | `-P`         | `precip`     | `--precipscale` / `-L` | Only used if `--interception` / `-i` is specified    |
| Daily temperature range      | K         | `--dtrvar`    | `-d`         | `dtr`        |                        | Only used if `--daynight` / `-D` is specified        |

In addition, the names of some other variables may be required.

`--timevar VARNAME` Specifies the name of the time variable in the input file. If not specified, the default value is `time`.

`--latvar VARNAME` Specifies the name of the latitude variable in the input file. If not specified, the default value is `lat`. This is only required if `--daynight` / `-D` is used.

`--lonvar VARNAME` Specifies the name of the longitude variable in the input file. If not specified, the default value is `lon`. This is only required if `--daynight` / `-D` is used.

`--gridmapvar VARNAME` Specifies the name of the grid mapping variable in the input file to ensure it is applied to the output netCDF file. If not specified, then the grid mapping variable is not used.

## Useage: `combine_peti_components`

The basic useage is

`combine_peti_components COMPONENTFILETEMPLATE PRECIPFILENAME OUTFILENAME`

where `COMPONENTFILETEMPLATE` contains monthly PET and PEI (usually an output file from `make_pet`), `PRECIPFILENAME` is the name of the file containing daily precipitation and `OUTFILENAME` is the name of the output file.

## Input variables: `combine_peti_components`

The PET and PEI components are assumed to be monthly mean values, in units of mm d⁻¹. The precipitation is specified as for `make_pet`, and may be scaled in the same way.

## Output files
The output of both `make_pet` and `combine_peti_components` is written into netCDF files, the structure of which is copied from the input netCDF files. Extra options may be specified.

`--outputcomponents` / `-O` When using `--interception` with `make_pet`, this outputs the potential evapotranspiration (PET) and potential interception (PEI) separately, without accounting for the precipitation, instead of outputting the combined PETI. These can later be used as input to `combine_peti_components` to calculate PETI using a precipitation data set.

`--fullout` / `-f` If `--daynight` is used, then this option will enable output of the day- and nighttime values, as well as the overall PET or PETI.

`--version VERSION` / `-V VERSION` Version number to be written in the metadata of the output netCDF file. Default is v1.0.

`--user USER` / `-u USER` User name to be written in the metadata of the output netCDF file. If not specified, then no user will be written.

`--email EMAIL` / `-e EMAIL` Email address to be written in the metadata of the output netCDF file. If not specified, then no email will be written.

# Acknowledgement 

This work was supported by the Natural Environment Research Council award number NE/S017380/1 as part of the Hydro-JULES programme <https://hydro-jules.org/>.

# References

[1] Robinson, E. L., Brown, M. J., Kay, A. L., Lane, R. A., Chapman, R., Bell, V. A., and Blyth, E. M.: Hydro-PE: gridded datasets of historical and future Penman-Monteith potential evaporation for the United Kingdom, Earth Syst. Sci. Data Discuss. [preprint], <https://doi.org/10.5194/essd-2022-288>, in review, 2022.

[2] Hough, M., Palmer, S., Weir, A., Lee, M., and Barrie, I.: The Meteorological Office Rainfall and Evaporation Calculation System: MORECS version 2.0, An update to Hydrological Memorandum 45, Met Office, <https://digital.nmla.metoffice.gov.uk/IO_9d68dec6-8ad2-420b-a971-806f7a6987d8/>, 1997.

[3] Stewart, J. B.: On the use of the Penman-Monteith equation for determining areal evapotranspiration, in: Estimation of Areal Evapotranspiration (Proceedings of a workshop held at Vancouver, B.C., Canada, August 1987), IAHS, Wallingford, Oxfordshire, UK, 1989.

[4] Richards, J. M.: A simple expression for the saturation vapour pressure of water in the range -50 to 140°C, Journal of Physics D: Applied Physics, 4, L15--L18, <https://doi.org/10.1088/0022-3727/4/4/101>, 1971

[5] Kruijt, B., Witte, J.-P. M., Jacobs, C. M., and Kroon, T.: Effects of rising atmospheric CO₂ on evapotranspiration and soil moisture: A practical approach for the Netherlands, Journal of Hydrology, 349, 257--267, <https://doi.org/https://doi.org/10.1016/j.jhydrol.2007.10.052>, 2008

[6] Williams, K. and Clark, D.: Disaggregation of daily data in JULES, Hadley Centre Technical Note 96, Met Office, <https://digital.nmla.metoffice.gov.uk/IO_8059be68-58cf-46b3-9492-1b1a0c290c0c/>, 2014.
