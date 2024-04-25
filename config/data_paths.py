"""
<1><2>_<3>_<4>_ROOT_PATH: Root path.
<1><2>_<3>_<4>_BASE_FILE_NAME: Base file name where a place holder like `date` is meant to be replaced using `.format(date=...)`.
    - <1>: Data source. (M = Meteorological, R = Remote sensing)
    - <2>: Duration. (H = Hourly, D = Daily, M = Monthly, A = Annualy)
    - <3>: Kind of data. (e.g. 'T_MIN' or 'PREC')
    - <4>: File extension. (e.g. 'ASC' or 'NC')
"""

# Meteorological data
## Monthly
### Temperature min
MM_T_MIN_ASC_ROOT_PATH = "/media/jtrvz/1tb/drought_data/temperature/dwd/min"
MM_T_MIN_ASC_BASE_FILE_NAME = "grids_germany_monthly_air_temp_min_{date}.asc"
### Temperature max
MM_T_MAX_ASC_ROOT_PATH = "/media/jtrvz/1tb/drought_data/temperature/dwd/max"
MM_T_MAX_ASC_BASE_FILE_NAME = "grids_germany_monthly_air_temp_max_{date}.asc"
### Temperature mean
MM_T_MEAN_ASC_ROOT_PATH = "/media/jtrvz/1tb/drought_data/temperature/dwd/avg"
MM_T_MEAN_ASC_BASE_FILE_NAME = "grids_germany_monthly_air_temp_mean_{date}.asc"
### Precipitation sum
MM_PREC_ASC_ROOT_PATH = "/media/jtrvz/1tb/drought_data/precipitation/dwd/avg"
MM_PREC_ASC_BASE_FILE_NAME = "grids_germany_monthly_precipitation_{date}.asc"
### Global radiation
MM_RS_KWH_ASC_ROOT_PATH = "/media/jtrvz/1tb/drought_data/radiation/dwd/global"
MM_RS_KWH_ASC_BASE_FILE_NAME = "grids_germany_monthly_radiation_global_{date}.asc"

# Remote sensing data
## Monthly
### Temperature min
RM_T_MIN_NC_ROOT_PATH = ""
RM_T_MIN_NC_BASE_FILE_NAME = ""
### Temperature max
RM_T_MAX_NC_ROOT_PATH = ""
RM_T_MAX_NC_BASE_FILE_NAME = ""
### Temperature mean
RM_T_MEAN_NC_ROOT_PATH = "/media/jtrvz/1tb/drought_data/temperature/era5/Global/monthly/netcdf/avg"
RM_T_MEAN_NC_BASE_FILE_NAME = "t2m_{date}.nc"
### Precipitation sum
RM_PREC_NC_ROOT_PATH = "/media/jtrvz/1tb/drought_data/precipitation/nasa_gpm/Global/monthly/netcdf/avg"
RM_PREC_NC_BASE_FILE_NAME = "3B-MO.MS.MRG.3IMERG.{full_date}-S000000-E235959.{month}.V07B.HDF5.nc4"
