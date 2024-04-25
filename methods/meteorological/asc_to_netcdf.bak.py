from osgeo import gdal, osr
import numpy as np
from netCDF4 import Dataset

# Given header information

# ASC file
header = {}
asc_file = '/home/jtrvz/Downloads/grids_germany_monthly_precipitation_202312.asc/RSMS_12_2023_01.asc'

# Extract header info
with open(asc_file, 'r') as file:
    counter = 0
    for line in file:
        if counter > 10:
            break
        # Split each line at whitespace
        parts = line.strip().split()
        # Check if the line contains one of the desired headers
        if parts[0] in ('NCOLS', 'NROWS', 'XLLCORNER', 'YLLCORNER', 'CELLSIZE', 'NODATA_VALUE'):
            # Store the value in the dictionary, converting it to float if necessary
            header[parts[0]] = float(parts[1]) if parts[0] != 'NODATA_VALUE' else int(parts[1])
        counter += 1

# Open ASC file
dataset = gdal.Open(asc_file)
band = dataset.GetRasterBand(1)
band_arr = band.ReadAsArray().astype(float)
band_arr[band_arr == header['NODATA_VALUE']] = np.nan

# Calculate latitude and longitude arrays based on the grid information
latitudes = np.linspace(header['YLLCORNER'], header['YLLCORNER'] + header['CELLSIZE'] * (header['NROWS'] - 1), header['NROWS'])
longitudes = np.linspace(header['XLLCORNER'], header['XLLCORNER'] + header['CELLSIZE'] * (header['NCOLS'] - 1), header['NCOLS'])

# Create netCDF file
nc_file = '/home/jtrvz/Downloads/grids_germany_monthly_precipitation_202312.asc/RSMS_12_2023_01.nc'
nc_dataset = Dataset(nc_file, 'w', format='NETCDF4')

# Define dimensions
lat_dim = nc_dataset.createDimension('lat', header['NROWS'])
lon_dim = nc_dataset.createDimension('lon', header['NCOLS'])

# Define variables
latitudes_var = nc_dataset.createVariable('lat', np.float32, ('lat',))
longitudes_var = nc_dataset.createVariable('lon', np.float32, ('lon',))
values = nc_dataset.createVariable('value', np.float32, ('lat','lon',), fill_value=header['NODATA_VALUE'])

# Assign values to latitude and longitude variables
latitudes_var[:] = latitudes
longitudes_var[:] = longitudes

# Write the data array
values[:] = np.flipud(band_arr)

# Close the netCDF file
nc_dataset.close()
