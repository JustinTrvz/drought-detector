# Lat and lon have an offset...
from osgeo import gdal, osr
import numpy as np
from netCDF4 import Dataset

# Given header information

# ASC file
header = {}
asc_file = 'grids_germany_monthly_precipitation_202312/RSMS_12_2023_01.asc'
prj_file = 'gk3.prj'
nc_file = asc_file.replace('.asc', '.nc')

# Read the spatial reference from the PRJ file
spatial_ref = osr.SpatialReference()
with open(prj_file, 'r') as file:
    prj_txt = file.read()
    spatial_ref.ImportFromWkt(prj_txt)

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
            header[parts[0]] = float(parts[1]) \
                if (parts[0] != 'NODATA_VALUE' or parts[0] != 'NROWS' or parts[0] != 'NCOLS' ) \
                    else int(parts[1])
        counter += 1

## Open ASC file with GDAL
dataset = gdal.Open(asc_file)
band = dataset.GetRasterBand(1)
band_arr = band.ReadAsArray().astype(float)
band_arr[band_arr == header['NODATA_VALUE']] = np.nan

# Assuming the dataset is in a projected CRS and needs to be transformed to EPSG:4326
target_cs = osr.SpatialReference()
target_cs.ImportFromEPSG(4326)
transform = osr.CoordinateTransformation(spatial_ref, target_cs)

# Transform the lower left corner
llx, lly, _ = transform.TransformPoint(header['XLLCORNER'], header['YLLCORNER'])
urx, ury, _ = transform.TransformPoint(header['XLLCORNER'] + header['CELLSIZE'] * header['NCOLS'], header['YLLCORNER'] + header['CELLSIZE'] * header['NROWS'])

# Calculate the new latitude and longitude arrays
latitudes = np.linspace(ury, lly, int(header['NCOLS'])) # Note: may need to adjust the order (ury to lly or lly to ury) based on the actual data orientation
longitudes = np.linspace(llx, urx, int(header['NROWS']))

# netCDF creation (as before), adjusted for lat/lon values
nc_dataset = Dataset(nc_file, 'w', format='NETCDF4')
lat_dim = nc_dataset.createDimension('lon', header['NCOLS'])
lon_dim = nc_dataset.createDimension('lat', header['NROWS'])
latitudes_var = nc_dataset.createVariable('lat', np.float32, ('lat',))
longitudes_var = nc_dataset.createVariable('lon', np.float32, ('lon',))
values = nc_dataset.createVariable('value', np.float32, ('lat','lon',), fill_value=np.nan)
latitudes_var[:] = longitudes + 0.18
longitudes_var[:] = latitudes - 0.38
values[:] = np.fliplr(np.flipud(band_arr))

nc_dataset.close()
