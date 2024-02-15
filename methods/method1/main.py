from spei_calc import *

# Define dates
day = 1
date_begin = datetime(2014, 4, day)
date_end = datetime(2014, 4, day)

# Precipitation
prec_file_paths = generate_imerg_filenames(
    date_begin, date_end, "/media/jtrvz/1tb/drought_data/precipitation/nasa_gpm/Global/monthly/netcdf/avg/")
# Temperature
temp_file_paths = generate_t2m_filenames(
    date_begin, date_end, "/media/jtrvz/1tb/drought_data/temperature/era5/Global/monthly/netcdf/avg/")

# Load nc files
prec_ds = load_nc_files(prec_file_paths)
temp_ds = load_nc_files(temp_file_paths)

# Preprocess datasets
prec_ds = preprocess_prec(prec_ds)
temp_ds = preprocess_temp(temp_ds)

# Lat and lon bounds for Germany
lat_bounds = [47.0, 55.0]
lon_bounds = [5.5, 15.0]

# Spatial subset Germany
prec_ds = spatial_subset(prec_ds, lat_bounds, lon_bounds)
temp_ds = spatial_subset(temp_ds, lat_bounds, lon_bounds)

# Calculate PET (potential evapotranspiration)
pet_ds = calc_pet_thornthwaite(temp_ds)
pet_ds, prec_ds = preprocess_pet(pet_ds, prec_ds)

# Calculate difference
diff_ds = calc_difference(prec_ds, pet_ds)

# Calculate CDF values
cdf_vals = calc_cdf_values(diff_ds)

spei_da = spei(diff_ds, cdf_vals)

# Create plot
shape_file_path = '/home/jtrvz/Downloads/vg2500_geo84/vg2500_krs.shp'
spei_plot(spei_da, shape_file_path, lon_bounds, lat_bounds,
          f"SPEI Germany ({date_begin.strftime('%Y-%m')}-{date_end.strftime('%Y-%m')})", "test.png")