import xarray as xr
import numpy as np
from os.path import exists
from os import remove
from scipy.stats import genlogistic
from scipy.special import ndtri


def load_nc_files(nc_paths: list | str) -> xr.DataArray:
    # Check input
    if not nc_paths:
        print("Path(s) to netCDF file(s) not provided.")  # TODO: logging
        return

    if isinstance(nc_paths, str):
        # Path was provided as string
        return xr.open_dataset(nc_paths)
    else:
        # Paths were provided as list
        return xr.open_mfdataset(nc_paths)


def save_as_nc(dataset: xr.Dataset, file_path: str) -> xr.DataArray:
    """ Stores xarray.DataArray as netCDF file at file path."""
    if exists(file_path):
        remove(file_path)
        print(f"Removed '{file_path}'")

    dataset.to_netcdf(file_path)
    print(f"Saved '{file_path}'")


def preproces_prec(prec_ds: xr.Dataset) -> xr.DataArray:
    prec_ds = prec_ds.transpose('time', 'lat', 'lon', 'latv', 'lonv', 'nv')
    prec_ds['time'] = prec_ds['time'].astype('datetime64[ns]')
    return prec_ds


def preprocess_temp(temp_ds: xr.Dataset) -> xr.DataArray:
    # Rename coordinate dimensions
    temp_ds = temp_ds.rename({"latitude": "lat", "longitude": "lon"})

    # Convert longitude coordinates from 0-360 to -180 to 180
    if (temp_ds.lon >= 180).any():
        print("Longitude conversion...")  # TODO: logging
        temp_ds.coords['lon'] = (temp_ds.coords['lon'] + 180) % 360 - 180
        temp_ds = temp_ds.sortby("lon")
        temp_ds = temp_ds.sortby("lat")
    temp_ds['time'] = temp_ds['time'].astype('datetime64[ns]')
    return temp_ds


def spatial_subset(ds: xr.Dataset, lat_bounds: list, lon_bounds: list) -> xr.DataArray:
    # Check input
    if len(lat_bounds) != 2:
        print(
            f"Latitude bounds should contain 2 values but you provided {len(lat_bounds)} values.")
        return

    if len(lon_bounds) != 2:
        print(
            f"Longitude bounds should contain 2 values but you provided {len(lon_bounds)} values.")
        return

    # Select spatial subset
    return ds.sel(lat=slice(*lat_bounds), lon=slice(*lon_bounds))


def calc_pet_thornthwaite(temp_ds: xr.Dataset, prec_ds: xr.Dataset) -> xr.DataArray:
    """
    Calculate PET using the Thornthwaite equation.
    temp: Monthly average temperature in Celsius.
    """
    # Convert temperature to Celsius
    temp_celsius = temp_ds - 273.15

    # Calculate I (heat index)
    # I = ((temp_celsius / 5.0) ** 1.514).sum('time')
    I = ((temp_celsius / 5.0) ** 1.514).groupby("time.month").sum('time')

    # Calculate a
    a = (6.75e-07) * I**3 - (7.71e-05) * I**2 + (1.792e-02) * I + 0.49239

    # Days in month
    days_in_month = temp_celsius.time.dt.days_in_month

    # Calculate PET
    pet_ds: xr.DataArray = (16 * (10 * temp_celsius / I)**a) * (days_in_month / 12)

    # Processing after calculation
    pet_ds: xr.DataArray = pet_ds.sel(expver=1)
    pet_ds = pet_ds.transpose('time', 'lat', 'lon', 'month')
    pet_ds['time'] = pet_ds['time'].astype('datetime64[ns]')
    pet_ds = pet_ds.interp(lat=prec_ds.lat, lon=prec_ds.lon, method='linear')

    return pet_ds


def calc_difference(prec_ds: xr.DataArray, pet_ds: xr.DataArray) -> xr.DataArray:
    return prec_ds - pet_ds


def calc_cdf_values2(diff_ds: xr.DataArray) -> xr.DataArray:
    diff_ds_vals = diff_ds.values
    diff_ds_flat = diff_ds_vals.flatten()
    diff_clean = diff_ds_flat[np.isfinite(diff_ds_flat)]
    params = genlogistic.fit(diff_clean)

    return genlogistic.cdf(diff_ds, *params)


def calc_cdf_values(diff_ds: xr.DataArray) -> xr.DataArray:
    # Convert the DataArray to a NumPy array and flatten it
    diff_ds_vals: np.ndarray = diff_ds.values
    print(type(diff_ds_vals))
    diff_ds_vals = diff_ds_vals.flatten()
    # Remove NaN or infinite values to clean the data
    diff_clean = diff_ds_vals[np.isfinite(diff_ds_vals)]
    # Fit the generalized logistic distribution to the cleaned data
    params = genlogistic.fit(diff_clean)

    # Apply the CDF to the original data using the fitted parameters
    # Note: We need to apply the CDF to each element in the original data array.
    # This involves creating a function that can be applied over the DataArray.
    cdf_applied = np.vectorize(lambda x: genlogistic.cdf(x, *params))

    # Since we can't directly apply np.vectorize result on xarray DataArray, we use .map() method of xarray if available,
    # or apply the function to the .values of the DataArray and then create a new DataArray from the result.
    # Ensure the output has the same shape as the original DataArray.
    cdf_values = xr.DataArray(cdf_applied(
        diff_ds.values), dims=diff_ds.dims, coords=diff_ds.coords)

    return cdf_values


def spei(diff_ds, cdf_values):
    return xr.DataArray(ndtri(cdf_values), coords=diff_ds.coords, dims=diff_ds.dims)
