import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
from os.path import exists
from os import remove
from scipy.stats import genlogistic, norm, fisk
from scipy.stats import genextreme as gev
from scipy.special import ndtri
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from os.path import exists
from typing import Literal


def generate_imerg_filenames(start_date: datetime, end_date: datetime, dir: str = ""):
    """
    Generate a list of IMERG filenames for each month between start_date and end_date.

    Parameters:
    - start_date (tuple): Start date as (year, month, day).
    - end_date (tuple): End date as (year, month, day).

    Returns:
    - list of strings: Filenames for each month in the specified date range.
    """

    # Ensure the start date is the first of the month
    start_date = start_date.replace(day=1)

    # Add trailing slash if needed
    if dir:
        if dir[-1] != "/":
            dir += "/"

    file_names = []

    current = start_date
    while current <= end_date:
        # Generate filename
        file_name = f"{dir}3B-MO.MS.MRG.3IMERG.{current.strftime('%Y%m%d')}-S000000-E235959.{current.strftime('%m')}.V07B.HDF5.nc4"
        if exists(file_name):
            file_names.append(file_name)
        else:
            print(f"File '{file_name}' does not exist.")
            raise Exception(f"File '{file_name}' does not exist.")

        # Move to the next month
        # Adding a month to a datetime object is tricky because not all months have the same number of days.
        # This method handles the end-of-month edge cases.
        next_month = current.replace(
            day=28) + timedelta(days=4)  # this will never fail
        current = next_month - timedelta(days=next_month.day - 1)

    return file_names


def generate_t2m_filenames(start_date: datetime, end_date: datetime, dir: str = "", file_prefix: str = "t2m", ):
    file_names = []

    # Ensure the start date is the first of the month
    start_date = start_date.replace(day=1)

    # Add trailing slash if needed
    if dir:
        if dir[-1] != "/":
            dir += "/"

    current = start_date
    while current <= end_date:
        # Generate filename
        file_name = f"{dir}{file_prefix}_{current.strftime('%Y%m')}.nc"
        if exists(file_name):
            file_names.append(file_name)
        else:
            print(f"File '{file_name}' does not exist.")
            raise Exception(f"File '{file_name}' does not exist.")

        # Move to the next month
        # Adding a month to a datetime object is tricky because not all months have the same number of days.
        # This method handles the end-of-month edge cases.
        next_month = current.replace(
            day=28) + timedelta(days=4)
        current = next_month - timedelta(days=next_month.day - 1)

    return file_names


def read_nc_files(nc_paths: list | str) -> xr.DataArray:
    # Check input
    if not nc_paths:
        print("Path(s) to netCDF file(s) not provided.")  # TODO: logging
        return

    if isinstance(nc_paths, str):
        # Path was provided as string
        return xr.open_dataset(nc_paths)
    elif len(nc_paths) == 1:
        # List with one path
        return xr.open_dataset(nc_paths[0])
    else:
        # Paths were provided as list
        return xr.open_mfdataset(nc_paths, combine='by_coords')


def save_as_nc(dataset: xr.Dataset, file_path: str) -> bool:
    """ Stores xarray.DataArray as netCDF file at file path."""
    if exists(file_path):
        remove(file_path)
        print(f"Removed '{file_path}'")

    print(f"Saving '{file_path}' ...")
    dataset.to_netcdf(file_path)

    file_exists = exists(file_path)
    if not file_exists:
        print(f"Failed to save '{file_path}'.")

    return file_exists


def preprocess_prec(prec_ds: xr.Dataset) -> xr.Dataset:
    prec_ds = prec_ds.transpose('time', 'lat', 'lon', 'latv', 'lonv', 'nv')
    prec_ds['time'] = prec_ds['time'].astype('datetime64[ns]')
    return prec_ds


def preprocess_temp(temp_ds: xr.Dataset) -> xr.DataArray:
    # Rename coordinate dimensions
    temp_ds = temp_ds.rename({"latitude": "lat", "longitude": "lon"})
    temp_ds = temp_ds.sel(expver=1)


    # Convert time dimensions type to datetime64[ns]
    temp_ds['time'] = temp_ds['time'].astype('datetime64[ns]')

    # Convert to Celsius
    temp_ds['t2m'] = temp_ds['t2m'] - 273.15

    # Convert longitude coordinates from 0-360 to -180 to 180
    if (temp_ds.lon >= 180).any():
        temp_ds.coords['lon'] = (temp_ds.coords['lon'] + 180) % 360 - 180
        temp_ds = temp_ds.sortby("lon")
        temp_ds = temp_ds.sortby("lat")

    return temp_ds


def spatial_subset(ds: xr.Dataset, lat_bounds: list, lon_bounds: list) -> xr.Dataset:
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

def calc_pet_thornthwaite(temp_ds: xr.Dataset, last_year_temp_ds: xr.Dataset):
    """
    Calculate PET using the Thornthwaite equation from a Dataset with
    monthly average 2m temperature in Celsius.
    """
    # Ensure all temperatures are non-negative, as PET calculation does not consider negative temperatures.
    nan_mask = temp_ds["t2m"].isnull()
    temp_ds['t2m'] = temp_ds['t2m'].where(cond=temp_ds['t2m'] >= 0, other=0)
    temp_ds['t2m'] = temp_ds['t2m'].where(cond=~nan_mask, other=np.nan)
    # Same for last year temperature dataset
    nan_mask = last_year_temp_ds["t2m"].isnull()
    last_year_temp_ds['t2m'] = last_year_temp_ds['t2m'].where(cond=last_year_temp_ds['t2m'] >= 0, other=0)
    last_year_temp_ds['t2m'] = last_year_temp_ds['t2m'].where(cond=~nan_mask, other=np.nan)

    # Calculate the heat index I from the mean temperature
    I = calc_heat_index(last_year_temp_ds['t2m'])
    # if I.all() == 0:
    #     print("I is zero.")
    #     return xr.full_like(temp_ds['t2m'], fill_value=np.nan), temp_ds  # Return NaN if I is zero

    # Calculate the coefficient 'a' from I
    a = calc_sensitivity(I)

    # Calculate PET using the Thornthwaite equation
    # PET = 16 * (10 * temp_ds['t2m'] / I) ** a * (temp_ds['t2m'].time.dt.daysinmonth / 30)
    # temp_ds = temp_ds.sel(expver=1)
    # I = I.sel(expver=1)
    # a = a.sel(expver=1)
    PET = 16 * (10 * temp_ds['t2m'] / I) ** a

    return PET, temp_ds

# def calc_pet_thornthwaite(temp_ds: xr.Dataset) -> tuple[xr.DataArray, xr.DataArray]:
#     """
#     Calculate PET using the Thornthwaite equation.
#     temp: Monthly average 2m temperature in Celsius.
#     """
#     # Select values from "temp_ds" that are negative and set them to 0.
#     nan_mask = temp_ds["t2m"].isnull()
#     temp_ds['t2m'] = temp_ds['t2m'].where(cond=temp_ds['t2m'] >= 0, other=0)
#     temp_ds['t2m'] = temp_ds['t2m'].where(cond=~nan_mask, other=np.nan)

#     # Calculate head index `I`
#     I = calc_heat_index(temp_ds['t2m'])
#     # Calculate the sensitivity of PET to the temperature `a`
#     a = calc_sensitivity(I)
#     # Sun light hours
#     N = 12  # assume 12 hours of sunlight
#     # Days in month
#     NDM = temp_ds['t2m'].time.dt.days_in_month

#     # # Calculate PET
#     # # Assume month's day avg. is 30.437
#     # NDM = 30.437  # average days in a month
#     pet_ds = 16 * ((10 * temp_ds['t2m'] / I) ** a) * (NDM / 30)
#     # pet_ds = 16 * ((10 * temp_ds['t2m'] / I) ** a) * (N/12) * (NDM / 30)

#     # # Modify the calculation of PET based on the Thornthwaite formulation
#     # pet_ds = xr.where(
#     #     temp_ds['t2m'] < 0,
#     #     0,
#     #     pet_ds,
#     # )
#     # pet_ds = xr.where(
#     #     (temp_ds['t2m'] >= 0) & (temp_ds['t2m'] <= 26.5),
#     #     16 * (N/12) * (NDM / 30) * (10 * temp_ds['t2m'] / I) ** a,
#     #     pet_ds,
#     # )
#     # pet_ds = xr.where(
#     #     temp_ds['t2m'] > 26.5,
#     #     N/360 * (-415.85 + 32.24 * temp_ds['t2m'] - 0.43 * temp_ds['t2m'] ** 2),
#     #     pet_ds,
#     # )

#     return pet_ds, temp_ds


# def calc_heat_index(last_year_temp_ds: xr.DataArray) -> xr.DataArray:
#     """
#     Calculate the heat index required for Thornthwaite PET calculation.

#     :param temp_ds: Monthly average 2m temperature in Celsius.
#     """
#     # I = ((temp_ds / 5) ** 1.514).groupby("time.month").mean('time')
#     # I = ((temp_ds / 5) ** 1.514).sum(
#     #     dim=["lat", "lon"],
#     #     skipna=True,
#     #     keep_attrs=True,
#     # )
#     last_year_temp_ds = last_year_temp_ds.where(last_year_temp_ds >= 0, 0)
#     last_year_temp_ds.to_netcdf("last_year_temp_ds.nc")
#     print("* * * * * * * * * * * * ")
#     print("last_year_temp_ds: ", last_year_temp_ds)
#     print("* * * * * * * * * * * * ")
#     I = ((last_year_temp_ds / 5) ** 1.514).sum(dim="time")
#     return I
def calc_heat_index(last_year_temp_ds: xr.DataArray) -> xr.DataArray:
    """
    Calculate the heat index required for Thornthwaite PET calculation.

    :param temp_ds: Monthly average 2m temperature in Celsius.
    """
    # Ensure no negative values; replace them with zero
    last_year_temp_ds = last_year_temp_ds.where(last_year_temp_ds >= 0, 0)
    last_year_temp_ds.to_netcdf("last_year_temp_ds.nc")

    # Calculate the heat index component
    I = ((last_year_temp_ds.sum(dim="time") / 5) ** 1.514)

    # Compute the final result as well to get the actual values
    I = I.compute()

    return I


def calc_sensitivity(I: xr.DataArray) -> xr.DataArray:
    a = (6.75e-07 * I**3) - (7.71e-05 * I**2) + (1.792e-02 * I) + 0.49239
    return a


def preprocess_pet(pet_ds: xr.DataArray, prec_ds: xr.DataArray):
    # Processing after calculation
    if 'expver' in pet_ds.dims:
        pet_ds = pet_ds.isel(expver=1)
    # pet_ds = pet_ds.transpose('time', 'lat', 'lon', 'month')
    pet_ds = pet_ds.transpose('time', 'lat', 'lon')

    # Convert time dimensions type to datetime64[ns]
    pet_ds['time'] = pet_ds['time'].astype('datetime64[ns]')

    pet = pet_ds.interp(lat=prec_ds['precipitation'].lat,
                        lon=prec_ds['precipitation'].lon,
                        method='linear')

    return pet


def calc_difference(prec_ds: xr.DataArray, pet_ds: xr.DataArray, time_scale: int = 1) -> xr.DataArray:
    # if time_scale > 1:
    #     prec_ds = prec_ds.rolling(time=time_scale, center=True).mean()
    return prec_ds['precipitation'] - pet_ds


def calc_spei(D: xr.DataArray, distribution: Literal['gev', 'gl', 'll'] = "ll") -> xr.DataArray:
    """
    D:              water balance (difference between precipitation and PET)
    distribution:   generalized extreme value = 'gev', generalized logistic distribution = 'gl', log-logistic distribution = 'll',
    """

    # Flatten D for fitting and remove NaN and finite values
    D_flat = D.values.flatten()
    D_flat_clean = D_flat[np.isfinite(D_flat)]

    if distribution == 'gev':
        spei: xr.DataArray = gev_distr(D_flat_clean, D)
    elif distribution == 'll':
        spei: xr.DataArray = ll_distr(D_flat_clean, D)
    elif distribution == 'gl':
        spei: xr.DataArray = gl_distr(D_flat_clean)
    else:
        raise ValueError(
            "Unsupported distribution. Choose 'gev', 'gl' or 'll'.")

    return spei


def gev_distr(D_flat_clean: np.ndarray, D_origin: xr.DataArray) -> xr.DataArray:
    # Fit GEV
    shape, loc, scale = gev.fit(D_flat_clean)

    # Calculate standardized value z
    z = (D_flat_clean - loc) / scale

    # Calculate SPEI
    condition = shape != 0
    x = (1 + shape * z) ** (-1/shape)
    y = np.exp(-z)
    t = np.where(condition, x, y)

    spei_values = -np.log(t)

    # Reshape the spei_values to match the original data shape
    spei_values_reshaped = spei_values.reshape(D_origin.shape)

    # Create a new xr.DataArray with the same coordinates as the original data
    spei = xr.DataArray(spei_values_reshaped,
                        coords=D_origin.coords, dims=D_origin.dims)

    return spei


def ll_distr(D_flat_clean: np.ndarray, D_origin: xr.DataArray) -> xr.DataArray:
    """
    Calculate SPEI using the Log-Logistic distribution fitted to the cleaned, flattened data.
    It reconstructs the full-sized data array based on the original's shape and coordinates.
    """
    # Fit the distribution to the cleaned data
    fitted_params = fisk.fit(D_flat_clean)

    # Compute the CDF using fitted parameters and then convert to z-scores
    cdf_values = fisk.cdf(D_flat_clean, *fitted_params)
    spei_values = norm.ppf(cdf_values)

    # Initialize a full-sized array of NaNs based on the shape of D_origin
    spei_full = np.full(D_origin.shape, np.nan, dtype=np.float64)

    # Create a mask of finite values from the original DataArray
    finite_mask = np.isfinite(D_origin.values.flatten())

    # Ensure the flattened data structures are used for assignments
    spei_full_flattened = spei_full.flatten()
    spei_full_flattened[finite_mask] = spei_values
    spei_full = spei_full_flattened.reshape(D_origin.shape)

    # Create an xarray.DataArray from the reshaped array
    return xr.DataArray(spei_full, coords=D_origin.coords, dims=D_origin.dims)


def gl_distr(D_flat_clean: np.ndarray):
    # Fit the generalized logistic distribution to the cleaned data
    params = genlogistic.fit(D_flat_clean)

    # Apply the CDF to the original data using the fitted parameters
    # Note: We need to apply the CDF to each element in the original data array.
    # This involves creating a function that can be applied over the DataArray.
    cdf_applied = np.vectorize(lambda x: genlogistic.cdf(x, *params))

    # Since we can't directly apply np.vectorize result on xarray DataArray, we use .map() method of xarray if available,
    # or apply the function to the .values of the DataArray and then create a new DataArray from the result.
    # Ensure the output has the same shape as the original DataArray.
    cdf_values = xr.DataArray(
        cdf_applied(D_flat_clean.values),
        dims=D_flat_clean.dims,
        coords=D_flat_clean.coords,
    )

    spei_da = xr.DataArray(
        ndtri(cdf_values),
        coords=D_flat_clean.coords,
        dims=D_flat_clean.dims,
    )

    return spei_da


def calc_cdf_values(D_ds: xr.DataArray) -> xr.DataArray:
    # Convert the DataArray to a NumPy array and flatten it
    D_vals_ds: np.ndarray = D_ds.values
    D_vals_ds = D_vals_ds.flatten()
    # Remove NaN or infinite values to clean the data
    D_clean_ds = D_vals_ds[np.isfinite(D_vals_ds)]
    # Fit the generalized logistic distribution to the cleaned data
    params = genlogistic.fit(D_clean_ds)

    # Apply the CDF to the original data using the fitted parameters
    # Note: We need to apply the CDF to each element in the original data array.
    # This involves creating a function that can be applied over the DataArray.
    cdf_applied = np.vectorize(lambda x: genlogistic.cdf(x, *params))

    # Since we can't directly apply np.vectorize result on xarray DataArray, we use .map() method of xarray if available,
    # or apply the function to the .values of the DataArray and then create a new DataArray from the result.
    # Ensure the output has the same shape as the original DataArray.
    cdf_values_ds = xr.DataArray(cdf_applied(
        D_ds.values), dims=D_ds.dims, coords=D_ds.coords)

    return cdf_values_ds


# def create_spei_ds(D_ds, cdf_values_ds):
#     return xr.DataArray(ndtri(cdf_values_ds), coords=D_ds.coords, dims=D_ds.dims)


def spei_plot(spei: xr.DataArray | np.ndarray, shape_file_path, lon_bounds, lat_bounds,
              title, save_path="", show=True):
    print(spei.ndim)
    print(spei.shape)
    # Access the underlying numpy array for plotting
    if isinstance(spei, xr.DataArray):
        spei = spei.values
        print(spei.shape)

    spei = spei.squeeze()

    # If the data_array still isn't 2D because it contains non-singleton dimensions that weren't squeezed,
    # you may need to select a specific slice. For example, if the first and last dimensions are time and channel:
    # data_2d = data_array.isel(time=0, channel=0)
    print(spei.ndim)
    print(spei.shape)

    # Ensure the data is 2D
    if spei.ndim != 2:
        if spei.ndim > 2:
            # spei = spei[-1, :, :]
            # spei = np.nansum(spei, axis=0)  # along the `time` dimension
            # Initialize an output array filled with NaNs of shape (lat, lon)
            spei_sum = np.full(spei.shape[1:], np.nan)

            # Iterate over each lat-lon position
            for lat in range(spei.shape[1]):
                for lon in range(spei.shape[2]):
                    # Select the time series for this lat-lon position
                    time_series = spei[:, lat, lon]
                    # Ignore NaN inf and -inf values
                    time_series = time_series[np.isfinite(time_series)]

                    # Check if all values are NaN
                    if np.all(np.isnan(time_series)):
                        spei_sum[lat, lon] = np.nan
                        # print("1: ", spei_sum[lat, lon])
                        if lat >= 7.0 and lat <= 8.0 and lon >= 47.0 and lon <= 48.0:
                            print(
                                f"1 - Lat:{lat},Lon:{lon} -> {spei_sum[lat, lon]}")
                    else:
                        # Use np.nansum to ignore NaNs if there's at least one non-NaN
                        spei_sum[lat, lon] = np.nansum(time_series)
                        # print("2: ", spei_sum[lat, lon])
                        if lat >= 7.0 and lat <= 8.0 and lon >= 47.0 and lon <= 48.0:
                            print(
                                f"2 - Lat:{lat},Lon:{lon} -> {spei_sum[lat, lon]}")
            spei = spei_sum
        else:
            mean_arr = np.mean(spei, axis=0)
            final_mean = np.mean(mean_arr, axis=-1)
            spei = final_mean

    assert spei.ndim == 2, "Data is not 2D after squeezing/calculation mean."

    colors = ['#8B1A1A', '#DE2929', '#F3641D', '#FDC404',
              '#9AFA94', '#03F2FD', '#12ADF3', '#1771DE', '#00008B']
    spei_cmap = LinearSegmentedColormap.from_list(
        'spei_cmap', colors, N=1000)
    spei_norm = Normalize(vmin=-2.0, vmax=2.0)  # Normalize from -2.0 to 2.0

    # Assuming 'data_2d' is your 2D data array and 'spei_cmap' is your colormap
    # Load the shapefile
    region_shp = gpd.read_file(shape_file_path)
    region_shp = region_shp.to_crs('EPSG:4326')

    # Plot the data
    plt.figure(figsize=(10, 6))
    plt.imshow(spei, extent=[lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]],
               origin='lower',
               cmap=spei_cmap,
               norm=spei_norm,
               )

    # Overlay the shapefile
    region_shp.boundary.plot(ax=plt.gca(), linewidth=0.5, edgecolor='black')

    # Add a colorbar and labels
    plt.colorbar(label='SPEI')
    plt.title(title)
    plt.xlabel('Longitude (Degrees East)')
    plt.ylabel('Latitude (Degrees North)')

    # Save plot
    if save_path != "":
        plt.savefig(save_path)

    # Show plot
    if show:
        plt.show()
