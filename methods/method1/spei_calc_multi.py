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


def calc_pet_thornthwaite(temp_ds: xr.Dataset, last_year_temp_ds: xr.Dataset, extended_version: bool = False):
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
    last_year_temp_ds['t2m'] = last_year_temp_ds['t2m'].where(
        cond=last_year_temp_ds['t2m'] >= 0, other=0)
    last_year_temp_ds['t2m'] = last_year_temp_ds['t2m'].where(
        cond=~nan_mask, other=np.nan)

    # Calculate the heat index I from the mean temperature
    I = calc_heat_index(last_year_temp_ds['t2m'])

    # Calculate the coefficient 'a' from I
    a = calc_sensitivity(I)

    if extended_version:
        # Calculate monthly correction factor for day length
        L = calc_correction_factor(temp_ds)
        # Calculate PET using the Thornthwaite equation
        PET = 16 * (10 * temp_ds['t2m'] / I) ** a * L
    else:
        PET = 16 * (10 * temp_ds['t2m'] / I) ** a

    # Set name of the xarray.DataArray variable
    PET.name = "pet"

    return PET, temp_ds


def calc_correction_factor(temp_ds):
    # Get month number from data
    months = temp_ds['time'].dt.month

    # Apply calculate_day_length for each latitude and month
    day_length = xr.apply_ufunc(
        calculate_day_length,
        temp_ds['lat'],
        months,
        vectorize=True,
    )

    # Calculate L
    L = (1 / 12) * (day_length / 12) * day_length**2

    return L


def calculate_day_length(lat, month):
    pi = np.pi
    degrees_to_radians = pi / 180.0

    # Approximate day of the year for the middle of each month
    middle_of_month = [15, 45, 74, 105, 135, 162, 198, 228, 258, 288, 319, 344]

    # Calculate solar declination using approximate day of the year
    delta = 23.45 * np.sin(degrees_to_radians *
                           (360 * (284 + np.array(middle_of_month)) / 365))

    # Convert latitude and declination to radians
    lat_rad = lat * degrees_to_radians
    # month-1 to adjust Python's 0-indexing
    delta_rad = delta[month-1] * degrees_to_radians

    # Calculate hour angle at sunrise/sunset
    cos_omega = -np.tan(lat_rad) * np.tan(delta_rad)
    # Ensure the values are within the valid range for acos
    cos_omega = np.clip(cos_omega, -1, 1)
    omega = np.arccos(cos_omega)

    # Calculate day length in hours
    day_length = (2 * omega * 24) / (2 * pi)

    return day_length


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


def calc_difference(prec_ds: xr.DataArray, pet_ds: xr.DataArray) -> xr.DataArray:
    # if time_scale > 1:
    #     prec_ds = prec_ds.rolling(time=time_scale, center=True).mean()
    D = prec_ds['precipitation'] - pet_ds

    # Set name of the xarray.DataArray variable
    # D.name = "diff"

    return D


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

    # Set name of the xarray.DataArray variable
    spei.name = "spei"

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
    # Select the first time slice
    spei_slice = spei.isel(time=0)

    # Define a custom colormap
    colors = ['#8B1A1A', '#DE2929', '#F3641D', '#FDC404',
              '#9AFA94', '#03F2FD', '#12ADF3', '#1771DE', '#00008B']
    spei_cmap = LinearSegmentedColormap.from_list('spei_cmap', colors, N=1000)
    spei_norm = Normalize(vmin=-2.0, vmax=2.0)

    # Load the shapefile
    region_shp = gpd.read_file(shape_file_path)
    region_shp = region_shp.to_crs('EPSG:4326')

    # Set up the plot
    plt.figure(figsize=(10, 6))
    # Use imshow with xarray's plotting interface for correct axis handling
    # You need to specify 'lon' and 'lat' correctly based on your DataArray's coordinate names
    ax = plt.gca()

    lon_range = spei_slice.lon.max() - spei_slice.lon.min()
    lat_range = spei_slice.lat.max() - spei_slice.lat.min()
    aspect_ratio = lon_range / lat_range

    im = ax.imshow(spei_slice.values, cmap=spei_cmap, norm=spei_norm,
                   extent=[spei_slice.lon.min(), spei_slice.lon.max(),
                           spei_slice.lat.min(), spei_slice.lat.max()],
                   origin='lower', interpolation='none', aspect=aspect_ratio)  # 'extent' helps align the axes

    # Overlay the shapefile
    region_shp.boundary.plot(ax=plt.gca(), linewidth=0.5, edgecolor='black')

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, label='SPEI')
    cbar.set_label('SPEI', rotation=270, labelpad=15)

    # Add titles and labels
    plt.title(title)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    # Save plot
    if save_path != "":
        plt.savefig(save_path)

    # Show plot
    if show:
        plt.show()

    # plt.close()
