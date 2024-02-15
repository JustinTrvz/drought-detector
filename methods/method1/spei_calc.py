import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
from os.path import exists
from os import remove
from scipy.stats import genlogistic
from scipy.special import ndtri
from datetime import datetime, timedelta
from os.path import exists


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

        # Move to the next month
        # Adding a month to a datetime object is tricky because not all months have the same number of days.
        # This method handles the end-of-month edge cases.
        next_month = current.replace(
            day=28) + timedelta(days=4)  # this will never fail
        current = next_month - timedelta(days=next_month.day - 1)

    return file_names


def generate_t2m_filenames(start_date: datetime, end_date: datetime, dir: str = ""):
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
        file_name = f"{dir}t2m_{current.strftime('%Y%m')}.nc"
        if exists(file_name):
            file_names.append(file_name)
        else:
            print(f"File '{file_name}' does not exist.")

        # Move to the next month
        # Adding a month to a datetime object is tricky because not all months have the same number of days.
        # This method handles the end-of-month edge cases.
        next_month = current.replace(
            day=28) + timedelta(days=4)
        current = next_month - timedelta(days=next_month.day - 1)

    return file_names


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

    # Convert time dimensions type to datetime64[ns]
    temp_ds['time'] = temp_ds['time'].astype('datetime64[ns]')

    # Convert longitude coordinates from 0-360 to -180 to 180
    if (temp_ds.lon >= 180).any():
        print("Longitude conversion...")  # TODO: logging
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


def calc_pet_thornthwaite(temp_ds: xr.Dataset) -> xr.DataArray:
    """
    Calculate PET using the Thornthwaite equation.
    temp: Monthly average temperature in Celsius.
    """
    # Select `t2m` values
    temp_ds = temp_ds['t2m']

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
    pet_ds: xr.DataArray = (16 * (10 * temp_celsius / I)
                            ** a) * (days_in_month / 12)

    return pet_ds


def preprocess_pet(pet_ds: xr.DataArray, prec_ds: xr.DataArray) -> xr.DataArray:
    # Processing after calculation
    pet_ds = pet_ds.sel(expver=1)
    pet_ds = pet_ds.transpose('time', 'lat', 'lon', 'month')

    # Select `precipitation` values
    prec_ds = prec_ds['precipitation']

    # Convert time dimensions type to datetime64[ns]
    pet_ds['time'] = pet_ds['time'].astype('datetime64[ns]')

    pet = pet_ds.interp(lat=prec_ds.lat, lon=prec_ds.lon, method='linear')

    return pet, prec_ds


def calc_difference(prec_ds: xr.DataArray, pet_ds: xr.DataArray) -> xr.DataArray:
    return prec_ds - pet_ds


# def calc_cdf_values2(diff_ds: xr.DataArray) -> xr.DataArray:
#     diff_ds_vals = diff_ds.values
#     diff_ds_flat = diff_ds_vals.flatten()
#     diff_clean = diff_ds_flat[np.isfinite(diff_ds_flat)]
#     params = genlogistic.fit(diff_clean)

#     return genlogistic.cdf(diff_ds, *params)


def calc_cdf_values(diff_ds: xr.DataArray) -> xr.DataArray:
    # Convert the DataArray to a NumPy array and flatten it
    diff_ds_vals: np.ndarray = diff_ds.values
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


def spei_plot(spei_da, shape_file_path, lon_bounds, lat_bounds,
              title, save_path=""):
    # Access the underlying numpy array for plotting
    data = spei_da.values
    print(data.shape)
    data_2d = data.squeeze()

    # If the data_array still isn't 2D because it contains non-singleton dimensions that weren't squeezed,
    # you may need to select a specific slice. For example, if the first and last dimensions are time and channel:
    # data_2d = data_array.isel(time=0, channel=0)

    # Ensure the data is 2D at this point
    print(data_2d.ndim)
    print(data_2d.shape)
    assert data_2d.ndim == 2, "Data is not 2D after squeezing/selecting."

    colors = ['#8B1A1A', '#DE2929', '#F3641D', '#FDC404',
              '#9AFA94', '#03F2FD', '#12ADF3', '#1771DE', '#00008B']
    custom_cmap = LinearSegmentedColormap.from_list(
        'custom_cmap', colors, N=len(colors))

    # Assuming 'data_2d' is your 2D data array and 'custom_cmap' is your colormap
    # Load the shapefile
    region_shp = gpd.read_file(shape_file_path)
    region_shp = region_shp.to_crs('EPSG:4326')

    # Plot the data
    plt.figure(figsize=(10, 6))
    plt.imshow(data_2d, extent=[lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]],
               origin='lower', cmap=custom_cmap, vmin=-3, vmax=3)

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
    plt.show()
