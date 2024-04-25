import spei_calc_multi as scm
from datetime import datetime
from dateutil.relativedelta import relativedelta
import xarray as xr
import pandas as pd


class SPEICalculator:
    """
    Initialize the SPEI calculator with the necessary parameters.
    Thornthwaite PET method is used for PET calculation.

    Parameters:
    begin_date (datetime): Begin date of the calculation
    end_date (datetime): End date of the calculation
    time_scale (int): n-month time scale (1, 3, 6, 9, 12, ...)
    lat_bounds (list): [min_lat, max_lat]
    lon_bounds (list): [min_lon, max_lon]
    prec_ds_dir (str): Directory with precipitation data
    temp_ds_dir (str): Directory with temperature data
    """
    begin_date = None
    end_date = None
    time_scale = None
    lat_bounds = []
    lon_bounds = []
    prec_ds_dir = ""
    temp_ds_dir = ""
    save_dir = ""
    prec_ds: xr.Dataset = None
    temp_ds: xr.Dataset = None
    diff_ds: xr.Dataset = None
    pet_ds: xr.Dataset = None
    spei_ds: xr.Dataset = None

    def __init__(self, begin_date: datetime, end_date: datetime,
                 time_scale: int, lat_bounds: list, lon_bounds: list,
                 prec_ds_dir: str, temp_ds_dir: str, save_dir: str):
        print("Initializing SPEI calculator...")
        self.begin_date: datetime = begin_date  # Begin date of the calculation
        self.end_date: datetime = end_date  # End date of the calculation
        # n-month time scale (1, 3, 6, 9, 12, ...)
        self.time_scale: int = time_scale
        self.lat_bounds: list = lat_bounds  # [min_lat, max_lat]
        self.lon_bounds: list = lon_bounds  # [min_lon, max_lon]
        self.prec_ds_dir: str = prec_ds_dir  # Directory with precipitation data
        self.temp_ds_dir: str = temp_ds_dir  # Directory with temperature data
        self.save_dir: str = save_dir  # Directory to save the output
        if self.save_dir[-1] != "/":  
            self.save_dir += "/"  # Add "/" to the end of the path
        print("Initialization done.")

    def calculate(self) -> xr.Dataset:
        # Read nc files
        print("Reading nc files...")
        self.prec_ds, self.temp_ds = self.read_nc_files(
            self.begin_date, self.end_date, self.prec_ds_dir, self.temp_ds_dir)
        print("Reading nc files done.")

        # Preprocess datasets
        print("Preprocessing datasets...")
        self.prec_ds, self.temp_ds = self.preprocess_datasets(
            self.prec_ds, self.temp_ds)
        print("Preprocessing datasets done.")

        # Spatial subset datasets
        print("Spatial subsetting datasets...")
        self.prec_ds, self.temp_ds = self.spatial_subset_datasets(
            self.prec_ds, self.temp_ds, self.lat_bounds, self.lon_bounds)
        print("Spatial subsetting datasets done.")

        # If time_scale is greater than 2, and we do not subtract 1 from time_scale
        # then e.g. time_scale = 2, "2014-03" - time_scale = "2014-01"
        # but we want to have two months, not three!
        if self.time_scale >= 2:
            time_scale_temp = self.time_scale - 1
        else:
            time_scale_temp = 0
        
        # Example: 3-month SPEI for 2014-08, window: 2014-06 to 2014-08 needed
        # e.g. begin_date = "2014-01", end_date = "2014-12", time_scale = 3
        #      windows_begin = "2014-01", window_end = "2014-03", target_date = "2014-04"
        print("Set window begin and end dates...")
        window_begin: datetime = self.begin_date
        window_end: datetime = window_begin + \
            relativedelta(months=time_scale_temp)
        print("Set window begin and end dates done.")

        print("Iterating over the time window...")
        counter = 0
        while window_end <= self.end_date:
            print(f"Calculating SPEI for {window_end.strftime('%Y-%m')}")
            print(f"Window: {window_begin.strftime('%Y-%m')} - {window_end.strftime('%Y-%m')}")
            print(f"Begin iteration #{counter}")
            
            # Slice datasets
            window_begin_last_year = window_begin - relativedelta(years=1)
            last_year_temp_ds = self.temp_ds.sel(time=slice(window_begin_last_year.strftime('%Y-%m-%d'), window_begin.strftime('%Y-%m-%d')))
            prec_windowed_ds = self.prec_ds.sel(
                time=slice(window_begin, window_end))
            temp_windowed_ds = self.temp_ds.sel(
                time=slice(window_begin, window_end))
            print("Slicing datasets done.")

            # Calculate PET
            pet_windowed_ds = self.calc_pet_thornthwaite(
                temp_ds=temp_windowed_ds,
                prec_ds=prec_windowed_ds,
                last_year_temp_ds=last_year_temp_ds,
            )
            print("PET calculation done.")
            # Save PET dataset
            save_path = self.save_dir + f"pet_{window_end.strftime('%Y-%m')}.nc"
            pet_windowed_ds.to_netcdf(save_path)
            print(f"Saved PET to '{save_path}'.")
            # Concatenate PET datasets
            if self.pet_ds is None:
                # First iteration
                self.pet_ds = pet_windowed_ds
            else:
                # Concatenate new PET dataset
                self.pet_ds = xr.concat([self.pet_ds, pet_windowed_ds], dim='time')

            # Calculate difference
            diff_windowed_ds = self.calculate_difference(
                prec_ds=prec_windowed_ds,
                pet_ds=pet_windowed_ds,
                time_scale=self.time_scale,
            )
            # Save difference dataset
            save_path = self.save_dir + f"diff_{window_end.strftime('%Y-%m')}.nc"
            diff_windowed_ds.to_netcdf(save_path)
            print(f"Saved D to '{save_path}'.")
            # Concatenate difference datasets
            if self.diff_ds is None:
                # First iteration
                self.diff_ds = diff_windowed_ds
            else:
                # Concatenate new difference dataset
                self.diff_ds = xr.concat([self.diff_ds, diff_windowed_ds], dim='time')

            # Calculate SPEI
            spei_windowed_ds = self.calc_spei(self.diff_ds)
            # Save SPEI dataset
            save_path = self.save_dir + f"spei_{window_end.strftime('%Y-%m')}.nc"
            spei_windowed_ds.to_netcdf(save_path)
            print(f"Saved SPEI to '{save_path}'.")
            # Concatenate SPEI datasets
            if self.spei_ds is None:
                # First iteration
                self.spei_ds = spei_windowed_ds
            else:
                # Concatenate new SPEI dataset
                self.spei_ds = xr.concat([self.spei_ds, spei_windowed_ds], dim='time')

            # Create plot
            shape_file_path = '/home/jtrvz/Downloads/vg2500_geo84/vg2500_krs.shp'
            title = f"{self.time_scale}-month SPEI Germany ({window_end.strftime('%Y-%m')})"
            plot_save_path = self.save_dir + f"spei_{self.time_scale}-month_{window_end.strftime('%Y-%m')}.png"
            scm.spei_plot(
                spei=spei_windowed_ds,
                shape_file_path=shape_file_path,
                lon_bounds=self.lon_bounds,
                lat_bounds=self.lat_bounds,
                title=title,
                save_path=plot_save_path,
                show=False,
            )
            print(f"Saved plot to '{save_path}'.")

            window_begin += relativedelta(months=1)
            window_end += relativedelta(months=1)
            print(f"End iteration #{counter}")
            counter += 1
        print("Iteration over the time window done. Number of iterations: ", counter)

        return self.spei_ds

    def read_nc_files(self, begin_date: datetime, end_date: datetime,
                      prec_ds_dir: str, temp_ds_dir: str):
        prec_file_paths = scm.generate_imerg_filenames(
            begin_date,
            end_date,
            prec_ds_dir,
        )
        print("Generating precipitation file paths done.")
        temp_file_paths = scm.generate_t2m_filenames(
            # last 12 month needed for heat index calculation
            begin_date - relativedelta(months=12),
            end_date,
            temp_ds_dir,
            "t2m",
        )
        print("Generating temperature file paths done.")

        if not prec_file_paths or not temp_file_paths:
            print("No precipitation or temperature files found.")
            raise Exception("No precipitation or temperature files found.")

        prec_ds = scm.read_nc_files(prec_file_paths)
        print("Reading precipitation files done.")
        temp_ds = scm.read_nc_files(temp_file_paths)
        print("Reading temperature files done.")
        return prec_ds, temp_ds

    def preprocess_datasets(self, prec_ds: xr.Dataset, temp_ds: xr.Dataset):
        prec_ds = scm.preprocess_prec(self.prec_ds)
        print("Preprocessing precipitation done.")
        temp_ds = scm.preprocess_temp(self.temp_ds)
        print("Preprocessing temperature done.")
        return prec_ds, temp_ds

    def spatial_subset_datasets(self, prec_ds: xr.Dataset, temp_ds: xr.Dataset,
                                lat_bounds: list, lon_bounds: list):
        prec_ds = scm.spatial_subset(
            prec_ds, lat_bounds, lon_bounds)
        print("Spatial subsetting precipitation done.")
        temp_ds = scm.spatial_subset(
            temp_ds, lat_bounds, lon_bounds)
        print("Spatial subsetting temperature done.")
        return prec_ds, temp_ds

    def calc_pet_thornthwaite(self, temp_ds: xr.Dataset, prec_ds: xr.Dataset,
                              last_year_temp_ds: xr.Dataset):
        pet_ds, temp_ds = scm.calc_pet_thornthwaite(temp_ds, last_year_temp_ds)
        print("Calculating PET Thornthwaite done.")
        pet_ds = scm.preprocess_pet(pet_ds, prec_ds)
        print("Preprocessing PET Thornthwaite done.")
        return pet_ds

    def calculate_difference(self, prec_ds, pet_ds, time_scale):
        diff_ds = scm.calc_difference(
            prec_ds, pet_ds, time_scale)
        print("Calculating difference done.")
        return diff_ds

    def calc_spei(self, diff_ds: xr.Dataset):
        spei_ds = scm.calc_spei(D=diff_ds, distribution="ll")
        print("Calculating SPEI done.")
        return spei_ds
