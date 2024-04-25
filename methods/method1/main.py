from spei_calculator import SPEICalculator
from datetime import datetime

begin = datetime(2014, 1, 1)
end = datetime(2023, 8, 1)
spei_calc = SPEICalculator(
    begin_date=begin,
    end_date=end,
    time_scale=1,
    lat_bounds=[47.0, 55.5],
    lon_bounds=[5.5, 15.5],
    prec_ds_dir="/media/jtrvz/1tb/drought_data/precipitation/nasa_gpm/Global/monthly/netcdf/avg",
    temp_ds_dir="/media/jtrvz/1tb/drought_data/temperature/era5/Global/monthly/netcdf/avg/",
    save_dir="/media/jtrvz/1tb/drought_data/spei/spei_1"
)
spei_ds = spei_calc.calculate()

# import spei_calc_multi as scm
# from datetime import datetime
# from dateutil.relativedelta import relativedelta

# # Define dates
# begin = datetime(2013, 9, 1)
# end = datetime(2023, 8, 31)
# curr = begin
# while curr <= end:
#     print(f"Current date: {curr}")
#     day = 1
#     time_scale = 1

#     date_target = curr
#     date_begin = date_target - relativedelta(months=time_scale)
#     print(date_begin)
#     date_end = date_begin
#     print("date_target: ", date_target)
#     print("date_begin: ", date_begin)
#     print("date_end: ", date_end)

#     # Precipitation
#     prec_file_paths = scm.generate_imerg_filenames(
#         date_begin, date_end, "/media/jtrvz/1tb/drought_data/precipitation/nasa_gpm/Global/monthly/netcdf/avg")
#     # Temperature
#     temp_file_paths = scm.generate_t2m_filenames(
#         date_end - relativedelta(month=12), date_end, "/media/jtrvz/1tb/drought_data/temperature/era5/Global/monthly/netcdf/avg/", "t2m")
#     # last 12 month needed for heat index calculation

#     # if len(prec_file_paths) != len(temp_file_paths):
#     #     print(f"Number of precipitation and temperature files do not match. {len(prec_file_paths)} (prec) != {len(temp_file_paths)} (temp)")
#     #     exit()

#     # Load nc files
#     prec_ds = scm.read_nc_files(prec_file_paths)
#     temp_ds = scm.read_nc_files(temp_file_paths)

#     # Preprocess datasets
#     prec_ds = scm.preprocess_prec(prec_ds)
#     temp_ds = scm.preprocess_temp(temp_ds)

#     # Lat and lon bounds for Germany
#     lat_bounds = [47.0, 55.5]
#     lon_bounds = [5.5, 15.5]

#     # Spatial subset Germany
#     prec_ds = scm.spatial_subset(prec_ds, lat_bounds, lon_bounds)
#     temp_ds = scm.spatial_subset(temp_ds, lat_bounds, lon_bounds)
#     prec_ds.to_netcdf(f"spei_{time_scale}/prec_{date_target.strftime('%Y%m')}.nc")

#     # Calculate PET (potential evapotranspiration)
#     pet_ds, temp_ds = scm.calc_pet_thornthwaite(temp_ds)
#     pet_ds = scm.preprocess_pet(pet_ds, prec_ds)
#     pet_ds.to_netcdf(f"spei_{time_scale}/pet_{date_target.strftime('%Y%m')}.nc")
#     temp_ds.to_netcdf(f"spei_{time_scale}/temp_{date_target.strftime('%Y%m')}.nc")
#     # print(temp_ds)

#     # Calculate difference
#     diff_ds = scm.calc_difference(prec_ds, pet_ds, time_scale)
#     diff_ds.to_netcdf(f"spei_{time_scale}/diff_{date_target.strftime('%Y%m')}.nc")

#     spei = scm.calc_spei(D=diff_ds,
#                     distribution="ll")
#     spei.to_netcdf(f"spei_{time_scale}/spei_{date_target.strftime('%Y%m')}.nc")

#     # Create plot
#     shape_file_path = '/home/jtrvz/Downloads/vg2500_geo84/vg2500_krs.shp'
#     title = f"{time_scale}-month SPEI Germany ({date_target.strftime('%Y-%m')})"

#     scm.spei_plot(
#         spei=spei,
#         shape_file_path=shape_file_path,
#         lon_bounds=lon_bounds,
#         lat_bounds=lat_bounds,
#         title=title,
#         save_path=f"spei_{time_scale}/spei_{date_target.strftime('%Y%m')}.png",
#         show=False)
    
#     curr = curr + relativedelta(months=1)
#     print(f"Increased date by 1 month -> {curr}")
