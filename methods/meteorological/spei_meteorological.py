import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
from datetime import datetime, timedelta
from pyproj import Transformer
from scipy import stats
from os.path import exists
from osgeo import osr


def read_asc_file(filepath):
    """Reads an .asc file with a variable header size and returns the data as a numpy array."""
    header_keys = ['NCOLS', 'NROWS', 'XLLCORNER',
                   'YLLCORNER', 'CELLSIZE', 'NODATA_VALUE']
    header = {}

    with open(filepath, 'r') as file:
        line = file.readline()
        while line:
            parts = line.strip().split()
            if parts[0].upper() in header_keys:
                header[parts[0].lower()] = float(parts[1]) if parts[0].upper(
                ) != 'NCOLS' and parts[0].upper() != 'NROWS' else int(parts[1])
                # Break the loop if all necessary header info is collected
                if len(header) == len(header_keys):
                    break
            line = file.readline()

        # No need to skip rows now, as the file pointer is already at the right position
        data = np.loadtxt(file, skiprows=0)
        # Convert NoData values to NaN
        data[data == header['nodata_value']] = np.nan

    return data, header


def write_asc_file(data, file_name, header):
    """Writes data to an .asc file using the provided header information."""

    try:
        with open(file_name, 'w') as file:
            file.write(f"NCOLS {header['ncols']}\n")
            file.write(f"NROWS {header['nrows']}\n")
            file.write(f"XLLCORNER {header['xllcorner']}\n")
            file.write(f"YLLCORNER {header['yllcorner']}\n")
            file.write(f"CELLSIZE {header['cellsize']}\n")
            file.write(f"NODATA_VALUE {header['nodata_value']}\n")
            np.savetxt(file, data, fmt="%.4f")
            return True
    except Exception as e:
        print(f"Failed to save file '{file_name}'. Exception: '{e}'")
        return False


def estimate_pet_tw(temperature):
    """Estimates PET using the Thornthwaite equation."""
    # This is a simplified example; actual implementation may vary
    # depending on your exact needs and temperature input
    I = (temperature / 5.0) ** 1.514
    a = (6.75e-7 * I ** 3) - (7.71e-5 * I ** 2) + (1.792e-2 * I) + 0.49239
    PET = 16 * ((10 * temperature / I) ** a)
    return PET


def calculate_spei(D, header, distribution="ll"):
    """Calculates SPEI using a log-logistic distribution"""
    if distribution == "ll":
        # Log-logistic distribution
        params = stats.genlogistic.fit(D[~np.isnan(D)])
        # Calculate the CDF values of the fitted distribution
        cdf = stats.genlogistic.cdf(D, *params)
    elif distribution == "g":
        # Gamma distribution
        params = stats.gamma.fit(D[~np.isnan(D)])
        # Calculate the CDF values of the fitted distribution
        cdf = stats.gamma.cdf(D, *params)
    else:
        print(f"ERROR: Distribution '{distribution}' is not supported.")
        exit(0)

    # Standardize the CDF values using the inverse of the standard normal distribution (z-scores)
    spei = stats.norm.ppf(cdf)
    # Replace NaN and inf values with 'nodata_value'
    spei = np.where(np.isnan(spei) | np.isinf(spei),
                    header['nodata_value'], spei)
    return spei


def estimate_svp(T):
    """Calculate the estimated saturation vapor pressure (in millibars) for a given temperature (in Celsius) using the August-Roche-Magnus formula."""
    return 6.1094 * np.exp((17.625 * T) / (T + 243.04))
    # return 6.112 * np.exp((17.67 * T) / (T + 243.5))
    # return 0.61078 * np.exp(17.269388 * (T) / (T - 35.86))
    # """Calculate the estimated saturation vapor pressure (in millibars) for a given temperature (in Celsius) using the Buck formula. (see https://en.wikipedia.org/wiki/Vapour_pressure_of_water)"""
    # return 6.1121 * np.exp((18.678 - (T / 234.5)) * (T / 257.14 + T))


def estimate_net_radiation(Rs_kwh, T_max, T_min, svp):
    """
    Estimate net radiation using the daily solar radiation, the temperature and the saturation vapor pressure .

    Parameters:
    - Rs_kwh: Daily solar radiation (kWh/m^2)
    - T_max, T_min: Daily maximum and minimum temperatures (degrees Celsius)
    - svp: Saturation vapor pressure (millibars)

    Returns:
    - Rn: Net radiation (MJ/m^2/day)
    """
    # Constants
    sigma = 4.903e-9  # Stefan-Boltzmann constant (MJ K⁻⁴ m⁻² day⁻¹)
    albedo = 0.23  # Reflectivity of the surface, common value for grass

    # Convert solar radiation from kWh/m² to MJ/m⁻²/day⁻¹
    Rs = Rs_kwh * 3.6  # TODO: use "days_in_month=30.4583" ?

    # Convert saturation vapor pressure from millibars to kPa
    es_kpa = svp * 0.1

    # Estimate net shortwave radiation
    Rns = (1 - albedo) * Rs

    # Estimate net longwave radiation using saturation vapor pressure as a proxy
    T_max_K = T_max + 273.15  # Convert to Kelvin
    T_min_K = T_min + 273.15  # Convert to Kelvin
    Rnl = sigma * ((T_max_K**4 + T_min_K**4) / 2) * (0.34 - 0.14 *
                                                     np.sqrt(es_kpa)) * (1.35 * (Rs / (Rs_kwh * 3.6)) - 0.35)

    # Net radiation
    Rn = Rns - Rnl

    return Rn


def estimate_dew_point(T_min, T_max):
    """
    Enhanced approximation of the dew point using the maximum and minimum temperatures.
    This method assumes that the dew point is more accurately estimated by considering
    the daily temperature range.

    Parameters:
    - T_max: Maximum temperature in Celsius
    - T_min: Minimum temperature in Celsius
    """
    # Calculate the average daily range
    daily_range = T_max - T_min
    # Approximation coefficent
    adjustment_factor = 0.2  # TODO: get real value
    # Calculate the enhanced dew point estimate
    T_dew = T_min + adjustment_factor * daily_range
    return T_dew


def estimate_relative_humidity(T_min, T_max, T_mean):
    """
    Calculate the relative humidity using the Magnus formula and an approximation for the dew point.

    Parameters:
    - T_max: Maximum temperature in Celsius
    - T_mean: Mean temperature in Celsius
    - T_min: Minimum temperature in Celsius (used to approximate dew point)
    """
    # Approximate the dew point from the minimum temperature
    T_dew = estimate_dew_point(T_min, T_max)

    # Calculate saturation vapor pressure at the dew point and at the mean temperature
    E_s_dew = estimate_svp(T_dew)  # SVP at the dew point temperature
    E_s_mean = estimate_svp(T_mean)  # SVP at the mean air temperature

    # Calculate and return the relative humidity as a percentage
    RH = (E_s_dew / E_s_mean) * 100
    return RH


def estimate_pet_fao56pm(Rn, T_mean, RH, es_mb, u2=2.0, G=0, days_in_month=30.4583, nodata_val=-999.0):
    """
    Estimate PET using the FAO-56 Penman-Monteith equation.

    Parameters:
    - Rn: Net radiation (MJ/m^2/day)
    - T_mean: Mean air temperature (°C)
    - RH: Relative Humidity (%)
    - es_mb: Saturation vapor pressure (millibars)
    - u2: Wind speed at 2m height (m/s), assuming an average if not provided

    Returns:
    - PET: Potential Evapotranspiration (mm/month)
    """
    # Constants
    # albedo = 0.23
    # sigma = 4.903e-9  # Stefan-Boltzmann constant (MJ/K^4/m^2/day)
    es_kpa = es_mb * 0.1  # Convert es from millibars to kPa
    ea_kpa = es_kpa * (RH / 100.0)  # Calculate actual vapor pressure

    # Calculate other required parameters
    delta_val = (4098 * (0.6108 * np.exp((17.27 * T_mean) / (T_mean + 237.3)))
                 ) / ((T_mean + 237.3) ** 2)  # Slope of vapor pressure curve
    gamma_val = 0.665 * 0.001 * 101.3  # Psychrometric constant, simplified assumption

    # FAO-56 Penman-Monteith equation
    PET = (0.408 * delta_val * (Rn - G) + gamma_val * (900 / (T_mean + 273))
           * u2 * (es_kpa - ea_kpa)) / (delta_val + gamma_val * (1 + 0.34 * u2))

    # Replace nan with NODATA_VALUE
    PET = np.where(np.isnan(PET), nodata_val, PET)

    return PET


def estimate_pet_pt(Rn, T_mean, G=0, gamma_val=0.066, alpha_val=1.26, days_in_month=30.4583, nodata_val=-999.0):
    """
    Estimate PET using the Priestley-Taylor equation.

    Parameters:
    - Rn: Net radiation (MJ/m^2/day)
    - T_mean: Mean air temperature (°C)
    - days_in_month: Total days count of the to be analyzed month (average 30.4583)
    - G: Ground heat flux (MJ/m⁻²/day⁻¹) (usually 0)
    - gamma_val: Psychrometric constant (usually ~0.066)
    - alpha_val: Empirical constant factor (usually 1.26)

    Returns:
    - PET: Potential Evapotranspiration (mm/month) 
    """

    # Delta: Slope of vapor pressure curve
    delta_val = (4098 * (0.6108 * np.exp((17.27 * T_mean) / (T_mean + 237.3)))) / \
        ((T_mean + 237.3) ** 2)

    # PET = alpha_val * (delta_val * (Rn - G) /
    #                    lambda_val * (delta_val + gamma_val))
    PET = alpha_val * (Rn - G) * (delta_val / delta_val + gamma_val)

    # Replace nan with NODATA_VALUE
    PET = np.where(np.isnan(PET), nodata_val, PET)

    return PET


def coords_to_latlon(x, y):
    """
    Convert coordinates from EPSG:31467 to latitude and longitude.
    """
    source = osr.SpatialReference()
    source.ImportFromEPSG(31467)
    target = osr.SpatialReference()
    target.ImportFromEPSG(4326)
    transform = osr.CoordinateTransformation(source, target)
    lon, lat, _ = transform.TransformPoint(x, y)
    return lat, lon


def calc_ra_for_lat(latitude, month):
    """
    Calculate monthly average extraterrestrial radiation (Ra) in MJ/m2/day for a given latitude and month.

    :param latitude: Latitude in degrees.
    :param month: Month as integer (1 = January, ..., 12 = December).
    :return: Monthly average Ra in MJ/m2/day.
    """
    # Constants
    Gsc = 0.0820  # Solar constant, MJ/m^2/min
    dr = [1.000110, 1.034420, 1.068710, 1.102990, 1.037390, 1.000880,
          0.968290, 0.939180, 0.921490, 0.896320, 0.868010, 0.833390]  # Earth-Sun distance reciprocals for each month

    # Convert latitude to radians
    phi = np.radians(latitude)

    # Calculate solar declination
    delta = 0.409 * \
        np.sin(((2 * np.pi / 365) * ((month - 1) * 30.4 + 10)) - 1.39)

    # Calculate sunset hour angle
    omega = np.arccos(-np.tan(phi) * np.tan(delta))

    # Calculate daily extraterrestrial radiation for the month
    Ra = ((24 * 60) / np.pi) * Gsc * dr[month - 1] * (omega * np.sin(phi) * np.sin(delta) +
                                                      np.cos(phi) * np.cos(delta) * np.sin(omega))

    return Ra


def calc_ra(header, month):
    """
    Calculate monthly average Ra for each cell in a grid based on header information and month.
    """
    ra_array = np.zeros((header['nrows'], header['ncols']))

    for row in range(header['nrows']):
        for col in range(header['ncols']):
            # Calculate the center of the cell
            x = header['xllcorner'] + (col + 0.5) * header['cellsize']
            y = header['yllcorner'] + (row + 0.5) * header['cellsize']

            # Convert cell center to latitude and longitude
            lat, _ = coords_to_latlon(x, y)

            # Calculate Ra for the cell based on latitude and month
            ra_array[row, col] = calc_ra_for_lat(lat, month)

    return ra_array


def estimate_pet_hg(T_min, T_max, T_mean, header, month):
    """
    Estimate PET using the Hargreaves equation.
    """
    Ra = calc_ra(header, month)
    PET = 0.0023 * Ra * (T_mean + 17.8) * (T_max - T_min)**0.5
    return PET


def calc_spei_n_months(start_date: datetime, end_date: datetime, pet_equation: str = "fao56pm", distribution="ll", calc_1_month=True):
    months_diff = (end_date.year - start_date.year) * 12 + \
        end_date.month - start_date.month + 1
    print(f"Calculating {months_diff}-month SPEI...\nStart date: '{start_date.strftime('%Y%m')}'\nEnd date: '{end_date.strftime('%Y%m')}'\nPET equation: '{pet_equation}'\nDistribution: '{distribution}'\nCalculate 1-month SPEI: {calc_1_month}\n\n")

    T_min_root_path = "/media/jtrvz/1tb/drought_data/temperature/dwd/min"
    T_min_base_file_name = "grids_germany_monthly_air_temp_min_{date}.asc"

    T_max_root_path = "/media/jtrvz/1tb/drought_data/temperature/dwd/max"
    T_max_base_file_name = "grids_germany_monthly_air_temp_max_{date}.asc"

    T_mean_root_path = "/media/jtrvz/1tb/drought_data/temperature/dwd/avg"
    T_mean_base_file_name = "grids_germany_monthly_air_temp_mean_{date}.asc"

    prec_root_path = "/media/jtrvz/1tb/drought_data/precipitation/dwd/avg"
    prec_base_file_name = "grids_germany_monthly_precipitation_{date}.asc"

    Rs_kWh_root_path = "/media/jtrvz/1tb/drought_data/radiation/dwd/global"
    Rs_kWh_base_file_name = "grids_germany_monthly_radiation_global_{date}.asc"

    T_min_asc_content, header = load_multi_asc_files(
        start_date, end_date, T_min_root_path, T_min_base_file_name)
    T_max_asc_content, _ = load_multi_asc_files(
        start_date, end_date, T_max_root_path, T_max_base_file_name)
    T_mean_asc_content, _ = load_multi_asc_files(
        start_date, end_date, T_mean_root_path, T_mean_base_file_name)
    prec_asc_content, _ = load_multi_asc_files(
        start_date, end_date, prec_root_path, prec_base_file_name)
    Rs_kWh_asc_content, _ = load_multi_asc_files(
        start_date, end_date, Rs_kWh_root_path, Rs_kWh_base_file_name)

    if pet_equation.lower() == "fao56pm":
        pet_results = {}
        D_results = {}

        curr_date = start_date
        while curr_date <= end_date:
            curr_date_str = curr_date.strftime('%Y%m')

            RH = estimate_relative_humidity(
                T_min=T_min_asc_content[curr_date_str],
                T_max=T_max_asc_content[curr_date_str],
                T_mean=T_mean_asc_content[curr_date_str],
            )
            svp = estimate_svp(T=T_mean_asc_content[curr_date_str])
            Rn = estimate_net_radiation(
                Rs_kwh=Rs_kWh_asc_content[curr_date_str],
                T_min=T_min_asc_content[curr_date_str],
                T_max=T_max_asc_content[curr_date_str],
                svp=svp,
            )
            pet = estimate_pet_fao56pm(
                Rn=Rn,
                T_mean=T_mean_asc_content[curr_date_str],
                RH=RH,
                es_mb=svp,
            )

            pet_results[curr_date.strftime('%Y%m')] = pet
            D_results[curr_date.strftime(
                '%Y%m')] = prec_asc_content[curr_date_str] - pet

            if months_diff == 1 or calc_1_month:
                # 1-month SPEI
                spei = calculate_spei(
                    D=D_results[curr_date_str], header=header, distribution=distribution)

                output_file_path = f"output/1month/spei_fao56pm_{distribution}_{curr_date_str}.asc"
                ok = write_asc_file(
                    data=spei, file_name=output_file_path, header=header)
                if ok:
                    print(
                        f"SUCCESS: Saved 1-month SPEI as asc file '{output_file_path}'.")
                else:
                    print(
                        f"ERROR: Failed to save 1-month SPEI as asc file '{output_file_path}'.")

                # Plot
                plot_spei(spei, header,
                          f"output/1month/spei_fao56pm_{distribution}_{curr_date_str}.png")

            # Jump to next month
            next_month = curr_date.replace(day=28) + timedelta(days=4)
            curr_date = next_month - timedelta(days=next_month.day - 1)
    elif pet_equation == "hg":
        pet_results = {}
        D_results = {}

        curr_date = start_date
        while curr_date <= end_date:
            curr_date_str = curr_date.strftime('%Y%m')

            pet = estimate_pet_hg(T_min=T_min_asc_content[curr_date_str],
                                  T_max=T_max_asc_content[curr_date_str],
                                  T_mean=T_mean_asc_content[curr_date_str],
                                  header=header,
                                  month=curr_date.month,
                                  )

            pet_results[curr_date.strftime('%Y%m')] = pet
            D_results[curr_date.strftime('%Y%m')] = prec_asc_content[curr_date_str] - pet

            if months_diff == 1 or calc_1_month:
                # 1-month SPEI
                spei = calculate_spei(
                    D=D_results[curr_date_str], header=header, distribution=distribution)

                output_file_path = f"output/1month/spei_{pet_equation}_{distribution}_{curr_date_str}.asc"
                ok = write_asc_file(
                    data=spei, file_name=output_file_path, header=header)
                if ok:
                    print(
                        f"SUCCESS: Saved 1-month SPEI as asc file '{output_file_path}'.")
                else:
                    print(
                        f"ERROR: Failed to save 1-month SPEI as asc file '{output_file_path}'.")

                # Plot
                plot_spei(spei, header,
                          f"output/1month/spei_{pet_equation}_{distribution}_{curr_date_str}.png")

            # Jump to next month
            next_month = curr_date.replace(day=28) + timedelta(days=4)
            curr_date = next_month - timedelta(days=next_month.day - 1)

        # n-month SPEI
        if months_diff > 1:
            # Sum D (water balance) of all months
            D_sum = None
            for date_str in D_results:
                if D_sum is None:
                    D_sum = D_results[date_str]
                else:
                    D_sum += D_results[date_str]

            # Calculate SPEI
            spei = calculate_spei(
                D=D_sum, header=header, distribution=distribution)

            # Save n-month SPEI as asc file
            total_months = len(D_results)
            output_file_path = f"output/{total_months}month/spei_{pet_equation}_{distribution}_{curr_date_str}.asc"
            ok = write_asc_file(
                data=spei,
                file_name=output_file_path,
                header=header,
            )
            if ok:
                print(
                    f"SUCCESS: Saved {total_months}-month SPEI as asc file '{output_file_path}'.")
            else:
                print(
                    f"ERROR: Failed to save {total_months}-month SPEI as asc file '{output_file_path}'.")

            # Save n-month Plot as png file
            plot_spei(spei, header,
                      f"output/{total_months}month/spei_{pet_equation}_{distribution}_{curr_date_str}.png")


def plot_spei(spei_data, header, save_path=""):
    # Mask the NoData values
    spei_data_masked = np.ma.masked_where(spei_data == -999, spei_data)

    colors = ['#8B1A1A', '#DE2929', '#F3641D', '#FDC404',
              '#9AFA94', '#03F2FD', '#12ADF3', '#1771DE', '#00008B']
    spei_cmap = LinearSegmentedColormap.from_list(
        'custom_cmap', colors, N=1000)
    spei_norm = Normalize(vmin=-2.0, vmax=2.0)  # Normalize from -2.0 to 2.0

    plt.figure(figsize=(10, 6))

    # Calculate the extent using the asc header info
    x_min, y_min, x_max, y_max = transform_coordinates(
        header['xllcorner'], header['yllcorner'],
        header['xllcorner'] +
        header['ncols'] * header['cellsize'],
        header['yllcorner'] +
        header['nrows'] * header['cellsize'],
        'EPSG:31467', 'EPSG:4326'
    )
    extent = [x_min, x_max, y_min, y_max]

    # Create plot
    plt.imshow(X=spei_data_masked,
               cmap=spei_cmap,
               norm=spei_norm,
               extent=extent,
               # extent=[0, header['ncols'], 0, header['nrows']]
               )
    plt.colorbar(label='SPEI Value')
    plt.title('SPEI Values for Germany')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    # Save plot
    if save_path != "":
        plt.savefig(save_path)
        print(f"SUCCESS: Saved SPEI plot as png file '{save_path}.")
    else:
        plt.show()


def transform_coordinates(x_min, y_min, x_max, y_max, src_crs, dest_crs):
    transformer = Transformer.from_crs(src_crs, dest_crs, always_xy=True)
    x_min_transformed, y_min_transformed = transformer.transform(x_min, y_min)
    x_max_transformed, y_max_transformed = transformer.transform(x_max, y_max)
    return x_min_transformed, y_min_transformed, x_max_transformed, y_max_transformed


def load_multi_asc_files(start_date: datetime, end_date: datetime, root_path, base_file_name: str):
    header = {}
    asc_file_content = {}
    current_date: datetime = start_date
    while current_date <= end_date:
        file_name = base_file_name.format(date=current_date.strftime('%Y%m'))
        full_path = f"{root_path}/{file_name}".replace("//", "/")
        if exists(full_path):
            if header:
                asc_file_content[current_date.strftime(
                    '%Y%m')], _ = read_asc_file(full_path)
            else:
                asc_file_content[current_date.strftime(
                    '%Y%m')], header = read_asc_file(full_path)

        next_month = current_date.replace(day=28) + timedelta(days=4)
        current_date = next_month - timedelta(days=next_month.day - 1)

    return asc_file_content, header


########################################
start = datetime(2022, 1, 1)
end = datetime(2022, 3, 1)
calc_spei_n_months(start_date=start, end_date=end,
                   pet_equation="hg",distribution="g",
                   calc_1_month=True)
calc_spei_n_months(start_date=start, end_date=end,
                   pet_equation="hg",distribution="ll",
                   calc_1_month=True)
