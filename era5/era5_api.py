import cdsapi

file_extensions = {
    "grib": ".grib",
    "netcdf.zip": ".netcdf.zip",
    "netcdf": ".nc"
}


def get_data(year: int, month: int, days: range | int, hours: range | int, variable: str, save_dir: str, area: list = None, format: str = 'grib'):
    """
    Retrieves reanalysis ERA5-Land data from the Copernicus Climate Data Store.

    area: sub-region extraction [north, west, south, east] (e.g. [90, -180, -90, 180])

    Examples: 
    
    get_data(year=2024, month=1, days=range(1, 20), hours=range(12, 15), variable="2m_temperature", save_dir="/home/user/Downloads/era5" format="netcdf")
    -> 2m temperature for the date range from 1st to 20th January of 2024 between the hours 12:00 until 15:00 (12:00 AM until 3:00 PM) in the "netcdf" format.

    get_data(year=2005, month=12, days=5, hours=7, variable="skin_temperature", save_dir="/home/user/Downloads/era5", area=[90, -180, -90, 180], format="grib")
    -> Skin temperature for the day 5th December of 2005 at 7:00 o' clock (7:00 AM) from the area north 90°, west -180°, south -90° and east 180° in the "grib" format.

    See https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form for more info.
    """
    if format not in file_extensions.keys():
        print(
            f"Format '{format}' does not exist. Please choose one of these: {file_extensions.keys()}")
        return
    
    year = str(year)
    month = parse_month(month)
    days = parse_days(days)
    hours = parse_hours(hours)

    c = cdsapi.Client()
    request_parameters = {
        'year': year,
        'month': month,
        'variable': variable,
        'format': format,
        
    }
    if area:
        request_parameters["area"] = area,
    if days:
        request_parameters['day'] = days
    if hours:
        request_parameters['time'] = hours

    # Create file name
    filename = create_filename(year, month, days, hours, variable, format)
    print(filename)

    c.retrieve('reanalysis-era5-land',
               request_parameters,
               target=f"{save_dir}/{filename}.{file_extensions[format]}")

def create_filename(year, month, days, hours, variable, format):
    if isinstance(days, int):
        day_string = f"{days:02d}"
    else:
        day_string = f"{days[0]}-{days[-1]}"
    
    if isinstance(hours, int):
        hour_string = f"{hours:02d}00"
    else:
        hour_start = hours[0].replace(":","")
        hour_end = hours[-1].replace(":","")
        hour_string = f"{hour_start}-{hour_end}"

    return f"{variable}_{year}{month}{day_string}_{hour_string}{file_extensions[format]}"


def parse_hours(hours):
    if isinstance(hours, int):
        # Integer
        if hours < 0 or hours > 23:
            raise ValueError("Hour must be between 0 and 23.")
        return [f"{hours:02d}:00"]
    else:
        # Range
        hour_range = list(hours)
        for hour in hour_range:
            if hour < 0 or hour > 23:
                raise ValueError("Hour must be between 0 and 23.")
        return hour_range

def parse_month(month):
    if not isinstance(month, int):
        raise ValueError("Input must be an integer")

    if month < 1 or month > 12:
        raise ValueError("Day must be between 1 and 31")

    return f"{month:02d}"

def parse_days(days):
    if isinstance(days, int):
        # Integer
        if days < 1 and days > 31:
            raise ValueError("Day must be between 1 and 31.")
        return [days]
    else:
        # Range
        day_range = list(days)
        for day in day_range:
            if day < 1 or day > 31:
                raise ValueError("Day must be between 1 and 31.")
        return day_range


# Example usage:
get_data(year=2023, month=1, days=1, hours=0, variable='2m_temperature', save_dir="/home/jtrvz/Documents/drought_data/temperature/era5/api")

