import json
import os
import urllib
import datetime

import requests

# Examples for TIF download URLs:
# EU https://pmmpublisher.pps.eosdis.nasa.gov/products/s3/r04/gpm_3hr_1d/2024/044/gpm_3hr_1d.20240213.235959.tif
# S. America https://pmmpublisher.pps.eosdis.nasa.gov/products/s3/r03/gpm_3hr_1d/2024/044/gpm_3hr_1d.20240213.235959.tif
# C. America https://pmmpublisher.pps.eosdis.nasa.gov/products/s3/r02/gpm_3hr_1d/2024/044/gpm_3hr_1d.20240213.235959.tif
# N. America https://pmmpublisher.pps.eosdis.nasa.gov/products/s3/r01/gpm_3hr_1d/2024/044/gpm_3hr_1d.20240213.235959.tif
# Africa https://pmmpublisher.pps.eosdis.nasa.gov/products/s3/r05/gpm_3hr_1d/2024/044/gpm_3hr_1d.20240213.235959.tif
# Middle East https://pmmpublisher.pps.eosdis.nasa.gov/products/s3/r06/gpm_1d/2024/044/gpm_1d.20240213.tif
# India https://pmmpublisher.pps.eosdis.nasa.gov/products/s3/r07/gpm_1d/2024/044/gpm_1d.20240213.tif
# Indian Ocean https://pmmpublisher.pps.eosdis.nasa.gov/products/s3/r10/gpm_1d/2024/044/gpm_1d.20240213.tif
# Asia https://pmmpublisher.pps.eosdis.nasa.gov/products/s3/r09/gpm_1d/2024/044/gpm_1d.20240213.tif
# Australia https://pmmpublisher.pps.eosdis.nasa.gov/products/s3/r08/gpm_1d/2024/044/gpm_1d.20240213.tif
region_id_map = {
    "global": "Global",
    "north america": "r01",
    "central america": "r02",
    "south america": "r03",
    "europe": "r04",
    "africa": "r05",
    "middle east": "r06",
    "india": "r07",
    "australia": "r08",
    "asia": "r09",
    "indian ocean": "r10"
}

# Examples for TIF download URLs:
# Early Run 1d (updated every 3 hours) https://pmmpublisher.pps.eosdis.nasa.gov/products/s3/r08/gpm_3hr_1d/2024/044/gpm_3hr_1d.20240213.235959.tif
# Late Run 3d https://pmmpublisher.pps.eosdis.nasa.gov/products/s3/r08/gpm_3d/2024/044/gpm_3d.20240213.tif
# Early Run 3h https://pmmpublisher.pps.eosdis.nasa.gov/products/s3/r08/gpm_3hr/2024/044/gpm_3hr.20240213.235959.tif
# Early Run 30m https://pmmpublisher.pps.eosdis.nasa.gov/products/s3/r08/gpm_30mn/2024/044/gpm_30mn.20240213.235959.tif
dataset_id_map = {
    "late run 1d": "gpm_1d",
    "late run 3d": "gpm_3d",
    "late run 7d": "gpm_7d",
    "early run 30m": "gpm_30mn",
    "early run 3h": "gpm_3hr",
    "early run 1d": "gpm_3hr_1d"
}

data_formats = [
    "geojson",
    "topojson",
    "tif",
    "arcjson"
]

nasa_gpm_dir = "data/precipitation/nasa_gpn"
checkpoint_filename = "checkpoint.json"

MAX_RETRIES = 25


def get_region_id(region: str):
    region_id = region_id_map.get(region, None)

    if region_id is None:
        print(f"Invalid region '{region}' entered.")  # TODO: error logging

    return region_id


def get_dataset_id(dataset: str):
    dataset_id = dataset_id_map.get(dataset, None)

    if dataset_id is None:
        print(f"Invalid dataset '{dataset}' entered.")  # TODO: error logging

    return dataset_id


def get_url(date: datetime.date, region: str, dataset: str, format: str):
    # Sanitize user input
    format = format.lower()
    if not format in data_formats:
        print(f"Output format '{format}' is not supported.")
        return
    region = region.lower()
    dataset = dataset.lower()

    # Region ID
    region_id = get_region_id(region=region)
    if region_id is None:
        print(f"Region id for '{region}' is 'None'.")  # TODO: error logging
        return

    # Dataset ID
    dataset_id = get_dataset_id(dataset=dataset)
    if dataset_id is None:
        print(f"Dataset id for '{region}' is 'None'.")  # TODO: error logging
        return

    # Date parser
    year = date.year
    month = date.strftime("%m")
    day = date.strftime("%d")
    day_of_year = date.timetuple().tm_yday
    # Append leading zero if < 100
    day_of_year = "{:03d}".format(
        day_of_year) if day_of_year < 100 else day_of_year

    # Get specific URL
    if format == "geojson":
        return get_geojson_url(region_id, dataset_id, year, month, day, day_of_year)
    elif format == "tif":
        return get_tif_url(region_id, dataset_id, year, month, day, day_of_year)
    else:
        print(f"Output format '{format}' is not supported yet.")
        return


def get_tif_url(region_id: str, dataset_id: str, year: int, month: str, day: str, day_of_year: str | int):
    return f"https://pmmpublisher.pps.eosdis.nasa.gov/products/s3/{region_id}/{dataset_id}/{year}/{day_of_year}/{dataset_id}.{year}{month}{day}.tif"


def get_geojson_url(region_id: str, dataset_id: str, year: int, month: str, day: str, day_of_year: str | int):
    return f"https://pmmpublisher.pps.eosdis.nasa.gov/products/{dataset_id}/export/{region_id}/{year}/{day_of_year}/{dataset_id}.{year}{month}{day}.geojson"


def save_file(url: str, save_dir: str, forward=False):
    response = requests.get(url)
    file_name = url.split("/")[-1]
    if response.status_code == 200:
        # Create file path
        file_path = f"{save_dir}/{file_name}"

        # Check if file already exists
        if not os.path.exists(file_path):
            # Create dirs if not exist
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
                # TODO :logging
                print(f"Created directory. save_dir='{save_dir}'")

            # Save file
            with open(file_path, "wb") as file:
                file.write(response.content)
            # Save checkpoint
            save_checkpoint(file_name, save_dir, forward)
        else:
            print(f"File already exists. file_path='{file_path}'")
            return

        # Check if file has been saved correctly
        if not os.path.exists(file_path):
            # TODO: logging
            print(f"Failed to write '{file_path}' to file system.")
            return

        print(f"Saved file. file_path='{file_path}'")  # TODO: logging
        return file_path
    else:
        print(f"GET request failed for URL '{response}'.")


def save_files_since_last_checkpoint(region: str, dataset: str, format: str, future: bool = True):
    region_id = get_region_id(region)
    dataset_id = get_dataset_id(dataset)

    dir = f"{nasa_gpm_dir}/{region_id}/{dataset_id}"
    checkpoint_path = f"{dir}/{checkpoint_filename}"

    # Read checkpoint file
    oldest_file, youngest_file = read_checkpoint_file()
    if oldest_file is None or youngest_file is None:
        print(
            f"Checkpoint file does not exist. Please use 'save_from_to()' first. The file will be create automatically afterwards. checkpoint_path='{checkpoint_path}'")
        return

    # Extract the date string
    oldest_file_date_str = oldest_file.split('.')[1]
    youngest_file_date_str = youngest_file.split('.')[1]

    # Convert the date string to a datetime object
    oldest_file_date = datetime.strptime(oldest_file_date_str, '%Y%m%d')
    youngest_file_date = datetime.strptime(youngest_file_date_str, '%Y%m%d')

    if future:
        # All newest files since download session until today's date
        save_files_date_range(region=region, dataset=dataset, format=format,
                              from_date=youngest_file_date, to_date=datetime.datetime.now())
    else:
        # All files from the past going backwards from the first file's date
        save_files_date_range(region=region, dataset=dataset, format=format,
                              to_date=youngest_file_date)



def read_checkpoint_file(checkpoint_path: str):
    # If checkpoint file does not exist
    if not os.path.exists(checkpoint_path):
        print(
            f"Checkpoint file does not exist. checkpoint_path='{checkpoint_path}'")
        return None, None

    # Read checkpoint file
    with open(checkpoint_path, "r") as f:
        checkpoint_data = f.read()

    # Load JSON as dict
    checkpoint_dict = json.loads(checkpoint_data)
    if not checkpoint_dict:
        print(f"Checkpoint file is empty. checkpoint_path='{checkpoint_path}'")
        return None, None

    # Check for keys
    if "oldest_file" in checkpoint_dict and "youngest_file" in checkpoint_dict:
        return checkpoint_dict["oldest_file"], checkpoint_dict["youngest_file"]
    else:
        print(
            f"Checkpoint files does not contain 'oldest_file' and 'youngest_file' key. checkpoint_path='{checkpoint_path}'")
        return


def save_checkpoint(file_name: str, save_dir: str, forward: bool):
    """
    Creates or updates checkpoint file with the oldest file and youngest file name.

    Side note: Old and young refers to the capture time.
    """
    checkpoint_path = f"{save_dir}/{checkpoint_filename}"
    if os.path.exists(checkpoint_path):
        # Read JSON data from file
        with open(checkpoint_path, "r") as file:
            json_data = json.load(file)

        if forward:
            json_data["youngest_file"] = file_name
        else:
            json_data["oldest_file"] = file_name
    else:
        json_data = {
            "oldest_file": file_name,
            "youngest_file": file_name
        }

    # Save checkpoint JSON file
    with open(checkpoint_path, "w") as file:
        json.dump(json_data, file, indent=4)

    print(f"Saved checkpoint JSON file. checkpoint_path='{checkpoint_path}'")
    return

# test


def save_files_date_range(region: str, dataset: str, format: str,
                          to_date: datetime.datetime, from_date: datetime.datetime = None, ):
    """
    Saves files from `from_date` to `to_date`.

    E.g. `from_date=datetime.datetime(2020, 12, 01)`, `to_date= datetime.datetime(2024, 10, 27)`.

    If you provide no `from_date` the `to_date` will counted down until no files are available anymore.

    E.g. Starting with `to_date=datetime.datetime(2024, 10, 27)`, the next save file will be from the 2024/10/26.

    Checkpoint file will be created so you can use the function `save_files_since_last_checkpoint()` next to update your
    file directory with the most recent files.
    """

    fail_counter = 0
    if from_date is None:
        date = to_date
        try:
            while date >= datetime.datetime(2000, 1, 1):
                url = get_url(date=date, region=region,
                              dataset=dataset, format=format)
                ok = save_file(url=url, save_dir=save_dir, forward=False)
                if not ok:
                    print(
                        f"Failed to save file. date='{date}', url='{url}', save_dir='{save_dir}'")

                    # If failed for more than `MAX_RETRIES` times abort
                    fail_counter += 1
                    if fail_counter >= MAX_RETRIES:
                        print(
                            f"Max retries of '{MAX_RETRIES}' reached. date='{date}', to_date='{to_date}', region='{region}', dataset='{dataset}', format='{format}'")
                        return
                else:
                    fail_counter = 0

                # Subtract one day from `to_date`
                date -= datetime.timedelta(days=1)
        except Exception as e:
            save_files_date_range(
                region=region, dataset=dataset, format=format, to_date=date)

    else:
        date = from_date
        while date <= to_date:
            url = get_url(date=date, region=region,
                          dataset=dataset, format=format)
            ok = save_file(url=url, save_dir=save_dir, forward=True)
            if not ok:
                print(
                    f"Failed to save file. url='{url}', save_dir='{save_dir}'")

                # If failed for more than `MAX_RETRIES` times abort
                fail_counter += 1
                if fail_counter >= MAX_RETRIES:
                    print(
                        f"Max retries of '{MAX_RETRIES}' reached. date='{date}', to_date='{to_date}', region='{region}', dataset='{dataset}', format='{format}'")
                    return
            else:
                fail_counter = 0

            date += datetime.timedelta(days=1)  # Add one day to `from_date`

    print(
        f"Saved files from '{date}' to '{to_date}'. region='{region}', dataset='{dataset}', format='{format}'")


date = datetime.datetime(2018, 2, 4)
format = "geojson"
region = "global"
region_id = get_region_id(region)
dataset = "late run 1d"
dataset_id = get_dataset_id(dataset)
save_dir = f"/media/jtrvz/chugchug/Git/drought-detector/data/precipitation/nasa_gpm/{region_id}/{format}/{dataset_id}"
save_files_date_range(region=region, dataset=dataset,
                      format=format, to_date=date)
# save_files_since_last_checkpoint(region=region, dataset=dataset, format=format)

# # Loop downwards from the specified date
# while date >= datetime.datetime(2000, 1, 1):
#     url = get_url(date=date, region=region, dataset=dataset, format=format)
#     ok = save_file(url=url, save_dir=save_dir)
#     if not ok:
#         print(f"Failed to save file. url='{url}', save_dir='{save_dir}'")
#         exit(0)

#     date -= datetime.timedelta(days=1)  # Move to the previous day

# date = datetime.datetime(2018, 12, 31)
# print(get_url(date=date, region="global", dataset="late run 1d", format="geojson"))
# https://pmmpublisher.pps.eosdis.nasa.gov/products/gpm_7d/export/Global/2020/001/gpm_7d.20200101.235959.geojson
# https://pmmpublisher.pps.eosdis.nasa.gov/products/gpm_7d/export/Global/2023/353/gpm_7d.20231219.geojson
