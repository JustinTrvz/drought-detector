import time
import unittest
import requests
from nasa_gpm import nasa_api
import datetime
import random

class TestNasaApi(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Initialize class variables in setUpClass method
        cls.start_date = datetime.datetime(2000, 1, 1)
        cls.end_date = datetime.datetime.now()
        cls.random_dates = [cls.random_date(cls.start_date, cls.end_date) for _ in range(5)]

    @staticmethod
    def random_date(start_date, end_date):
        # Calculate the difference in days between start_date and end_date
        delta = end_date - start_date
        
        # Generate a random number of days within the range
        random_days = random.randint(0, delta.days)
        
        # Add the random number of days to start_date
        random_date = start_date + datetime.timedelta(days=random_days)
        
        return random_date

    def test_valid_tif_url(self):
        format = "tif"
        for date in self.random_dates:
            for region in nasa_api.region_id_map:
                for dataset in nasa_api.dataset_id_map:
                    tif_url = nasa_api.get_url(date=date, region=region, dataset=dataset, format=format)
                    response = requests.get(tif_url)
                    print(f"Date: '{date}', Region: '{region}', Dataset: '{dataset}', Format: '{format}'")
                    self.assertEqual(response.status_code, 200)
                    time.sleep(100)  # wait 1 sec until next request

    def test_valid_geojson_urls(self):
        for date in self.random_dates:
            for region in nasa_api.region_id_map:
                for dataset in nasa_api.dataset_id_map:
                    tif_url = nasa_api.get_url(date=date, region=region, dataset=dataset, format="geojson")
                    response = requests.get(tif_url)
                    self.assertEqual(response.status_code, 200)
            
        
