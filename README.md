# LandCoverDownloader
A script to download NLCD 2021 land cover data (https://www.mrlc.gov/data/nlcd-2021-land-cover-conus) for the extent of the provided shp file.
Takes an shp file and optionally the target spatial resolution in degrees (epsg:4326). The output is a tif in epsg:4326.

The script requires GDAL!

Please make sure you have necessary packages by using
```
pip install -r requirements.txt  # Install dependencies
```
or by creating an Anaconda environment with the provided GDAL_python.yml file

Use example:
python get_NLCD_2021.py D:/path/to/the/shapefile/my_area.shp 0.002

This will produce a my_area.tif downloaded for the extent of the input shape file, with 0.002 degrees cell size, saved to the same folder as the input shape file.
