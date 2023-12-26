# LandCoverDownloader
A script to download NLCD 2021 land cover data (https://www.mrlc.gov/data/nlcd-2021-land-cover-conus) for the extent of the provided shp file.
Takes an shp file and optionally the target spatial resolution in degrees (4326). The output is a tif in epsg:4326.

Use example:
python get_NLCD_2021.py D:/path/to/the/shapefile/my_area.shp 0.002
