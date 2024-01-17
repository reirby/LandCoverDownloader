# LandCoverDownloader

## Description
LandCoverDownloader is a Python script designed to streamline the download of [NLCD 2021 land cover](https://www.mrlc.gov/data/nlcd-2021-land-cover-conus) data for the contiguous US. It uses the extent of a provided shape file to generate a TIF file in epsg:4326.   
The parameters are:   
  -a path to the target shape file    
  -target spatial resolution in epsg:4326 degrees (optional) 

## Installation
To use LandCoverDownloader, ensure you have GDAL installed and available in python. If not, you can try installing the necessary packages using:
```
pip install -r requirements.txt  # Install dependencies
```
Or create an Anaconda environment with the provided [GDAL_python.yml](GDAL_python.yml) file.

## Usage
Run the script by providing the path to the target shape file and the optional target spatial resolution:
```
python get_NLCD_2021.py D:/path/to/the/shapefile/my_area.shp 0.002
```
This will produce a my_area.tif downloaded for the extent of the input shape file, with 0.002 degrees cell size, saved to the same folder as the input shape file.

## License
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.