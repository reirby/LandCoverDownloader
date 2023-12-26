"""
Download NLCD_2021_Land_Cover_L48 (US) landcover from www.mrlc.gov WCS 
for the extent of the input SHP file. The output is a tif image in 4326 coordinate system.
The target spatial resolution can be specified in degrees (4326) as an optional parameter.
Images for large extents are downloaded in tiles and then merged together

"""
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from osgeo_utils import gdal_calc
import os
import math
import sys

def main():
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3 and len(sys.argv) != 2:
        print("Usage: python get_NLCD_2021.py shp-file  cell-size-in-degrees(optional)")
        sys.exit(1)  # Exit with a non-zero status to indicate an error
   
   # Extract command-line arguments
    shapefile_path = sys.argv[1]
    cell_size = sys.argv[2]    
   
    # Check if the second argument is a valid float
    try:
        cell_size = float(cell_size)
    except ValueError:
        print("Error: The second parameter must be a valid float.")
        sys.exit(1)


    get_LC(shapefile_path, cell_size)



def get_LC(shapefile_path, spatial_resolution=None):
    min_x, max_x, min_y, max_y = get_bbox_4326(shapefile_path)
    spanX = abs(max_x - min_x)
    spanY = abs(max_y - min_y)
    area = spanX*spanY
    folder_path = os.path.dirname(shapefile_path)
    file_name = os.path.basename(shapefile_path)
    
    if area <= 9:
        #just get the whole image
        outImg = os.path.join(folder_path, file_name[:-4]+'.tif')
        wcs_url = f'https://www.mrlc.gov/geoserver/mrlc_display/NLCD_2021_Land_Cover_L48/wcs?service=WCS&version=2.0.1&request=getcoverage&coverageid=NLCD_2021_Land_Cover_L48&subset=Lat({min_y},{max_y})&subset=Long({min_x},{max_x})&SubsettingCRS=http://www.opengis.net/def/crs/EPSG/0/4326'
        get_IMG(wcs_url,outImg,spatial_resolution)
    else:
        #split extent to tiles of ~ 1sq deg
        tileExtents = split_extent(min_x, max_x, min_y, max_y)
        failedImg = []
        tilePaths = []
        for i, ext in enumerate(tileExtents):
            print(f'Getting image {i+1} out of {len(tileExtents)}')
            min_x, max_x, min_y, max_y = ext
            if spatial_resolution:
                min_x, max_x, min_y, max_y = min_x - spatial_resolution, max_x+spatial_resolution, min_y-spatial_resolution, max_y+spatial_resolution
                
            wcs_url = f'https://www.mrlc.gov/geoserver/mrlc_display/NLCD_2021_Land_Cover_L48/wcs?service=WCS&version=2.0.1&request=getcoverage&coverageid=NLCD_2021_Land_Cover_L48&subset=Lat({min_y},{max_y})&subset=Long({min_x},{max_x})&SubsettingCRS=http://www.opengis.net/def/crs/EPSG/0/4326'
            outFolder = os.path.join(folder_path, file_name[:-4])
            if not os.path.exists(outFolder):
                os.makedirs(outFolder)
            outImg = os.path.join(outFolder, file_name[:-4]+f'_{i}.tif')
            result = get_IMG(wcs_url, outImg, spatial_resolution)
            if result == 'fail':
                failedImg.append(i)
            else:
                tilePaths.append(outImg)
        #try one more time
        if len(failedImg) > 0:
            for i in failedImg:
                min_x, max_x, min_y, max_y = tileExtents[i]
                if spatial_resolution:
                    min_x, max_x, min_y, max_y = min_x - spatial_resolution, max_x+spatial_resolution, min_y-spatial_resolution, max_y+spatial_resolution
                outImg = os.path.join(outFolder, file_name[:-4]+f'_{i}.tif')
                wcs_url = f'https://www.mrlc.gov/geoserver/mrlc_display/NLCD_2021_Land_Cover_L48/wcs?service=WCS&version=2.0.1&request=getcoverage&coverageid=NLCD_2021_Land_Cover_L48&subset=Lat({min_y},{max_y})&subset=Long({min_x},{max_x})&SubsettingCRS=http://www.opengis.net/def/crs/EPSG/0/4326'
                result = get_IMG(wcs_url, outImg,spatial_resolution)
                if result == 'success':
                    tilePaths.append(outImg)
                    
        #mosaic tiles back together
        output_mosaic = os.path.join(outFolder, file_name[:-4]+'_mosaic.tif')
        mosaic_tifs(tilePaths, output_mosaic)
            
def split_extent(min_x, max_x, min_y, max_y):
    # Size of each subextent (approximately 1 square degree)
    subextent_size = 1.0
    tileBuffer = 0.0002
    subextents = []

    # Calculate the number of rows and columns
    num_rows = math.ceil((max_y - min_y) / subextent_size)
    num_cols = math.ceil((max_x - min_x) / subextent_size)

    stepX = (max_x - min_x)/num_cols
    stepY = round((max_y - min_y)/num_rows,8)
    
    for row in range(num_rows):
        for col in range(num_cols):
            # Calculate the bounds of the subextent
            subextent_min_x = round((min_x + col * stepX) - tileBuffer,8)
            subextent_max_x = round((min_x + (col + 1) * stepX) + tileBuffer,8)
            subextent_min_y = round((min_y + row * stepY) - tileBuffer,8)
            subextent_max_y = round((min_y + (row + 1) * stepY) + tileBuffer,8)

            # Add the subextent to the list
            subextents.append((subextent_min_x, subextent_max_x, subextent_min_y, subextent_max_y))

    return subextents

def get_IMG(wcs_url, fname, spatial_resolution=None):

    # Open the WCS dataset
    ds = gdal.Open(wcs_url)

    # Check if the dataset is successfully opened
    if ds is None:
        print('Failed to open WMS dataset')
        return 'fail'
    else:
        # Define the output file format
        output_format = 'GTiff'

        # Create output dataset
        driver = gdal.GetDriverByName(output_format)

        if spatial_resolution:
            # Use specified spatial resolution if provided
            output_ds = gdal.Warp(fname, ds, xRes=spatial_resolution, yRes=-spatial_resolution, resampleAlg='near')
            # output_ds = driver.CreateCopy(fname, ds)
        else:
            # Use default spatial resolution
            output_ds = driver.CreateCopy(fname, ds)

        # Close datasets
        ds = None
        output_ds = None

        print(f'Image saved to {fname}')
        return 'success'

def mosaic_tifs(input_IMGs, output_mosaic):
    # List of input raster files
    input_files = input_IMGs
    # Output mosaic file
    output_file = output_mosaic
    # Warp options - here we set the format to GeoTIFF
    warp_options = gdal.WarpOptions(format='GTiff', resampleAlg='near', creationOptions  = "COMPRESS=LZW")

    # Mosaic the input files
    gdal.Warp(output_file, input_files, options=warp_options)
    print(f"Mosaic created: {output_mosaic}")

def reclassify_tif(input_tif, output_tif):
    # reclass_table (dict): A dictionary representing the reclassification table. originial value: new value

    reclass_table = {11:1,21:2,22:2,23:2,24:2,31:3,32:3,33:3,12:3,42:4,41:5,43:4,90:4,51:6,52:6,71:7,85:7,95:7,61:8,81:8,82:8,83:8,84:8}
    palette = {1: (85,102,170,255),2: (117,107,177,255),3: (221,204,102,255),4: (102,136,34,255),5: (85,170,34,255),6: (153,187,85,255),7: (68,170,136,255),8: (17,119,51,255)}
    color_table = gdal.ColorTable()

    for key, val in palette.items():
        color_table.SetColorEntry(key, val)

    # Convert the reclassification table to a string
    reclass_formula = "+".join([f"(A=={key})*{value}" for key, value in reclass_table.items()]) 
    # ds = gdal_calc.Calc(calc="(A-B)/(A+B)", A="nir.tif", A_band=1, B="red.tif", B_band = 1, outfile="ndvi.tif") 

    ds = gdal_calc.Calc(reclass_formula, A=input_tif, outfile=output_tif, NoDataValue=0, quiet=True,color_table=color_table)
    ds = None

def get_bbox_4326(shapefile_path):
    target_srs = osr.SpatialReference()
    target_srs.ImportFromEPSG(4326)
    target_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    srid_t = target_srs.GetAuthorityCode(None)
    # print(srid_t)  
    driver = ogr.GetDriverByName('ESRI Shapefile')
    in_memory = ogr.GetDriverByName('Memory')

    
    # Open the shapefile
    dataset = driver.Open(shapefile_path)

    if dataset is None:
        print(f'Failed to open shapefile: {shapefile_path}')
    else:
        # Get the layer
        layer = dataset.GetLayer()
        
        # Get the spatial reference and reproject to 4326 if it is different
        layer_srs = layer.GetSpatialRef()
        srid_s = layer_srs.GetAuthorityCode(None)
        # print(srid_s)
        if srid_s != srid_t:
            # Create the target shapefile
            target_ds = in_memory.CreateDataSource('')

            
            # Create the target layer
            target_layer = target_ds.CreateLayer('reprojected_layer', geom_type=layer.GetGeomType())
            
            # Create a transformation from the layer's spatial reference to EPSG:4326
            transform = osr.CoordinateTransformation(layer_srs, target_srs)

            # Apply the coordinate transformation to each feature and add to the target layer
            for feature in layer:
                geometry = feature.GetGeometryRef()
                # print(geometry)
                geometry.Transform(transform)
                # print(geometry)
                
                target_feature = ogr.Feature(target_layer.GetLayerDefn())
                target_feature.SetGeometry(geometry)
                target_layer.CreateFeature(target_feature)
                target_feature = None  # Release the target feature

            # Get the extent
            minX, maxX, minY, maxY = target_layer.GetExtent()
            # print(target_layer.GetExtent())
        else:
            # print('no trans needed')
            minX, maxX, minY, maxY = layer.GetExtent()
            # print(layer.GetExtent())

        # Close the dataset
        dataset = None
        target_ds = None
        return  minX, maxX, minY, maxY

if __name__ == "__main__":
    main()