print('\n-----------------------------------------------------------------------------------------------')
print('( ) Loading libraries')
print('-----------------------------------------------------------------------------------------------')
import pandas as pd
import numpy as np
import os
import json
import geopandas as gpd
import shutil
import matplotlib.pyplot as plt
from glob import glob
import rioxarray as rxr
import richdem as rd
import wget
import zipfile
import subprocess
import sys
# Earth Engine
import ee
import requests

def check_projection(shp_file_path, main_dir,temporary_dir):
  print('\n-----------------------------------------------------------------------------------------------')
  print('( ) Check projection of shapefile')
  print('-----------------------------------------------------------------------------------------------')
  # find name of shapefile
  for shp_file in os.listdir(os.path.join(shp_file_path)):
      if shp_file.endswith(".shp"):
        shp_file_name = shp_file
  print('Shapefile name: ', shp_file_name)
  print('Shapefile path: ', os.path.join(main_dir, "shapefile",shp_file_name))
  shp_lyr_check = gpd.read_file(os.path.join(main_dir, "shapefile",shp_file_name))
  print('Shapefile CRS: ', shp_lyr_check.crs)

  if shp_lyr_check.crs !=  'EPSG:4326':
    # reproject
    shp_lyr_crs = shp_lyr_check.to_crs(epsg=4326)
    shp_lyr_crs.to_file(os.path.join(main_dir, "shapefile",shp_file))
    print('Shapefile layer has been reprojected to match shapefile')
  else:
    shp_lyr_crs = shp_lyr_check
    print('Coordinate systems match!')

def download_DEM(shp_file_path, data_source, band_name, scale, main_dir, temporary_dir):
    """
    Download a Digital Elevation Model (DEM) based on the provided shapefile and parameters.

    Args:
        shp_file_path (str): Path to the directory containing the shapefile.
        data_source: Earth Engine data source.
        band_name (str): Name of the band in the DEM.
        scale (float): Scale for the DEM download.
        main_dir (str): Main directory path.
        temporary_dir (str): Temporary directory path.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Download DEM')
    print('-----------------------------------------------------------------------------------------------')

    # Define buffer size
    buffer_size = 0.01

    # Find name of shapefile
    shp_file_path = os.path.join(main_dir, 'shapefile')
    for shp_file in os.listdir(os.path.join(shp_file_path)):
        if shp_file.endswith(".shp"):
            shp_file_name = shp_file

    # Determine the boundary of the provided shapefile
    bounds = gpd.read_file(os.path.join(main_dir, 'shapefile', shp_file_name)).bounds
    west, south, east, north = bounds = bounds.loc[0]
    west -= buffer_size * (east - west)
    east += buffer_size * (east - west)
    south -= buffer_size * (north - south)
    north += buffer_size * (north - south)

    img = ee.Image(data_source)
    region = ee.Geometry.BBox(west, south, east, north)

    # Multi-band GeoTIFF file.
    url = img.getDownloadUrl({
        'bands': [band_name],
        'region': region,
        'scale': scale,
        'format': 'GEO_TIFF'
    })

    # Define output directory for DEM
    dem_dir = os.path.join(temporary_dir, 'DEM')
    if not os.path.exists(dem_dir):
        os.makedirs(dem_dir)

    # Download DEM using requests
    response = requests.get(url)
    with open(os.path.join(dem_dir, 'dem.tif'), 'wb') as fd:
        fd.write(response.content)

def format_and_visualize_dem(shp_file_path, temporary_dir, main_dir):
    """
    Clip DEM to the study area, reproject if necessary, and save the clipped DEM to drive.

    Args:
        shp_file_path (str): Path to the directory containing the shapefile.
        temporary_dir (str): Temporary directory path.
        main_dir (str): Main directory path.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Clip DEM to study area')
    print('-----------------------------------------------------------------------------------------------')

    # Find name of shapefile
    for shp_file in os.listdir(os.path.join(shp_file_path)):
        if shp_file.endswith(".shp"):
            shp_file_name = shp_file

    # Find name of DEM file
    for dem_file in os.listdir(os.path.join(os.path.join(temporary_dir, 'DEM'))):
        if dem_file.endswith(".tif"):
            dem_file_name = dem_file

    # Open raster DEM layer
    dem_lyr = rxr.open_rasterio(os.path.join(temporary_dir, 'DEM', dem_file_name), masked=True).squeeze()

    # Load shapefile (crop extent)
    crop_extent = gpd.read_file(os.path.join(main_dir, 'shapefile', shp_file_name))

    print('Shapefile CRS:', crop_extent.crs)
    print('DEM CRS:', dem_lyr.rio.crs)

    # Check if CRS match and reproject if necessary
    if crop_extent.crs != dem_lyr.rio.crs:
        dem_lyr = dem_lyr.rio.reproject(crop_extent.crs)
        print('\nDEM layer has been reprojected to match shapefile')
    else:
        print('\nCoordinate systems match')

    # Open crop extent (study area extent boundary)
    crop_extent_buffered = crop_extent.buffer(0.001)

    # Clip the DEM layer
    lidar_clipped = dem_lyr.rio.clip(crop_extent_buffered, crop_extent_buffered.crs)
    print('\nDEM layer has been clipped')

    # Define output directory path
    dem_output_dir = os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'DEM')

    if not os.path.exists(dem_output_dir):
        os.makedirs(dem_output_dir)

    # Save clipped DEM layer to drive
    path_to_tif_file = os.path.join(dem_output_dir, 'dem.tif')
    lidar_clipped.rio.to_raster(path_to_tif_file)
    print('\nDEM layer has been saved to drive')


def visualize_dem(DEM_visualize, Slope_visualize, Aspect_visualize, save_aspect, main_dir, temporary_dir):
    """
    Visualize DEM layer, produce slope and aspect layers, and save them to drive.

    Args:
        DEM_visualize (str): 'yes' to visualize DEM layer, 'no' otherwise.
        Slope_visualize (str): 'yes' to visualize slope layer, 'no' otherwise.
        Aspect_visualize (str): 'yes' to visualize aspect layer, 'no' otherwise.
        save_shp_of_aspect (str): 'yes' to save shapefile of aspect layer, 'no' otherwise.
        main_dir (str): Main directory path.
        temporary_dir (str): Temporary directory path.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Visualize DEM layer and produce slope and aspect layers')
    print('-----------------------------------------------------------------------------------------------')

    # Define output directory path for slope and aspect
    slope_output_dir = os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'Slope')
    aspect_output_dir = os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'Aspect')

    # Create output directories if they don't exist
    for output_dir in [slope_output_dir, aspect_output_dir]:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    if DEM_visualize == 'yes':
        for dem_file in os.listdir(os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'DEM')):
            if dem_file.endswith(".tif"):
                print('---------------------------------------------- DEM Layer ----------------------------------------------')
                # Open and visualize DEM layer
                dem_lyr = rxr.open_rasterio(os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'DEM', dem_file),
                                            masked=True).squeeze()
                f, ax = plt.subplots(figsize=(8, 10))
                dem_lyr.plot(ax=ax)
                ax.set_axis_off()
                plt.show()

                # Calculate and print statistics of DEM layer
                dem_min = int(dem_lyr.min())
                dem_max = int(dem_lyr.max())
                dem_mean = int(dem_lyr.mean())
                print('\nMinimum Elevation:', dem_min)
                print('Maximum Elevation:', dem_max)
                print('Mean Elevation:', dem_mean)

    if Slope_visualize == 'yes':
        for dem_file in os.listdir(os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'DEM')):
            if dem_file.endswith(".tif"):
                print('---------------------------------------------- Slope Layer ----------------------------------------------')
                # Using richdem to calculate slope and visualize
                dem_lyr = rd.LoadGDAL(os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'DEM', dem_file),
                                     no_data=-9999)
                slope_lyr = rd.TerrainAttribute(dem_lyr, attrib="slope_degrees")
                rd.SaveGDAL(os.path.join(slope_output_dir, "slope.tif"), slope_lyr)
                rd.rdShow(slope_lyr, axes=False, cmap='magma', figsize=(8, 9))
                plt.show()

                # Open and clean slope data
                slope_data = rxr.open_rasterio(os.path.join(slope_output_dir, "slope.tif"), masked=True).squeeze()
                slope_clean = slope_data.where(slope_data != -9999, np.nan)
                slope_min = int(slope_clean.min())
                slope_max = int(slope_clean.max())
                slope_mean = int(slope_clean.mean())
                print('\nMinimum Slope:', slope_min)
                print('Maximum Slope:', slope_max)
                print('Mean Slope:', slope_mean)

    if Aspect_visualize == 'yes':
        for dem_file in os.listdir(os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'DEM')):
            if dem_file.endswith(".tif"):
                print('---------------------------------------------- Aspect Layer ----------------------------------------------')
                # Using richdem to calculate aspect and visualize
                dem_lyr = rd.LoadGDAL(os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'DEM', dem_file),
                                     no_data=-9999)
                aspect_lyr = rd.TerrainAttribute(dem_lyr, attrib="aspect")
                rd.SaveGDAL(os.path.join(aspect_output_dir, "aspect.tif"), aspect_lyr)
                rd.rdShow(aspect_lyr, axes=False, cmap='jet', figsize=(8, 9))
                plt.show()

                # Open aspect data and calculate statistics
                aspect_data = rxr.open_rasterio(os.path.join(aspect_output_dir, 'aspect.tif'), masked=True).squeeze()
                aspect_min = int(aspect_data.min())
                aspect_max = int(aspect_data.max())
                aspect_mean = int(aspect_data.mean())
                print('\nMinimum Aspect:', aspect_min)
                print('Maximum Aspect:', aspect_max)
                print('Mean Aspect:', aspect_mean)

    if save_aspect == "yes":
        # Save shapefile of aspect layer
        aspect_rast = os.path.join(aspect_output_dir, 'aspect.tif')
        aspect_shp = os.path.join(aspect_output_dir, 'aspect.shp')

        # Use gdal_polygonize to convert raster to vector (shapefile)
        with open(os.path.join(temporary_dir, 'polygon.sh'), 'w') as f3:
            print(f'gdal_polygonize.py "{aspect_rast}" "{aspect_shp}" -b 1 -f "ESRI Shapefile"', file=f3)

        sh_file = os.path.join(temporary_dir, 'polygon.sh')
        subprocess.run(['bash', sh_file])

        # Open and format the shapefile
        aspect_shp_gdf = gpd.read_file(aspect_shp)
        aspect_shp_gdf["O_ID_2"] = list(range(1, (len(aspect_shp_gdf.index) + 1)))
        aspect_shp_gdf["Aspect"] = aspect_shp_gdf.DN

        # Save the final elevation band shapefile to drive
        aspect_shp_gdf.to_file(aspect_shp)
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) DEM Complete!')
    print('-----------------------------------------------------------------------------------------------')

def remove_temp_data(temporary_dir):
  if os.path.exists(os.path.join(temporary_dir)):
    shutil.rmtree(os.path.join(temporary_dir))

"""
# define main dir
main_dir = sys.argv[1]

# read in information from configuration file
with open(os.path.join(main_dir,"configuration_file.json"), "r") as f:
    config_file_info = json.load(f)

# define temporary dir
temporary_dir = config_file_info['temporary_dir']


if config_file_info['generate_DEM'] == "no":
    dem_temp_dir = os.path.join(temporary_dir,'DEM')
    if not os.path.exists(dem_temp_dir):
      os.makedirs(dem_temp_dir)
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Upload DEM')
    print('-----------------------------------------------------------------------------------------------')
    print(f'drag-and-drop DEM file into following folder: {dem_temp_dir}')
    response = input("Have you uploaded the DEM (.tif) file (yes or no): ")
    if response == "yes":
      shp_file_path = os.path.join(main_dir, 'shapefile')
      # check projection
      check_projection(shp_file_path, main_dir,temporary_dir)
      # format and visualize
      format_and_visualize_dem(shp_file_path,temporary_dir,main_dir)

elif config_file_info['generate_DEM'] == 'yes':
    # Trigger the authentication flow.
    service_account = 'magpie-developer@magpie-id-409519.iam.gserviceaccount.com'
    credentials = ee.ServiceAccountCredentials(service_account, os.path.join(main_dir,'extras','magpie-key.json'))

    # Initialize the library.
    ee.Initialize(credentials)
    
    # define input variables
    shp_file_path = os.path.join(main_dir, 'shapefile')
    data_source = config_file_info['data_source_dem']
    band_name = config_file_info['band_name_dem']
    scale = int(config_file_info['scale_dem'])
    DEM_visualize = config_file_info['DEM_visualize']
    Slope_visualize = config_file_info['Slope_visualize']
    Aspect_visualize = config_file_info['Aspect_visualize']
    save_aspect = config_file_info['save_aspect']

    # check projection
    check_projection(shp_file_path, main_dir,temporary_dir)
    # download DEM
    download_DEM(shp_file_path, data_source, band_name, scale, main_dir, temporary_dir)
    # clip DEM
    format_and_visualize_dem(shp_file_path,temporary_dir,main_dir)
    # visualize final output
    visualize_dem(DEM_visualize, Slope_visualize, Aspect_visualize, save_aspect, main_dir, temporary_dir)
    # delete temorary folder
    remove_temp_data(temporary_dir)
"""