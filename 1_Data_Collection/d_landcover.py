import pandas as pd
import numpy as np
import os
import sys
import subprocess
import json
from IPython.display import display
import geopandas as gpd
import shutil
import rioxarray as rxr
import matplotlib.pyplot as plt
from glob import glob
# Earth Engine
import ee
import requests

def check_projection(shp_file_path, main_dir, temporary_dir):
    """
    Check the projection of a shapefile and reproject if necessary.

    Parameters:
    - shp_file_path: Path to the shapefile.
    - main_dir: Main directory containing the shapefile.
    - temporary_dir: Temporary directory for storing reprojected shapefiles.

    Returns:
    - shp_file_name: Name of the shapefile.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Check projection of shapefile')
    print('-----------------------------------------------------------------------------------------------')

    # Find name of shapefile
    for shp_file in os.listdir(os.path.join(shp_file_path)):
        if shp_file.endswith(".shp"):
            shp_file_name = shp_file

    print('Shapefile name: ', shp_file_name)
    print('Shapefile path: ', os.path.join(main_dir, "shapefile", shp_file_name))

    shp_lyr_check = gpd.read_file(os.path.join(main_dir, "shapefile", shp_file_name))
    print('Shapefile CRS: ', shp_lyr_check.crs)

    temp_dir = os.path.join(temporary_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    if shp_lyr_check.crs != 'EPSG:4326':
        # Reproject
        shp_lyr_crs = shp_lyr_check.to_crs(epsg=4326)
        shp_lyr_crs.to_file(os.path.join(temporary_dir, shp_file_name))
        print('Shapefile layer has been reprojected to match shapefile')
    else:
        shp_lyr_crs = shp_lyr_check
        print('Coordinate systems match!')
        shp_lyr_crs.to_file(os.path.join(temporary_dir, shp_file_name))

    return shp_file_name

def download_landcover(shp_file_path, shp_file_name, data_source, year_of_interest, band_name, scale, temporary_dir):
    """
    Download landcover data based on the given shapefile bounding box.

    Parameters:
    - shp_file_path: Path to the shapefile.
    - shp_file_name: Name of the shapefile.
    - data_source: Source of landcover data.
    - year_of_interest: Year of the landcover data.
    - band_name: Name of the band.
    - scale: Scale of the download.
    - temporary_dir: Temporary directory for storing downloaded landcover data.

    Returns:
    - landcov_dir: Directory path where the landcover data is stored.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Download Landcover')
    print('-----------------------------------------------------------------------------------------------')


    # Define buffer
    buffer_size = 0.3

    # Determine the boundary of the provided shapefile
    bounds = gpd.read_file(os.path.join(shp_file_path, shp_file_name)).bounds
    west, south, east, north = bounds = bounds.loc[0]
    west -= buffer_size
    south -= buffer_size

    print('Bounding box: ', west, south, east, north)

    # Concatenate input data to generate full path
    full_data_source = f"{data_source}/{year_of_interest}_01_01"
    img = ee.Image(full_data_source)
    region = ee.Geometry.BBox(west, south, east, north)

    # Multi-band GeoTIFF file.
    url = img.getDownloadUrl({
        'bands': band_name,
        'region': region,
        'scale': scale,
        'format': 'GEO_TIFF'
    })

    # Define output directory
    landcov_dir = os.path.join(temporary_dir, 'Landcover')
    if not os.path.exists(landcov_dir):
        os.makedirs(landcov_dir)

    response = requests.get(url)
    with open(os.path.join(landcov_dir, 'study_area_landcov.tif'), 'wb') as fd:
        fd.write(response.content)

    # Path to clipped output file
    clipped_bounds = os.path.join(landcov_dir, 'study_area_landcov.tif')
    return landcov_dir

def upload_landcover(land_temp_dir, main_dir, temporary_dir):
    # Find the name of the study area shapefile
    for shp_file in os.listdir(os.path.join(main_dir, 'shapefile')):
        if shp_file.endswith(".shp"):
            shp_file_name = shp_file

    # Find the name of the landcover shapefile
    for landcov_file in os.listdir(land_temp_dir):
        if landcov_file.endswith(".shp"):
            landcov_file_name = landcov_file

    # Print section header
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Check Projection')
    print('-----------------------------------------------------------------------------------------------')

    # Load study area shapefile
    shp_extent = gpd.read_file(os.path.join(main_dir, 'shapefile', shp_file_name))
    # Load landcover shapefile
    landcov_extent = gpd.read_file(os.path.join(land_temp_dir, landcov_file_name))

    # Check if the landcover shapefile has the required column 'Landuse_ID'
    if 'Landuse_ID' in landcov_extent:
        # Check if coordinate systems match, reproject if necessary
        if shp_extent.crs != landcov_extent.crs:
            # Reproject landcover layer to match the study area shapefile
            landcov_lyr = landcov_extent.to_crs(shp_extent.crs)
            print('\nLandcover layer has been reprojected to match the shapefile')
        else:
            landcov_lyr = landcov_extent
            print('\nCoordinate systems match')

        # Print section header
        print('\n-----------------------------------------------------------------------------------------------')
        print('( ) Visualize and save landcover layer')
        print('-----------------------------------------------------------------------------------------------')

        # Display landcover GeoDataFrame table
        display(landcov_lyr)

        # Visualize landcover
        f, ax = plt.subplots(figsize=(8, 9))
        landcov_lyr.plot(column='Landuse_ID', categorical=True, legend=True, ax=ax)
        ax.set_axis_off()
        plt.show()

        # Calculate and display area for each landcover type
        landcov_lyr['area'] = landcov_lyr.area
        df_L1 = pd.DataFrame(landcov_lyr.drop(columns='geometry'))
        df_landcov = df_L1.groupby('Landuse_ID').sum()
        df_landcov.reset_index(inplace=True)

        # Collect unique landcover IDs
        unique_landcov = df_landcov['Landuse_ID'].unique()

        # Compute general area percentage for each landcover type
        for val_L in unique_landcov:
            area_poly_L = df_landcov.loc[df_landcov['Landuse_ID'] == val_L, 'area']
            landcov_area = float((area_poly_L / df_landcov['area'].sum()) * 100)
            print(f'ID ({val_L}): {round(landcov_area, 2)}%')

        # Save the final landcover shapefile to drive
        landcov_output_dir = os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'Landcover')
        if not os.path.exists(landcov_output_dir):
            os.makedirs(landcov_output_dir)
        landcov_lyr.to_file(os.path.join(landcov_output_dir, 'studyArea_landcover.shp'))
    else:
        # Print an error message if the required column name is not found
        print('--- Invalid landuse column name ---\n')
        print('Please change the landuse column name to Landuse_ID and run again')


def overlay_shp_on_landcov(shp_file_path, shp_file_name, temporary_dir):
    """
    Overlay shapefile on landcover data and visualize the result.

    Parameters:
    - shp_file_path: Path to the shapefile.
    - shp_file_name: Name of the shapefile.
    - temporary_dir: Temporary directory for storing landcover data.

    Returns:
    - landcov_lyr: Landcover layer.
    - crop_extent: Crop extent layer.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Overlay shapefile on landcover')
    print('-----------------------------------------------------------------------------------------------')

    # Find name of raster
    for tif_file in os.listdir(os.path.join(temporary_dir, 'Landcover')):
        if tif_file.endswith(".tif"):
            tif_file_name = tif_file

    # Open raster
    landcov_lyr = rxr.open_rasterio(os.path.join(temporary_dir, 'Landcover', tif_file_name), masked=True).squeeze()
    # Load shapefile
    crop_extent = gpd.read_file(os.path.join(temporary_dir, shp_file_name))

    print('Shapefile CRS: ', crop_extent.crs)
    print('Landcover CRS: ', landcov_lyr.rio.crs)

    if crop_extent.crs != landcov_lyr.rio.crs:
        # Reproject
        landcov_lyr = landcov_lyr.rio.reproject(crop_extent.crs)
        print('Landcover layer has been reprojected to match the shapefile')
    else:
        print('Coordinate systems match!')

    f, ax = plt.subplots(figsize=(10, 5))
    landcov_lyr.plot.imshow(ax=ax)

    crop_extent.plot(ax=ax, alpha=.8, color="black")
    ax.set(title="Raster Layer with Shapefile Overlayed")
    ax.set_axis_off()

    landcov_shapfile_visualization = True

    if landcov_shapfile_visualization:
        plt.show()

    return landcov_lyr, crop_extent

def clip_and_format(landcov_dir, landcov_lyr, crop_extent, temporary_dir):
    """
    Clip and format landcover into a shapefile.

    Parameters:
    - landcov_dir: Directory to save the processed data.
    - landcov_lyr: Landcover layer to be clipped.
    - crop_extent: Study area extent boundary.
    - temporary_dir: Temporary directory to store intermediate files.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Clip and format landcover into shapefile')
    print('-----------------------------------------------------------------------------------------------')

    # Open crop extent (the study area extent boundary)
    crop_extent1 = crop_extent.buffer(.001)

    # Clip the landcover layer
    lidar_clipped = landcov_lyr.rio.clip(crop_extent1, crop_extent1.crs)
    print('Landcover layer has been clipped')

    # Save clipped landcover layer to drive
    path_to_tif_file = os.path.join(landcov_dir, 'clipped_lyr.tif')
    lidar_clipped.rio.to_raster(path_to_tif_file)
    print('Layer has been saved to temporary folder')

    # Define pathways necessary for GDAL Commands
    clip_name = 'studyArea_outline'
    clipped = os.path.join(landcov_dir, 'clipped.tif')
    clipped_lyr = os.path.join(landcov_dir, 'clipped_lyr.tif')
    landcov_shp = os.path.join(landcov_dir, 'studyArea_landcov.shp')
    bash_dir = os.path.join(temporary_dir, 'bash_scripts')

    # Create the bash directory if it doesn't exist
    if not os.path.exists(bash_dir):
        os.makedirs(bash_dir)
        print("Created folder: ", bash_dir)

    # GDAL polygonize
    with open(os.path.join(bash_dir, 'polygon.sh'), 'w') as f3:
        print(f'gdal_polygonize.py "{clipped_lyr}" "{landcov_shp}" -b 1 -f "ESRI Shapefile"', file=f3)

    # Define bash command path
    polygon_sh = os.path.join(bash_dir, 'polygon.sh')

    # Run bash command
    subprocess.run(['bash', polygon_sh])

def format_google_earth_data(main_dir, temporary_dir):
    """
    Format landcover attribute names.

    Parameters:
    - main_dir: Main directory where the workflow is located.
    - temporary_dir: Temporary directory to store intermediate files.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Format landcover attribute names')
    print('-----------------------------------------------------------------------------------------------')

    landcov_dir = os.path.join(temporary_dir, 'Landcover')

    # Find the name of the shapefile
    for land_shp_file in os.listdir(landcov_dir):
        if land_shp_file.endswith(".shp"):
            land_shp_file_name = land_shp_file

    landcov_shp = gpd.read_file(os.path.join(landcov_dir, land_shp_file_name))
    landcov_dissolve = landcov_shp.dissolve(by='DN')
    landcov_dissolve["Landuse_ID"] = landcov_dissolve.index
    landcov_final = landcov_dissolve.reset_index()
    landcov_final = landcov_final.drop('DN', axis=1)

    land_output_dir = os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'Landcover')
    if not os.path.exists(land_output_dir):
        os.makedirs(land_output_dir)

    # Save the final landcover shapefile to drive
    landcov_final.to_file(os.path.join(land_output_dir, 'studyArea_landcover.shp'))

def visualize_landcover(main_dir):
    # Print section header
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Visualize landcover outputs')
    print('-----------------------------------------------------------------------------------------------')

    # Create output directory if it doesn't exist
    landcov_dir_out = os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'Landcover')
    if not os.path.exists(landcov_dir_out):
        os.makedirs(landcov_dir_out)

    # Find name of shapefile in the output directory
    for land_shp_file in os.listdir(landcov_dir_out):
        if land_shp_file.endswith(".shp"):
            land_shp_file_name1 = land_shp_file

    # Read the shapefile into a GeoDataFrame
    landcov_final = gpd.read_file(os.path.join(landcov_dir_out, land_shp_file_name1))

    # Display GeoDataFrame table
    display(landcov_final)

    # Plot landcover using GeoDataFrame
    f, ax = plt.subplots(figsize=(8, 9))
    landcov_final.plot(column='Landuse_ID', categorical=True, legend=True, ax=ax)
    ax.set_axis_off()
    plt.show()

    # Calculate and display area for each landcover type
    landcov_final['area'] = landcov_final.area
    df_L1 = pd.DataFrame(landcov_final.drop(columns='geometry'))
    df_landcov = df_L1.groupby('Landuse_ID').sum()
    df_landcov.reset_index(inplace=True)

    # Collect unique landcover IDs
    unique_landcov = df_landcov['Landuse_ID'].unique()

    # Compute general area percentage for each landcover type
    for val_L in unique_landcov:
        area_poly_L = df_landcov.loc[df_landcov['Landuse_ID'] == val_L, 'area']
        landcov_area = float((area_poly_L / df_landcov['area'].sum()) * 100)
        print(f'ID ({val_L}): {round(landcov_area, 2)}%')

    # Print section footer
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Landcover is complete!')
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

if config_file_info['generate_landcover'] == 'no':
    # generate temporary folder
    land_temp_dir = os.path.join(temporary_dir,'Landcover')
    if not os.path.exists(land_temp_dir):
      os.makedirs(land_temp_dir)

    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Upload Landcover')
    print('-----------------------------------------------------------------------------------------------')
    print(f'drag-and-drop landcover file into following folder: {land_temp_dir}')
    response = input("Have you uploaded the landcover file (yes or no): ")
    if response == "yes":
      for land_file in os.listdir(land_temp_dir):
        if land_file.endswith(".tif"):
          shp_file_path = os.path.join(main_dir, 'shapefile')
          # step 1
          shp_file_name = check_projection(shp_file_path,main_dir,temporary_dir)
          # step 2
          landcov_lyr, crop_extent = overlay_shp_on_landcov(shp_file_path, shp_file_name,temporary_dir)
          # step 3
          # define output directory
          landcov_dir = os.path.join(temporary_dir,'Landcover')
          if not os.path.exists(landcov_dir):
            os.makedirs(landcov_dir)
          clip_and_format(landcov_dir, landcov_lyr, crop_extent,temporary_dir)
          format_google_earth_data(main_dir,temporary_dir)
          # step 4
          visualize_landcover(main_dir)
          remove_temp_data(temporary_dir)
        elif land_file.endswith(".shp"):
          upload_landcover(land_temp_dir,main_dir,temporary_dir)
          remove_temp_data(temporary_dir)

if config_file_info['generate_landcover'] == 'yes':
    # Trigger the authentication flow
    service_account = 'magpie-developer@magpie-id-409519.iam.gserviceaccount.com'
    credentials = ee.ServiceAccountCredentials(service_account, os.path.join(main_dir,'extras','magpie-key.json'))

    # Initialize the library
    ee.Initialize(credentials)

    # define paths
    shp_file_path = os.path.join(main_dir, 'shapefile')
    shp_file_name = check_projection(shp_file_path,main_dir,temporary_dir)

    # define input variables
    data_source = config_file_info['data_source_landcover']
    year_of_interest = config_file_info['year_of_interest']
    band_name = config_file_info['band_name_landcover']
    scale = int(config_file_info['scale_landcover'])

    # download landcover
    landcov_dir = download_landcover(shp_file_path, shp_file_name, data_source, year_of_interest, band_name, scale, temporary_dir)
    landcov_lyr, crop_extent = overlay_shp_on_landcov(shp_file_path, shp_file_name,temporary_dir)
    # clip and format
    clip_and_format(landcov_dir, landcov_lyr, crop_extent,temporary_dir)
    format_google_earth_data(main_dir,temporary_dir)
    # visualize final product
    visualize_landcover(main_dir)
    # delete temorary folder
    remove_temp_data(temporary_dir)
"""