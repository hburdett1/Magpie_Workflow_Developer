import pandas as pd
import numpy as np
import os
import sys
import subprocess
import json
import geopandas as gpd
import shutil
import rioxarray as rxr
import matplotlib.pyplot as plt
from IPython.display import display
from glob import glob
# Earth Engine
import ee
import requests

def check_projection(shp_file_path, main_dir, temporary_dir):
    # Print section header
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Check projection of shapefile')
    print('-----------------------------------------------------------------------------------------------')

    # Find the name of the shapefile in the given path
    for shp_file in os.listdir(shp_file_path):
        if shp_file.endswith(".shp"):
            shp_file_name = shp_file

    # Print shapefile information
    print('Shapefile name: ', shp_file_name)
    print('Shapefile path: ', os.path.join(main_dir, "shapefile", shp_file_name))

    # Read the shapefile into a GeoDataFrame to check its CRS
    shp_lyr_check = gpd.read_file(os.path.join(main_dir, "shapefile", shp_file_name))
    print('Shapefile CRS: ', shp_lyr_check.crs)

    # Create a temporary directory if it doesn't exist
    temp_dir = os.path.join(temporary_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # Check if the CRS is EPSG:4326, if not, reproject the shapefile
    if shp_lyr_check.crs != 'EPSG:4326':
        # Reproject the shapefile layer to match EPSG:4326
        shp_lyr_crs = shp_lyr_check.to_crs(epsg=4326)
        shp_lyr_crs.to_file(os.path.join(temporary_dir, shp_file_name))
        print('Shapefile layer has been reprojected to match EPSG:4326')
    else:
        shp_lyr_crs = shp_lyr_check
        print('Coordinate systems match!')
        shp_lyr_crs.to_file(os.path.join(temporary_dir, shp_file_name))

    # Return the name of the reprojected shapefile
    return shp_file_name

def download_soil(shp_file_path, shp_file_name, data_source, band_name, scale, temporary_dir):
    # Print section header
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Download Soil')
    print('-----------------------------------------------------------------------------------------------')

    # Define buffer size
    buffer_size = 0.3

    # Determine the boundary of the provided shapefile
    bounds = gpd.read_file(os.path.join(shp_file_path, shp_file_name)).bounds
    west, south, east, north = bounds = bounds.loc[0]
    west -= buffer_size
    south -= buffer_size
    print('Bounding box: ', west, south, east, north)

    # Concatenate input data to generate the full path
    #full_data_source = f'{data_source}/{year_of_interest}_01_01'
    full_data_source = f'{data_source}'

    # Create an Earth Engine image
    img = ee.Image(full_data_source)
    region = ee.Geometry.BBox(west, south, east, north)

    # Get the download URL for the image
    url = img.getDownloadUrl({
        'bands': [band_name],
        'region': region,
        'scale': scale,
        'format': 'GEO_TIFF'
    })

    # Define the output directory
    soil_dir = os.path.join(temporary_dir, 'Soil')
    if not os.path.exists(soil_dir):
        os.makedirs(soil_dir)

    # Download the image using the generated URL
    response = requests.get(url)
    with open(os.path.join(soil_dir, 'study_area_soil.tif'), 'wb') as fd:
        fd.write(response.content)

    # Define paths
    clipped_bounds = os.path.join(soil_dir, 'study_area_soil.tif')
    soil_filled = os.path.join(soil_dir, 'soil_filled.tif')

    # Define the directory for bash scripts
    bash_dir = os.path.join(temporary_dir, 'bash_scripts')
    if not os.path.exists(bash_dir):
        os.makedirs(bash_dir)
        print("Created folder:", bash_dir)

    # GDAL fill no data
    with open(os.path.join(bash_dir, 'fillnodata.sh'), 'w') as f1:
        print(f'gdal_fillnodata.py -md 10 -b 1 -of GTiff "{clipped_bounds}" "{soil_filled}"', file=f1)

    # Format and run the bash command
    fill_data = os.path.join(bash_dir, 'fillnodata.sh')
    subprocess.run(['bash', fill_data])

    # Remove the old unfilled soil layer
    os.remove(clipped_bounds)

    return soil_dir

def overlay_shp_on_soil(temporary_dir, shp_file_name):
    # Print section header
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Overlay shapefile on soil')
    print('-----------------------------------------------------------------------------------------------')

    # Find the name of the raster file in the temporary directory
    for tif_file in os.listdir(os.path.join(temporary_dir, 'Soil')):
        if tif_file.endswith(".tif"):
            tif_file_name = tif_file

    # Open the soil raster
    soil_lyr = rxr.open_rasterio(os.path.join(temporary_dir, 'Soil', tif_file_name), masked=True).squeeze()

    # Load the shapefile
    crop_extent = gpd.read_file(os.path.join(temporary_dir, shp_file_name))

    # Print CRS information
    print('Shapefile CRS: ', crop_extent.crs)
    print('Soil CRS: ', soil_lyr.rio.crs)

    # Check if coordinate systems match, reproject if necessary
    if crop_extent.crs != soil_lyr.rio.crs:
        soil_lyr = soil_lyr.rio.reproject(crop_extent.crs)
        print('Soil layer has been reprojected to match the shapefile')
    else:
        print('Coordinate systems match!')

    # Plot the overlay
    f, ax = plt.subplots(figsize=(10, 5))
    soil_lyr.plot.imshow(ax=ax)
    crop_extent.plot(ax=ax, alpha=0.8, color="black")
    ax.set(title="Raster Layer with Shapefile Overlayed")
    ax.set_axis_off()

    soil_shapfile_visualization = True

    if soil_shapfile_visualization:
        plt.show()

    # Return the soil layer and crop extent
    return soil_lyr, crop_extent

def clip_and_format(soil_dir, soil_lyr, crop_extent, temporary_dir):
    # Print section header
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Clip and format soil into shapefile')
    print('-----------------------------------------------------------------------------------------------')

    # Clip soil product to the study area
    # Open the crop extent (the study area extent boundary)
    crop_extent_buffered = crop_extent.buffer(.001)

    # Clip the soil layer
    soil_clipped = soil_lyr.rio.clip(crop_extent_buffered, crop_extent_buffered.crs)
    print('Soil layer has been clipped')

    # Save the clipped soil layer to the temporary folder
    path_to_tif_file = os.path.join(soil_dir, 'clipped_lyr.tif')
    soil_clipped.rio.to_raster(path_to_tif_file)
    print('Layer has been saved to the temporary folder')

    # Define necessary paths
    clip_name = 'studyArea_outline'
    clipped_lyr = os.path.join(soil_dir, 'clipped_lyr.tif')
    soil_shp = os.path.join(soil_dir, 'studyArea_soil.shp')

    # Define the bash script directory
    bash_dir = os.path.join(temporary_dir, 'bash_scripts')
    if not os.path.exists(bash_dir):
        os.makedirs(bash_dir)
        print("Created folder: ", bash_dir)

    # GDAL polygonize
    with open(os.path.join(bash_dir, 'polygon.sh'), 'w') as f3:
        print(f'gdal_polygonize.py "{clipped_lyr}" "{soil_shp}" -b 1 -f "ESRI Shapefile"', file=f3)

    # Define the path to the bash command
    polygon_sh = os.path.join(bash_dir, 'polygon.sh')

    # Run the bash command
    subprocess.run(['bash', polygon_sh])

def visualize_soil(soil_dir, main_dir):
    # Print section header
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Visualize soil outputs')
    print('-----------------------------------------------------------------------------------------------')

    # Find the name of the shapefile in the soil directory
    for soil_shp_file in os.listdir(soil_dir):
        if soil_shp_file.endswith(".shp"):
            soil_shp_file_name = soil_shp_file

    # Read the soil shapefile
    soil_shp = gpd.read_file(os.path.join(soil_dir, soil_shp_file_name))

    # Dissolve the soil shapefile by 'DN' field
    soil_dissolve = soil_shp.dissolve(by='DN')
    soil_dissolve["Soil_ID"] = soil_dissolve.index
    soil_final = soil_dissolve.reset_index()
    soil_final = soil_final.drop('DN', axis=1)

    # Display the dissolved soil shapefile
    display(soil_final)

    # Create the output directory for soil shapefile
    soil_output_dir = os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'Soil')
    if not os.path.exists(soil_output_dir):
        os.makedirs(soil_output_dir)

    # Save the final soil shapefile to the output directory
    soil_final.to_file(os.path.join(soil_output_dir, 'studyArea_soil.shp'))

    print('------------------------------------- Soil -------------------------------------')

    # Read the final soil shapefile
    soil_final = gpd.read_file(os.path.join(soil_output_dir, 'studyArea_soil.shp'))

    # Plot the soil shapefile
    f, ax = plt.subplots(figsize=(8, 9))
    soil_final.plot(column='Soil_ID', categorical=True, legend=True, ax=ax)
    ax.set_axis_off()
    plt.show()

    # Calculate and print the percentage area for each soil type
    soil_final['area'] = soil_final.area
    df_S1 = pd.DataFrame(soil_final.drop(columns='geometry'))
    df_soil = df_S1.groupby('Soil_ID').sum()
    df_soil.reset_index(inplace=True)

    # Collect unique soil IDs
    unique_soil = df_soil['Soil_ID'].unique()
    un_soi_len_lst = list(range(unique_soil.size))

    for val_s in unique_soil:
        area_poly_S = df_soil.loc[df_soil['Soil_ID'] == val_s, 'area']
        soil_area = float((area_poly_S / df_soil['area'].sum()) * 100)
        for num in un_soi_len_lst:
            if val_s == unique_soil[num]:
                print(f'ID ({val_s}): {round(soil_area, 2)}%')

    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Soil is complete!')
    print('-----------------------------------------------------------------------------------------------')

def upload_soil(soil_temp_dir, main_dir):
    # Find the name of the shapefile in the shapefile directory
    for shp_file in os.listdir(os.path.join(main_dir, 'shapefile')):
        if shp_file.endswith(".shp"):
            shp_file_name = shp_file

    # Find the name of the shapefile in the soil temporary directory
    for soil_file in os.listdir(soil_temp_dir):
        if soil_file.endswith(".shp"):
            soil_file_name = soil_file

    # Print section header
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Check Projection')
    print('-----------------------------------------------------------------------------------------------')

    # Load the study area shapefile
    shp_extent = gpd.read_file(os.path.join(main_dir, 'shapefile', shp_file_name))

    # Load the soil shapefile
    soil_extent = gpd.read_file(os.path.join(soil_temp_dir, soil_file_name))

    # Check if 'Soil_ID' column exists in soil shapefile
    if 'Soil_ID' in soil_extent:
        # Check if coordinate systems match, reproject if necessary
        if shp_extent.crs != soil_extent.crs:
            # Reproject soil layer to match the study area shapefile
            soil_lyr = soil_extent.to_crs(shp_extent.crs)
            print('\nSoil layer has been reprojected to match the shapefile')
        else:
            soil_lyr = soil_extent
            print('\nCoordinate systems match')

        # Print section header for visualization
        print('\n-----------------------------------------------------------------------------------------------')
        print('( ) Visualize and save soil layer')
        print('-----------------------------------------------------------------------------------------------')

        # Display the soil layer table
        display(soil_lyr)

        # Visualize the soil layer
        f, ax = plt.subplots(figsize=(8, 9))
        soil_lyr.plot(column='Soil_ID', categorical=True, legend=True, ax=ax)
        ax.set_axis_off()
        plt.show()

        # Calculate the area for each soil type
        soil_lyr['area'] = soil_lyr.area
        df_S1 = pd.DataFrame(soil_lyr.drop(columns='geometry'))
        df_soil = df_S1.groupby('Soil_ID').sum()
        df_soil.reset_index(inplace=True)

        # Collect unique soil IDs
        unique_soil = df_soil['Soil_ID'].unique()
        un_soi_len_lst = list(range(unique_soil.size))

        # Print the percentage area for each soil type
        for val_s in unique_soil:
            area_poly_S = df_soil.loc[df_soil['Soil_ID'] == val_s, 'area']
            soil_area = float((area_poly_S / df_soil['area'].sum()) * 100)
            for num in un_soi_len_lst:
                if val_s == unique_soil[num]:
                    print(f'ID ({val_s}): {round(soil_area, 2)}%')

        # Save the final soil shapefile to the output directory
        soil_output_dir = os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'Soil')
        if not os.path.exists(soil_output_dir):
            os.makedirs(soil_output_dir)
        soil_lyr.to_file(os.path.join(soil_output_dir, 'studyArea_soil.shp'))
    else:
        print('--- Invalid soil column name ---\n')
        print('Please change the soil column name to Soil_ID and run again')

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

if config_file_info['generate_soil'] == 'no':
    # define the output directory
    soil_dir = os.path.join(temporary_dir, 'Soil')
    if not os.path.exists(soil_dir):
        os.makedirs(soil_dir)
    # define temporary directory
    soil_temp_dir = os.path.join(temporary_dir,'Soil')
    if not os.path.exists(soil_temp_dir):
      os.makedirs(soil_temp_dir)

    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Upload Soil')
    print('-----------------------------------------------------------------------------------------------')
    print(f'drag-and-drop soil file into following folder: {soil_temp_dir}')
    response = input("Have you uploaded the soil file (yes or no): ")
    if response == "yes":
      for soil_file in os.listdir(soil_temp_dir):
          if soil_file.endswith(".tif"):
            shp_file_path = os.path.join(main_dir, 'shapefile')
            # step 1
            shp_file_name = check_projection(shp_file_path, main_dir, temporary_dir)
            # step 2
            soil_lyr, crop_extent = overlay_shp_on_soil(temporary_dir, shp_file_name)
            # step 3
            clip_and_format(soil_dir, soil_lyr, crop_extent, temporary_dir)
            # step 4
            visualize_soil(soil_dir, main_dir)
            remove_temp_data(temporary_dir)
          elif soil_file.endswith(".shp"):
            upload_soil(soil_temp_dir, main_dir)
            remove_temp_data(temporary_dir)

if config_file_info['generate_soil'] == 'yes':
    # Trigger the authentication flow
    service_account = 'magpie-developer@magpie-id-409519.iam.gserviceaccount.com'
    credentials = ee.ServiceAccountCredentials(service_account, os.path.join(main_dir,'extras','magpie-key.json'))

    # Initialize the library
    ee.Initialize(credentials)

    # define input data sources
    data_source = config_file_info['data_source_landcover']
    band_name = config_file_info['band_name_landcover']
    scale = int(config_file_info['scale_landcover'] )

    # define path
    shp_file_path = os.path.join(main_dir, 'shapefile')

    # determine shapefile name
    shp_file_name = check_projection(shp_file_path, main_dir, temporary_dir)
    # download soil data
    soil_dir = download_soil(shp_file_path, shp_file_name, data_source, band_name, scale, temporary_dir)
    soil_lyr, crop_extent = overlay_shp_on_soil(temporary_dir, shp_file_name)
    # clip and format
    clip_and_format(soil_dir, soil_lyr, crop_extent, temporary_dir)
    # visualize output
    visualize_soil(soil_dir, main_dir)
    # delete temorary folder
    remove_temp_data(temporary_dir)
"""

