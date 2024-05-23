print('\n-----------------------------------------------------------------------------------------------')
print('( ) Loading libraries')
print('-----------------------------------------------------------------------------------------------')
import os
import sys
import subprocess
import pandas as pd
import geopandas as gpd
from IPython.display import display
import wget
import json
import zipfile
import rioxarray as rxr
import shutil
import matplotlib.pyplot as plt

def upload_elev_bands(reclassify_dir, main_dir):
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Upload Elevation Band Shapefile')
    print('-----------------------------------------------------------------------------------------------')

    # Find the name of the shapefile in the 'shapefile' directory
    shp_file_name = next((shp_file for shp_file in os.listdir(os.path.join(main_dir, 'shapefile')) if shp_file.endswith(".shp")), None)

    # Find the name of the elevation band shapefile in the reclassify directory
    elevB_file_name = next((elevB_file for elevB_file in os.listdir(reclassify_dir) if elevB_file.endswith(".shp")), None)

    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Check Projection')
    print('-----------------------------------------------------------------------------------------------')

    # Load the study area shapefile
    shp_extent = gpd.read_file(os.path.join(main_dir, 'shapefile', shp_file_name))

    # Load the elevation band shapefile
    elevB_extent = gpd.read_file(os.path.join(reclassify_dir, elevB_file_name))

    if shp_extent.crs != elevB_extent.crs:
        # Reproject the elevation band layer to match the study area shapefile
        elevB_extent = elevB_extent.to_crs(shp_extent.crs)
        print('\nElevation Band layer has been reprojected to match the shapefile')
    else:
        print('\nCoordinate systems match')

    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Visualize and save elevation bands')
    print('-----------------------------------------------------------------------------------------------')

    # Display the table
    display(elevB_extent)

    # Visualize elevation bands
    f, ax = plt.subplots(figsize=(9, 10))
    elevB_extent.plot(categorical=False, legend=True, ax=ax)
    ax.set(title="Elevation Bands")
    ax.set_axis_off()
    plt.show()

    # Save the final elevation band shapefile to drive
    elev_band_output_dir = os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'Elevation_band')
    os.makedirs(elev_band_output_dir, exist_ok=True)
    elevB_extent.to_file(os.path.join(elev_band_output_dir, 'studyArea_elev_band.shp'))

def min_max_elev(reclassify_dir, main_dir):
    """
    Determine the minimum and maximum elevation of a Digital Elevation Model (DEM).

    Args:
        reclassify_dir (str): Directory for temporary files, including the downloaded package.
        main_dir (str): Main directory path.

    Returns:
        int: Maximum elevation of the DEM.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Determine minimum and maximum elevation')
    print('-----------------------------------------------------------------------------------------------')

    # Create reclassify_dir if it doesn't exist
    if not os.path.exists(reclassify_dir):
        os.makedirs(reclassify_dir)

    # Download and extract the gdal_reclassify package
    reclassify_url = 'https://github.com/hburdett1/Magpie_Workflow_Developer/archive/refs/heads/main.zip'
    wget.download(reclassify_url, out=reclassify_dir)

    # Unzip the downloaded package
    extension = ".zip"
    os.chdir(reclassify_dir)
    for item in os.listdir(reclassify_dir):
        if item.endswith(extension):
            file_name = os.path.abspath(item)
            zip_ref = zipfile.ZipFile(file_name)
            zip_ref.extractall(reclassify_dir)
            zip_ref.close()
            os.remove(file_name)

    # Find name of DEM file
    dem_path = os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'DEM')
    for dem_file in os.listdir(os.path.join(dem_path)):
        if dem_file.endswith(".tif"):
            dem_file_name = dem_file

    # Open DEM layer using rasterio and rioxarray
    src = rxr.open_rasterio(os.path.join(dem_path, dem_file_name), masked=True).squeeze()

    # Determine and print minimum elevation of DEM layer
    dem_min = int(src.min())
    print('\nMinimum Elevation: ', dem_min)

    # Determine and print maximum elevation of DEM layer
    dem_max = int(src.max())
    print('Maximum Elevation: ', dem_max)

    return dem_max

def reclassify_dem(increment_val, min_elev, dem_max, reclassify_dir, main_dir, temporary_dir):
    """
    Reclassify a Digital Elevation Model (DEM) based on elevation bands and convert it to a shapefile.

    Args:
        increment_val (float): Increment value for elevation bands.
        min_elev (float): Minimum elevation value.
        dem_max (float): Maximum elevation value.
        reclassify_dir (str): Directory for temporary files, including the reclassify script.
        main_dir (str): Main directory path.
        temporary_dir (str): Temporary directory path.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Reclassify DEM')
    print('-----------------------------------------------------------------------------------------------')

    # Initialize variables for elevation bands
    min_val = min_elev
    inc_lst = ['<=0', f'<{min_elev}']

    # Generate elevation bands
    while min_elev < (dem_max - increment_val):
        min_elev = min_elev + increment_val
        inc_lst.append(f'<{min_elev}')

    # Join elevation bands for gdal_reclassify command
    band_range = ','.join(str(x) for x in inc_lst)
    print('Elevation Band Increments: ', band_range)

    # Create band IDs
    id_vals = list(range(len(inc_lst)))
    band_id = ','.join(str(i) for i in id_vals)
    print('Number of Elevation Bands: ', band_id)
    print('\n')

    # Define paths for gdal commands
    gdal_reclassify = os.path.join(reclassify_dir, 'Magpie_Workflow_Developer-main', '1_Data_Collection','gdal_reclassify.py')
    elev_band_rast_path = os.path.join(reclassify_dir, 'elev_bands_rast.tif')
    elev_band_shp_path = os.path.join(reclassify_dir, 'elev_bands.shp')

    # Find name of DEM file
    dem_path = os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'DEM')
    for dem_file in os.listdir(dem_path):
        if dem_file.endswith(".tif"):
            dem_file_name = dem_file
    dem_path_full = os.path.join(dem_path, dem_file_name)

    # Write bash command to reclassify
    with open(os.path.join(temporary_dir, 'reclass.sh'), 'w') as f1:
        print(f'python "{gdal_reclassify}" "{dem_path_full}" "{elev_band_rast_path}" -c "{band_range}" -r "{band_id}" -d 0 -n true -p "COMPRESS=LZW"', file=f1)

    reclass_sh_file = os.path.join(temporary_dir, 'reclass.sh')
    subprocess.run(['bash', reclass_sh_file])

    # Write bash command to convert raster to shapefile
    with open(os.path.join(temporary_dir, 'polygonize.sh'), 'w') as f2:
        print(f'gdal_polygonize.py "{elev_band_rast_path}" "{elev_band_shp_path}" -b 1 -f "ESRI Shapefile"', file=f2)

    poly_sh_file = os.path.join(temporary_dir, 'polygonize.sh')
    subprocess.run(['bash', poly_sh_file])

def visualize_elev_bands_func(increment_val, min_elev, temporary_dir, main_dir, dem_max):
    """
    Visualize and save elevation bands based on a shapefile with elevation band information.

    Args:
        increment_val (float): Increment value for elevation bands.
        min_elev (float): Minimum elevation value.
        temporary_dir (str): Temporary directory path.
        main_dir (str): Main directory path.
        dem_max (float): Maximum elevation value.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Visualize and save elevation bands')
    print('-----------------------------------------------------------------------------------------------')

    # Load polygon shapefile
    elev_band_shp_path = os.path.join(temporary_dir, 'elev_band', 'elev_bands.shp')
    elev_band_shp = gpd.read_file(elev_band_shp_path)

    # Aggregate elevation bands based on ID values
    elev_band_dissolve = elev_band_shp.dissolve(by='DN')
    elev_band_dissolve.to_file(os.path.join(temporary_dir, 'elev_bands_dissolved.shp'))

    # Create lists of minimum and maximum contour values
    inc_min = [min_elev]
    while min_elev < (dem_max - increment_val):
        min_elev = min_elev + increment_val
        inc_min.append(min_elev)
    inc_min_lst = inc_min[:-1]
    max_inc = [max_val + increment_val for max_val in inc_min_lst]

    print('\nMinimum Contour Lines: ', inc_min_lst)
    print('Maximum Contour Lines: ', max_inc)

    # Format shapefile with contour values
    elev_band_final = gpd.read_file(os.path.join(temporary_dir, 'elev_bands_dissolved.shp'))
    elev_band_final["O_ID_1"] = list(range(1, (len(elev_band_final["DN"]) + 1)))  # Define new ID column
    elev_band_final["ContourMin"] = inc_min_lst  # Assign minimum contour values
    elev_band_final["ContourMax"] = [max_val + increment_val for max_val in inc_min_lst]  # Create list of maximum contour values
    del elev_band_final["DN"]  # Remove old ID column

    # Display the table
    display(elev_band_final)

    # Visualize elevation bands
    f, ax = plt.subplots(figsize=(9, 10))
    elev_band_shp.plot(column='DN', categorical=False, legend=True, ax=ax)
    ax.set(title="Elevation Bands")
    ax.set_axis_off()
    plt.show()

    # Save the final elevation band shapefile to drive
    elev_band_output_dir = os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'Elevation_band')
    if not os.path.exists(elev_band_output_dir):
        os.makedirs(elev_band_output_dir)
    elev_band_final.to_file(os.path.join(elev_band_output_dir, 'studyArea_elev_band.shp'))

    # Remove temporary directory
    shutil.rmtree(os.path.join(temporary_dir))

    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Elevation Bands Complete!')
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

if config_file_info['generate_elev_bands'] == 'no':
    # set directory/path
    reclassify_dir = os.path.join(temporary_dir, 'elev_band')
    if not os.path.exists(reclassify_dir):
      os.makedirs(reclassify_dir)
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Upload elevation band layer')
    print('-----------------------------------------------------------------------------------------------')
    print(f'drag-and-drop elevation band = file into following folder: {reclassify_dir}')
    response = input("Have you uploaded the elevation band (.shp) file (yes or no): ")
    if response == "yes":
      upload_elev_bands(reclassify_dir, main_dir)
      remove_temp_data(temporary_dir)

elif config_file_info['generate_elev_bands'] == 'yes':
    min_elev = int(config_file_info['min_elev'])
    increment_val = int(config_file_info['increment_val'])

    # set directory/path
    reclassify_dir = os.path.join(temporary_dir, 'elev_band')
    if not os.path.exists(reclassify_dir):
      os.makedirs(reclassify_dir)

    # determine min and max elevation
    dem_max = min_max_elev(reclassify_dir, main_dir)
    # reclassify DEM to generate elevation bands
    reclassify_dem(increment_val, min_elev, dem_max, reclassify_dir, main_dir, temporary_dir)
    # visualize final product
    visualize_elev_bands_func(increment_val, min_elev, temporary_dir, main_dir, dem_max)
    # delete temorary folder
    remove_temp_data(temporary_dir)
"""
