import pandas as pd
import numpy as np
import os
import sys
import json
import geopandas as gpd
import shutil
import matplotlib.pyplot as plt
from IPython.display import display
from glob import glob
#import wget
from pathlib import Path
# getting coordinated of study area
from geopy.geocoders import Nominatim
# BasinMaker
from basinmaker import basinmaker
import time


def define_hrus(landcover, soil, elevation_bands, dem, aspect,
                area_ratio_thresholds, pixel_size, version_number, main_dir):
    """
    Define Hydrologic Response Units (HRUs) using BasinMaker.

    Parameters:
    - landcover (str): Flag indicating whether to include landcover data.
    - soil (str): Flag indicating whether to include soil data.
    - elevation_bands (str): Flag indicating whether to include elevation band data.
    - dem (str): Flag indicating whether to include Digital Elevation Model (DEM) data.
    - aspect (str): Flag indicating whether to include aspect data.
    - area_ratio_thresholds (list): List of area ratio thresholds.
    - pixel_size (float): Pixel size for processing.
    - version_number (str): Version number for the output files.
    """

    # Define working directory paths
    hru_dir = os.path.join(main_dir, 'workflow_outputs', '1_HRU_data')
    routing_dir = os.path.join(main_dir, 'workflow_outputs', 'routing_product')

    # Generate maps folder
    HRU_output_folder = os.path.join(main_dir, 'workflow_outputs', 'RavenInput', 'maps')

    # Create the HRU output folder if it doesn't exist
    if not os.path.isdir(HRU_output_folder):
        os.makedirs(HRU_output_folder)
        print("created folder:", HRU_output_folder)

    delin_list = []

    # Paths for DEM data
    if dem == 'yes':
      for dem_file in os.listdir(os.path.join(hru_dir, 'DEM')):
        if dem_file.endswith(".tif"):
          dem_file_name = dem_file
      # Adjust the following paths based on your data locations
      path_to_dem = os.path.join(hru_dir, 'DEM',dem_file_name)
      print(path_to_dem)
    else:
      path_to_dem = '#'

    # Paths for Elevation Bands data
    if elevation_bands == 'yes':
      for elev_file in os.listdir(os.path.join(hru_dir, 'Elevation_band')):
        if elev_file.endswith(".shp"):
          elev_file_name = elev_file
      # Adjust the following paths based on your data locations
      path_other_polygon_1 = os.path.join(hru_dir, 'Elevation_band',elev_file_name)
      other_polygon = gpd.read_file(path_other_polygon_1)
      print(path_other_polygon_1)
      delin_list.append('O_ID_1')
    else:
      path_other_polygon_1 = '#'

    # Paths for Aspect data
    if aspect == 'yes':
        for aspect_file in os.listdir(os.path.join(hru_dir, 'Aspect')):
          if aspect_file.endswith(".shp"):
            aspect_file_name = aspect_file
        # Adjust the following paths based on your data locations
        path_to_aspect = os.path.join(hru_dir, 'Aspect',aspect_file_name)
        print(path_to_aspect)
    else:
      path_to_aspect = '#'

    # Paths for Landcover data
    if landcover == 'yes':
      for landcov_file in os.listdir(os.path.join(hru_dir, 'Landcover')):
        if landcov_file.endswith(".shp"):
          landcov_file_name = landcov_file
      # Adjust the following paths based on your data locations
      path_landuse_polygon = os.path.join(hru_dir,'Landcover',landcov_file_name)
      path_landuse_info = os.path.join(hru_dir,'Landcover','landcover_info.csv')
      path_veg_info = os.path.join(hru_dir, 'Landcover', "veg_info.csv")
      #landuse = gpd.read_file(path_landuse_polygon)
      landuse_info = pd.read_csv(path_landuse_info)
      print(path_landuse_polygon)
      delin_list.append('Landuse_ID')
      delin_list.append('Veg_ID')
    else:
      path_landuse_polygon = '#'
      path_landuse_info = os.path.join(hru_dir,'Landcover','landcover_info.csv')
      path_veg_info = os.path.join(hru_dir, 'Landcover', "veg_info.csv")

    # Paths for Soil data
    if soil == 'yes':
      for soil_file in os.listdir(os.path.join(hru_dir, 'Soil')):
        if soil_file.endswith(".shp"):
          soil_file_name = soil_file
      # Adjust the following paths based on your data locations
      path_soil_polygon = os.path.join(hru_dir,'Soil',soil_file_name)
      path_soil_info = os.path.join(hru_dir,'Soil',"soil_info.csv")
      #soil = gpd.read_file(path_soil_polygon)
      soil_info = pd.read_csv(path_soil_info)
      print(path_soil_polygon)
      delin_list.append('Soil_ID')
    else:
      path_soil_polygon = '#'
      path_soil_info = os.path.join(hru_dir,'Soil',"soil_info.csv")

    # Define the input folder path
    input_routing_product_folder = routing_dir

    # Define paths for connected and non-connected lake polygons
    path_connect_lake_polygon_dir = Path(os.path.join(routing_dir, f"sl_connected_lake_{version_number}.shp"))
    path_non_connect_lake_polygon_dir = Path(os.path.join(routing_dir, f"sl_non_connected_lake_{version_number}.shp"))

    # Check if connected lake polygon exists
    path_connect_lake_polygon = path_connect_lake_polygon_dir if path_connect_lake_polygon_dir.exists() else '#'

    # Check if non-connected lake polygon exists
    path_non_connect_lake_polygon = path_non_connect_lake_polygon_dir if path_non_connect_lake_polygon_dir.exists() else '#'

    # Run BasinMaker
    bm = basinmaker.postprocess()
    start = time.time()

    bm.Generate_HRUs(
        path_output_folder=HRU_output_folder,
        path_subbasin_polygon=os.path.join(routing_dir, f"finalcat_info_{version_number}.shp"),
        path_connect_lake_polygon=path_connect_lake_polygon,
        path_non_connect_lake_polygon=path_non_connect_lake_polygon,
        path_landuse_polygon=path_landuse_polygon,
        path_soil_polygon=path_soil_polygon,
        path_other_polygon_1=path_other_polygon_1,
        path_other_polygon_2=path_to_aspect,
        path_landuse_info=path_landuse_info,
        path_soil_info=path_soil_info,
        path_veg_info=path_veg_info,
        path_to_dem=path_to_dem,
        area_ratio_thresholds=area_ratio_thresholds,
        gis_platform="purepy",
        projected_epsg_code='EPSG:3161',
        pixel_size=pixel_size
    )

    end = time.time()
    print("This section took ", end - start, " seconds\n")

    # Visualize HRUs
    hru_polygon = gpd.read_file(os.path.join(HRU_output_folder, "finalcat_hru_info.shp"))

    # Define paths for connected and non-connected lake polygons
    path_connect_lake_polygon_dir = Path(os.path.join(routing_dir, f"sl_connected_lake_{version_number}.shp"))
    path_non_connect_lake_polygon_dir = Path(os.path.join(routing_dir, f"sl_non_connected_lake_{version_number}.shp"))

    ax = hru_polygon.plot(linewidth=1, edgecolor='black', facecolor="none", zorder=0, figsize=(11, 12))

    # Plot connected lake polygons if exists
    if path_connect_lake_polygon_dir.exists():
        sl_lake_ply = gpd.read_file(path_connect_lake_polygon_dir).to_crs(hru_polygon.crs)
        sl_lake_ply.plot(ax=ax, linewidth=0.00001, edgecolor='black', alpha=0.6, zorder=1)

    # Plot non-connected lake polygons if exists
    if path_non_connect_lake_polygon_dir.exists():
        nsl_lake_ply = gpd.read_file(path_non_connect_lake_polygon_dir).to_crs(hru_polygon.crs)
        nsl_lake_ply.plot(ax=ax, linewidth=0.00001, edgecolor='black', alpha=0.6, zorder=1)

    plt.show()

"""
# define main dir
main_dir = sys.argv[1]

# read in information from configuration file
with open(os.path.join(main_dir,"configuration_file.json"), "r") as f:
    config_file_info = json.load(f)

# define temporary dir
temporary_dir = config_file_info['temporary_dir']

if config_file_info['upload_HRU_files'] == 'yes':
    # define version number
    version_number = config_file_info['version_num_hru']

    # define working directory path
    HRU_output_folder = os.path.join(main_dir, 'workflow_outputs', 'RavenInput', 'maps')
    if not os.path.exists(HRU_output_folder):
      os.makedirs(HRU_output_folder)

    routing_dir = os.path.join(main_dir,'workflow_outputs','routing_product')

    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Upload HRUs')
    print('-----------------------------------------------------------------------------------------------')
    print(f'drag-and-drop HRU files into following folder: {HRU_output_folder}')
    response = input("Have you uploaded the related HRU shapefiles (.shp) file (yes or no): ")
    if response == "yes":
      # -- Visualize --
      print('\n-----------------------------------------------------------------------------------------------')
      print('( ) Visualize Uploaded HRUs')
      print('-----------------------------------------------------------------------------------------------')
      hru_polygon = gpd.read_file(os.path.join(HRU_output_folder, "finalcat_hru_info.shp"))
      # define paths
      path_connect_lake_polygon_dir = Path(os.path.join(routing_dir, "sl_connected_lake_"+version_number+".shp"))
      path_non_connect_lake_polygon_dir = Path(os.path.join(routing_dir, "sl_non_connected_lake_"+version_number+".shp"))

      ax = hru_polygon.plot(linewidth = 1,edgecolor='black',facecolor="none",zorder=0,figsize=(11,12))

      if path_connect_lake_polygon_dir.exists():
        sl_lake_ply = gpd.read_file(os.path.join(routing_dir, "sl_connected_lake_"+version_number+".shp")).to_crs(hru_polygon.crs)
        sl_lake_ply.plot(ax = ax, linewidth = 0.00001,edgecolor='black',alpha=0.6,zorder=1)

      if path_non_connect_lake_polygon_dir.exists():
        nsl_lake_ply = gpd.read_file(os.path.join(routing_dir, "sl_non_connected_lake_"+version_number+".shp")).to_crs(hru_polygon.crs)
        nsl_lake_ply.plot(ax = ax, linewidth = 0.00001,edgecolor='black',alpha=0.6,zorder=1)

if config_file_info['upload_HRU_files'] == 'no':
    landcover = config_file_info['landcover_hru']
    soil = config_file_info['soil_hru']
    elevation_bands = config_file_info['elevation_bands_hru']
    dem = config_file_info['dem_hru']
    aspect = config_file_info['aspect_hru']

    area_ratio_thresholds = config_file_info['area_ratio_thresholds']
    pixel_size = int(config_file_info['pixel_size'])
    product_name = config_file_info['product_name_hru']
    version_number = config_file_info['version_num_hru']

    define_hrus(landcover,soil,elevation_bands,dem,aspect,area_ratio_thresholds,pixel_size,version_number)
"""