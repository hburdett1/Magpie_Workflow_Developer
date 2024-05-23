import pandas as pd
import numpy as np
import os
import sys
import subprocess
import json
import geopandas as gpd
from IPython.display import display
import shutil
import matplotlib.pyplot as plt
from glob import glob
import wget
from pathlib import Path
# getting coordinated of study area
from geopy.geocoders import Nominatim
# BasinMaker
from basinmaker import basinmaker
from basinmaker.postprocessing.plotleaflet import plot_routing_product_with_ipyleaflet
from basinmaker.postprocessing.downloadpd import Download_Routing_Product_For_One_Gauge
from basinmaker.postprocessing.downloadpdptspurepy import Download_Routing_Product_From_Points_Or_LatLon
import time
import rasterstats as rs

def download_routing_product(product_name, gauge_name, define_lat, define_lon):
    """
    Download routing product based on provided information.

    Parameters:
    - product_name (str): Name of the routing product.
    - gauge_name (str): Gauge name for downloading by gauge.
    - city_name (str): City name for downloading by city coordinates.
    - define_lat (str): Latitude for downloading by specified coordinates.
    - define_lon (str): Longitude for downloading by specified coordinates.

    Returns:
    - product_path (str): Path to the downloaded routing product.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Download routing product')
    print('-----------------------------------------------------------------------------------------------')

    if gauge_name != "NA":
        subid, product_path = Download_Routing_Product_For_One_Gauge(gauge_name=gauge_name, product_name=product_name)
        print('Successfully downloaded routing product using the gauge name!')
    elif define_lat != 'NA':
        lat, lon = float(define_lat), float(define_lon)
        print("Study area coordinates:", lat, lon)
        coords = pd.DataFrame({'lat': [lat], 'lon': [lon]})
        subid, product_path = Download_Routing_Product_From_Points_Or_LatLon(product_name=product_name,
                                                                             Lat=coords['lat'], Lon=coords['lon'])
        print('Successfully downloaded routing product using the specified coordinates!')

    return product_path

def extract_drainage_area(product_path, subid_of_interested_gauges, most_up_stream_subbasin_ids, temporary_dir):
    """
    Extract drainage area based on provided subbasin IDs.

    Parameters:
    - product_path (str): Path to the downloaded and unzipped lake-river routing product folder.
    - subid_of_interested_gauges (str): Subbasin ID where the gauge is situated.
    - most_up_stream_subbasin_ids (str): Most upstream subbasin IDs for extraction.

    Returns:
    - folder_product_for_interested_gauges (str): Path to the folder containing the extracted drainage area.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Extract drainage area')
    print('-----------------------------------------------------------------------------------------------')

    # Define subbasin ID lists
    subid_of_interested_gauges_lst = [subid_of_interested_gauges]
    most_up_stream_subbasin_ids_lst = [most_up_stream_subbasin_ids]

    # Define the folder path for downloaded and unzipped lake-river routing product folder
    unzip_routing_product_folder = product_path

    # Define another folder that will save the outputs
    folder_product_for_interested_gauges = os.path.join(temporary_dir, 'catchment_extraction')

    # Initialize the basinmaker
    start = time.time()
    bm = basinmaker.postprocess()

    # Extract subregion of the routing product
    bm.Select_Subregion_Of_Routing_Structure(
        path_output_folder=folder_product_for_interested_gauges,
        routing_product_folder=unzip_routing_product_folder,
        most_down_stream_subbasin_ids=subid_of_interested_gauges_lst,
        most_up_stream_subbasin_ids=most_up_stream_subbasin_ids_lst,
        gis_platform="purepy",
    )
    end = time.time()
    print("This section took ", end - start, " seconds")

    return folder_product_for_interested_gauges

def remove_small_lakes(lake_size, temporary_dir):
    """
    Remove small lakes based on the provided lake size threshold.

    Parameters:
    - lake_size (float): Lake size threshold (unit: km^2).

    Returns:
    - folder_product_after_filter_lakes (str): Path to the folder containing the drainage network after removing small lakes.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Remove small lakes')
    print('-----------------------------------------------------------------------------------------------')

    # Print all lake ID file paths in the temporary directory
    if os.path.exists(os.path.join(temporary_dir, 'catchment_extraction*', 'sl_connected_lake_v*.shp')):
      for lake_ID_file_path in glob(os.path.join(temporary_dir, 'catchment_extraction*', 'sl_connected_lake_v*.shp')):
          print(lake_ID_file_path)
          # Read the lake ID file
          lake_ID_file = gpd.read_file(lake_ID_file_path)

          # Filter lakes based on size
          small_lakes = lake_ID_file[lake_ID_file["Lake_area"] > lake_size]

          # Extract lake IDs
          lake_list = (small_lakes['Hylak_id'].unique()).tolist()
          lake_IDs = [i for i in lake_list if i != 0]
          print("Lake IDs: ", lake_IDs)

          # Define the output folder for interested gauges
          folder_product_for_interested_gauges = os.path.join(temporary_dir, 'catchment_extraction')

          # Define the input product folder path, which is the output folder of the previous section
          input_routing_product_folder = folder_product_for_interested_gauges

          # Define a list containing HyLakeId IDs of lakes of interest
          interested_lake_ids = lake_IDs

          # Update the variable inside the if block
          folder_product_after_filter_lakes = os.path.join(temporary_dir, 'filter_lakes')

          start = time.time()

          # Remove small lakes using BasinMaker
          bm = basinmaker.postprocess()

          bm.Remove_Small_Lakes(
              path_output_folder=folder_product_after_filter_lakes,
              routing_product_folder=input_routing_product_folder,
              connected_lake_area_thresthold=lake_size,
              non_connected_lake_area_thresthold=lake_size,
              selected_lake_ids=interested_lake_ids,
              gis_platform="purepy",
          )
          end = time.time()
          print("This section took  ", end - start, " seconds")
    else:
      # Define the output folder for interested gauges
      folder_product_for_interested_gauges = os.path.join(temporary_dir, 'catchment_extraction')
      folder_product_after_filter_lakes = folder_product_for_interested_gauges

    return folder_product_after_filter_lakes

def simplify_drainage_area(minimum_subbasin_drainage_area, temporary_dir):
    """
    Simplify the drainage product by increasing the size of subbasins.

    Parameters:
    - minimum_subbasin_drainage_area (float): Minimum drainage area of subbasins (unit: km^2).

    Returns:
    - folder_product_after_increase_catchment_drainage_area (str): Path to the folder containing the drainage network after simplification.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Simplify drainage product')
    print('-----------------------------------------------------------------------------------------------')

    # Print information about the process
    print(f"Simplifying drainage product with a minimum subbasin drainage area of {minimum_subbasin_drainage_area} km^2")

    # Define the input folder path, which is the output folder of the previous section
    if os.path.exists(os.path.join(temporary_dir, 'catchment_extraction*', 'sl_connected_lake_v*.shp')):
      folder_product_after_filter_lakes = os.path.join(temporary_dir, 'filter_lakes')
    else:
      folder_product_after_filter_lakes = os.path.join(temporary_dir, 'catchment_extraction')
    input_routing_product_folder = folder_product_after_filter_lakes

    # Define the output folder after increasing catchment drainage area
    folder_product_after_increase_catchment_drainage_area = os.path.join(temporary_dir, 'drainage_area')

    # Initialize the basinmaker
    start = time.time()
    bm = basinmaker.postprocess()

    # Remove river reaches and increase the size of subbasins
    bm.Decrease_River_Network_Resolution(
        path_output_folder=folder_product_after_increase_catchment_drainage_area,
        routing_product_folder=input_routing_product_folder,
        minimum_subbasin_drainage_area=minimum_subbasin_drainage_area,
        gis_platform="purepy",
    )
    end = time.time()
    print("This section took  ", end - start, " seconds")

    return folder_product_after_increase_catchment_drainage_area

def save_routing_product(version_number, temporary_dir, main_dir):
    """
    Save the simplified routing product to the drive.

    Parameters:
    - version_number (str): Version number of the routing product.

    Returns:
    - None
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Save routing product')
    print('-----------------------------------------------------------------------------------------------')

    # Define the routing directory
    routing_dir = os.path.join(main_dir, 'workflow_outputs', 'routing_product')

    # Create the routing directory if it doesn't exist
    if not os.path.exists(routing_dir):
        os.makedirs(routing_dir)
        print("created folder:", routing_dir)

    # Define paths for different routing product files
    finalcat_info_riv_path = Path(os.path.join(temporary_dir, 'drainage_area', f'finalcat_info_riv_{version_number}.shp'))
    finalcat_info_path = Path(os.path.join(temporary_dir, 'drainage_area', f'finalcat_info_{version_number}.shp'))
    sl_connected_lake_path = Path(os.path.join(temporary_dir, 'drainage_area', f'sl_connected_lake_{version_number}.shp'))
    sl_non_connected_lake_path = Path(os.path.join(temporary_dir, 'drainage_area', f'sl_non_connected_lake_{version_number}.shp'))

    # Helper function to reproject and save GeoDataFrame
    def reproject_and_save(gdf, output_path):
        if gdf.exists():
            gdf_data = gpd.read_file(gdf)
            # Reproject if needed
            if gdf_data.crs != 'EPSG:4326':
                gdf_data = gdf_data.to_crs(epsg=4326)
            # Save the GeoDataFrame
            gdf_data.to_file(output_path)

    # Reproject and save each routing product file
    reproject_and_save(finalcat_info_riv_path, os.path.join(routing_dir, f'finalcat_info_riv_{version_number}.shp'))
    reproject_and_save(finalcat_info_path, os.path.join(routing_dir, f'finalcat_info_{version_number}.shp'))
    reproject_and_save(sl_connected_lake_path, os.path.join(routing_dir, f'sl_connected_lake_{version_number}.shp'))
    reproject_and_save(sl_non_connected_lake_path, os.path.join(routing_dir, f'sl_non_connected_lake_{version_number}.shp'))

def remove_temp_data(product_path,temporary_dir):
    if os.path.exists(temporary_dir):
      shutil.rmtree(temporary_dir)
    if os.path.exists(product_path):
      shutil.rmtree(product_path)
    zip_files_rm = ((glob(os.path.join(product_path+"*.zip"))))
    for files_rm in zip_files_rm:
      os.remove(files_rm)

"""
# define main dir
main_dir = sys.argv[1]

# read in information from configuration file
with open(os.path.join(main_dir,"configuration_file.json"), "r") as f:
    config_file_info = json.load(f)

# define temporary dir
temporary_dir = config_file_info['temporary_dir']

if config_file_info['upload_routing_product'] == 'yes':
    # define input variables
    version_num = config_file_info['version_num']
    product_path = 'None'

    routing_temp_dir = os.path.join(main_dir, 'workflow_outputs', 'routing_product')
    if not os.path.exists(routing_temp_dir):
        os.makedirs(routing_temp_dir)

    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Upload routing product')
    print('-----------------------------------------------------------------------------------------------')
    print(f'drag-and-drop routing product files into following folder: {routing_temp_dir}')
    response = input("Have you uploaded the related subbasin shapefiles (.shp) file (yes or no): ")

    if response.lower() == "yes":
        # -- Visualize --
        print('\n-----------------------------------------------------------------------------------------------')
        print('( ) Visualize Simplified Routing Product')
        print('-----------------------------------------------------------------------------------------------')
        # define routing number
        routing_product_version_number = version_num
        # plot product
        display(plot_routing_product_with_ipyleaflet(path_to_product_folder=routing_temp_dir, version_number=routing_product_version_number))

        # Remove temporary data after visualization
        remove_temp_data(product_path,temporary_dir)

if config_file_info['upload_routing_product'] == 'no':
    product_name = config_file_info['product_name']
    version_number = config_file_info['version_num']
    gauge_name = config_file_info['gauge_name'] 
    define_lat = config_file_info['define_lat']
    define_lon = config_file_info['define_lon']
    subid_of_interested_gauges = config_file_info['most_down_stream_subbasin_ids']
    most_up_stream_subbasin_ids = config_file_info['most_up_stream_subbasin_ids']
    lake_size = int(config_file_info['lake_size'])
    minimum_subbasin_drainage_area = config_file_info['minimum_subbasin_drainage_area']

    # download routing product
    product_path = download_routing_product(product_name, gauge_name, define_lat, define_lon)
    # extract drainage area 
    folder_product_for_interested_gauges = extract_drainage_area(product_path, subid_of_interested_gauges, most_up_stream_subbasin_ids, temporary_dir)
    # define the path to the routing product folder
    path_to_input_routing_product_folder = folder_product_for_interested_gauges
    # plot product
    display(plot_routing_product_with_ipyleaflet(path_to_product_folder = path_to_input_routing_product_folder,version_number = version_number))

    # remove small lakes
    folder_product_after_filter_lakes_derived = remove_small_lakes(lake_size, temporary_dir)
    # Define the path to the routing product folder
    path_to_input_routing_product_folder = folder_product_after_filter_lakes_derived 
    # plot product
    display(plot_routing_product_with_ipyleaflet(path_to_product_folder = path_to_input_routing_product_folder,version_number = version_number))

    # run function to simplify drainage basin
    folder_product_after_increase_catchment_drainage_area = simplify_drainage_area(minimum_subbasin_drainage_area, temporary_dir)
    # plot product
    display(plot_routing_product_with_ipyleaflet(path_to_product_folder = folder_product_after_increase_catchment_drainage_area,version_number = version_number))

    # save routing product
    save_routing_product(version_number, temporary_dir, main_dir)
    # delete temorary folder
    remove_temp_data(product_path,temporary_dir)
"""