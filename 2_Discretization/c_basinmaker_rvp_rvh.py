import pandas as pd
import numpy as np
import os
import sys
import json
import geopandas as gpd
import shutil
import matplotlib.pyplot as plt
from glob import glob
#import wget
from pathlib import Path
# getting coordinated of study area
from geopy.geocoders import Nominatim
# BasinMaker
from basinmaker import basinmaker
from basinmaker.postprocessing.plotleaflet import plot_routing_product_with_ipyleaflet
from basinmaker.postprocessing.downloadpd import Download_Routing_Product_For_One_Gauge
from basinmaker.postprocessing.downloadpdptspurepy import Download_Routing_Product_From_Points_Or_LatLon
import time

def rvp_rvh_generate(model_name,main_dir,temporary_dir):
  # Generate Raven RVP and RVH input files
  # save final HRU shapefile to Geojson format for RavenView
  HRU_output_folder = os.path.join(main_dir, 'workflow_outputs', 'RavenInput', 'maps')
  hru_polygon = gpd.read_file(os.path.join(HRU_output_folder, "finalcat_hru_info.shp"))
  hru_polygon.to_file(os.path.join(main_dir, 'shapefile', 'myshpfile.geojson'), driver='GeoJSON')

  # generate Raven files
  bm = basinmaker.postprocess()

  bm.Generate_Raven_Model_Inputs(
      path_hru_polygon         = os.path.join(main_dir, 'workflow_outputs', 'RavenInput', 'maps', 'finalcat_hru_info.shp'),
      model_name            = model_name,                         # This is used for naming the output files.
      subbasingroup_names_channel   =["Allsubbasins"],                        # A subbasin group will be created in the rvh file for simultaneous manipulation in Raven modeling.
      subbasingroup_length_channel   =[-1],
      subbasingroup_name_lake      =["AllLakesubbasins"],
      subbasingroup_area_lake      =[-1],
      path_output_folder         = os.path.join(temporary_dir), # define temporary folder
      aspect_from_gis          = 'purepy',
  )

  # save to drive
  for file in glob(os.path.join(temporary_dir, 'RavenInput','*')):
      shutil.move(file, os.path.join(main_dir, 'workflow_outputs', 'RavenInput'))

# remove temporary data
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

if config_file_info['upload_RVP_RVH'] == 'yes':
    # upload RVH file
    rvp_rvh_temp_dir = os.path.join(main_dir,'workflow_outputs','RavenInput')
    if not os.path.exists(rvp_rvh_temp_dir):
      os.makedirs(rvp_rvh_temp_dir)
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Upload RVP and RVH files')
    print('-----------------------------------------------------------------------------------------------')
    print(f'drag-and-drop routing product files into following folder: {rvp_rvh_temp_dir}')

if config_file_info['upload_RVP_RVH'] == 'no':
    model_name = config_file_info['model_name']
    # generate RVH file
    rvp_rvh_generate(model_name,main_dir,temporary_dir)
    # remove temporary data
    remove_temp_data(temporary_dir)
"""