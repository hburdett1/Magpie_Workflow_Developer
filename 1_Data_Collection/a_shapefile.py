print('\n-----------------------------------------------------------------------------------------------')
print('( ) Load libraries')
print('-----------------------------------------------------------------------------------------------')
import pandas as pd
import numpy as np
import os
import json
from IPython.display import display
import geopandas as gpd
import shutil
import matplotlib.pyplot as plt
from glob import glob
import zipfile
import sys
import time
# getting coordinated of study area
from geopy.geocoders import Nominatim
# interactive map
import folium
import branca
import rasterstats as rs
# BasinMaker
from basinmaker import basinmaker
from basinmaker.postprocessing.plotleaflet import plot_routing_product_with_ipyleaflet
from basinmaker.postprocessing.downloadpd import Download_Routing_Product_For_One_Gauge
from basinmaker.postprocessing.downloadpdptspurepy import Download_Routing_Product_From_Points_Or_LatLon


def studyArea_location(city_name, define_lat, define_lon):
    """
    Collects study area location information using either city name or defined coordinates.

    Parameters:
    - city_name (str): The name of the city to determine coordinates.
    - define_lat (str): Defined latitude (if available).
    - define_lon (str): Defined longitude (if available).

    Returns:
    Tuple of latitude and longitude.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Collect study area location information')
    print('-----------------------------------------------------------------------------------------------')

    # Determine coordinates based on city name using geopy.geocoders library
    if city_name != 'NA' and define_lat == 'NA':
        geolocator = Nominatim(user_agent="Magpie")
        location = geolocator.geocode(city_name)
        print(location)
        found_lat, found_lon = location.latitude, location.longitude
        lat, lon = found_lat, found_lon
        print("Study area coordinates:", lat, lon)
        return lat, lon

    # Use defined coordinates
    elif city_name == 'NA' and define_lat != 'NA':
        lat, lon = float(define_lat), float(define_lon)
        print("Study area coordinates:", lat, lon)
        return lat, lon

    # If neither coordinates nor city name are defined, encourage the user to define a gauge name
    elif city_name and define_lat == "NA":
        lat = None
        print("Define gauge name")

def interactive_gauge_map(lat, lon, main_dir):
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Interactive gauge plot')
    print('-----------------------------------------------------------------------------------------------')

    # Read in CSV with gauge information
    gauge_info = pd.read_csv(os.path.join(main_dir, 'extras', 'subbasin_plots', 'obs_gauges_NA_v2-1.csv'))

    # Format fancy folium pop-up
    def fancy_html(row):
        subId_val, obs_gauge, lat_info, lon_info, sub_Reg = (
            gauge_info['SubId'].iloc[row],
            gauge_info['Obs_NM'].iloc[row],
            gauge_info['POINT_Y'].iloc[row],
            gauge_info['POINT_X'].iloc[row],
            gauge_info['Sub_Region'].iloc[row]
        )
        html = f"""<!DOCTYPE html>
            <html>
            <p>SubID: {subId_val}</p>
            <p>Obs Gauge: {obs_gauge}</p>
            <p>Lat: {lat_info}</p>
            <p>Lon: {lon_info}</p>
            <p>Sub Region: {sub_Reg}</p>
            </html>
            """
        return html

    # Generate map with rectangular guide to assist in identifying the ideal subbasin/gauges to use
    if lat is None:
        print("Define gauge name in cell above")
    else:
        grid_pt = (lat, lon)
        W, E, N, S = grid_pt[1] - 0.5, grid_pt[1] + 0.5, grid_pt[0] + 0.5, grid_pt[0] - 0.5
        upper_left, upper_right, lower_right, lower_left = (N, W), (N, E), (S, E), (S, W)
        line_color, fill_color, weight, text = 'red', 'red', 2, 'text'
        edges = [upper_left, upper_right, lower_right, lower_left]

        map_osm = folium.Map(location=[lat, lon], zoom_start=9)
        folium.LatLngPopup().add_to(map_osm)

        for i in range(len(gauge_info)):
            html = fancy_html(i)
            iframe = branca.element.IFrame(html=html, width=200, height=200)
            popup = folium.Popup(iframe, parse_html=True)
            # Adds markers for each gauge
            folium.Marker([gauge_info['POINT_Y'].iloc[i], gauge_info['POINT_X'].iloc[i]], popup=popup).add_to(map_osm)

        # Displays interactive map
        display(map_osm.add_child(folium.vector_layers.Polygon(locations=edges, color=line_color, fill_color=fill_color,
                                                               weight=weight, popup=folium.Popup(text))))

def download_routing_product_lat_lon(lat, lon, product_name):
    """
    Download BasinMaker routing product using latitude and longitude.

    Parameters:
    - lat (float): Latitude coordinate.
    - lon (float): Longitude coordinate.
    - product_name (str): Name of the routing product to be downloaded.

    Returns:
    - str: Path to the downloaded routing product.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Download routing product')
    print('-----------------------------------------------------------------------------------------------')
    print('product_name: ', product_name)

    # Create a DataFrame with the given lat and lon coordinates
    coords = pd.DataFrame({'lat': [lat], 'lon': [lon]})

    # Download routing product using provided coordinates
    subid, product_path = Download_Routing_Product_From_Points_Or_LatLon(
        product_name=product_name, Lat=coords['lat'], Lon=coords['lon']
    )

    print('Successfully downloaded routing product using lat and lon!')
    return product_path

def download_routing_product_gauge(product_name, gauge_name):
    """
    Download BasinMaker routing product using a gauge name.

    Parameters:
    - product_name (str): Name of the routing product to be downloaded.
    - gauge_name (str): Name of the gauge for which the routing product is downloaded.

    Returns:
    - str: Path to the downloaded routing product.
    """
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Download routing product')
    print('-----------------------------------------------------------------------------------------------')

    # Download routing product using the provided gauge name
    subid, product_path = Download_Routing_Product_For_One_Gauge(gauge_name=gauge_name, product_name=product_name)

    print('Successfully downloaded routing product using the gauge name!')
    return product_path


# extracts drainage area
def extract_drainage_area(product_path,most_down_stream_subbasin_ids,
                          most_up_stream_subbasin_ids,temporary_dir,version_num,main_dir):
    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Extract drainage area and simplify drainage product')
    print('-----------------------------------------------------------------------------------------------')
    most_down_stream_subbasin_ids_lst = [most_down_stream_subbasin_ids]
    most_up_stream_subbasin_ids_lst = [most_up_stream_subbasin_ids]

    # define the folder path for downloaded and unziped lake river routing prodcut folder, where several GIS files exist
    unzip_routing_product_folder = product_path

    # define another folder that will save the outputs
    folder_product_for_interested_gauges=os.path.join(temporary_dir,f'catchment_extraction_{most_down_stream_subbasin_ids}')
    # if folder doesn's exist
    if not os.path.exists(folder_product_for_interested_gauges):
        os.makedirs(folder_product_for_interested_gauges)

    # Initialize the basinmaker
    start = time.time()
    bm = basinmaker.postprocess()

    # extract subregion of the routing product
    bm.Select_Subregion_Of_Routing_Structure(
        path_output_folder = folder_product_for_interested_gauges,
        routing_product_folder = unzip_routing_product_folder,
        most_down_stream_subbasin_ids=most_down_stream_subbasin_ids_lst,
        most_up_stream_subbasin_ids=most_up_stream_subbasin_ids_lst,               # -1: extract to the most-upstream (headwater) subbasin; other subbasin ID: extract the areas from the outlet to the provided subbasin.
        gis_platform="purepy",
    )
    end = time.time()
    print("This section took  ", end - start, " seconds")

    # read in study area
    studyArea_bound = gpd.read_file(os.path.join(folder_product_for_interested_gauges,f'catchment_without_merging_lakes_{version_num}.shp'))
    studyArea_bound["dissolve"] = 1

    boundary = studyArea_bound[['dissolve', 'geometry']]
    cont_studyArea = boundary.dissolve(by='dissolve')

    # generate drive folder
    shp_dir = os.path.join(main_dir, 'shapefile')
    shp_path = os.path.isdir(shp_dir)
    if not shp_path:
      os.makedirs(shp_dir)
      print("created  folder: ", shp_dir)

    # save to drive
    cont_studyArea.to_file(os.path.join(shp_dir,'studyArea_outline.shp'))

    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Shapefile complete!')
    print('-----------------------------------------------------------------------------------------------')

def format_shapefile(shp_file_path):
  print('\n-----------------------------------------------------------------------------------------------')
  print('( ) Format Shapefile')
  print('-----------------------------------------------------------------------------------------------')
  # find name of shapefile
  for shp_file in os.listdir(os.path.join(shp_file_path)):
      if shp_file.endswith(".shp"):
        shp_file_name = shp_file

  studyArea_bound = gpd.read_file(os.path.join(shp_file_path,shp_file_name))
  studyArea_bound["dissolve"] = 1

  boundary = studyArea_bound[['dissolve', 'geometry']]
  cont_studyArea = boundary.dissolve(by='dissolve')

  # remove contents in folder
  for f in glob (os.path.join(shp_file_path,'*')):
    os.remove(f)

  cont_studyArea.to_file(os.path.join(main_dir, 'shapefile','studyArea_outline.shp'))

def visualize_shapefile(shp_file_path):
  print('\n-----------------------------------------------------------------------------------------------')
  print('( ) Visualize Shapefile')
  print('-----------------------------------------------------------------------------------------------')

  shp_boundary = gpd.read_file(os.path.join(shp_file_path))

  # check projection
  if shp_boundary.crs !=  'EPSG:4326':
    # reproject
    shp_lyr_crs = shp_boundary.to_crs(epsg=4326)
    print('Shapefile layer has been reprojected to match shapefile')
  else:
    shp_lyr_crs = shp_boundary
    print('Coordinate systems match!')

  # determine the boundary of the provided shapefile
  bounds = shp_lyr_crs.bounds
  west, south, east, north = bounds = bounds.loc[0]
  shp_bounds = [south,west]

  map = folium.Map(location=shp_bounds, zoom_start=10)
  folium.GeoJson(data=shp_boundary["geometry"]).add_to(map)
  display(map)

  print('\n-----------------------------------------------------------------------------------------------')
  print('( ) Shapefile complete!')
  print('-----------------------------------------------------------------------------------------------')

def remove_temp_data(main_dir, temporary_dir,product_path):
  print('\n-----------------------------------------------------------------------------------------------')
  print('( ) Remove unnecessary files')
  print('-----------------------------------------------------------------------------------------------')

  if os.path.exists(temporary_dir):
    shutil.rmtree(temporary_dir)
  if os.path.exists(product_path):
    shutil.rmtree(product_path)
  zip_files_rm = ((glob(os.path.join(product_path+"*.zip"))))
  for files_rm in zip_files_rm:
    os.remove(files_rm)

"""
main_dir = sys.argv[1]

# read in information from configuration file
with open(os.path.join(main_dir,"configuration_file.json"), "r") as f:
    config_file_info = json.load(f)

# define temporary dir
temporary_dir = config_file_info['temporary_dir']

if config_file_info['generate_shapefile'] == 'no':
  # generate drive folder
  shp_dir = os.path.join(main_dir, 'shapefile')
  shp_path = os.path.isdir(shp_dir)
  if not shp_path:
    os.makedirs(shp_dir)
    print("created  folder: ", shp_dir)
  print('\n-----------------------------------------------------------------------------------------------')
  print('Upload Shapefile')
  print('-----------------------------------------------------------------------------------------------')
  print(f'drag-and-drop shapefile into the following folder: {shp_dir}')
  input_response = input("Have you uploaded the study area (.shp) file (yes or no): ")
  if input_response == 'yes':
    shp_file_path = os.path.join(main_dir, 'shapefile')
    # format shapefile
    format_shapefile(shp_file_path)
    # visualize final shapefile
    visualize_shapefile(shp_file_path)

elif config_file_info['generate_shapefile'] == 'yes':
  # define input variables
  city_name = config_file_info['city_name_shp']
  define_lat = config_file_info['define_lat_shp']
  define_lon = config_file_info['define_lon_shp']
  interactive_map_response = config_file_info['interactive_map_response']
  product_name = config_file_info['product_name_shp']
  gauge_name = config_file_info['gauge_name_shp']
  most_down_stream_subbasin_ids = config_file_info['most_down_stream_subbasin_ids_shp'] 
  most_up_stream_subbasin_ids = config_file_info['most_up_stream_subbasin_ids_shp']

  # assign lat and lon to either found or identified coordinates
  if city_name != 'NA' and define_lat == 'NA':
      lat, lon = studyArea_location(city_name,define_lat,define_lon)
  elif city_name == 'NA' and define_lat != 'NA':
      lat,lon = float(define_lat), float(define_lon),
      print("Study area coordinates:", lat,lon)
  elif city_name and define_lat == "NA":
      lat = None

  # generate map with rectangular guide to assist in identifying the ideal subbasin/gauges to use
  if interactive_map_response == 'yes':
      interactive_gauge_map(lat,lon, main_dir)
  else:
    interactive_map_response = 'no'

  # download routing product
  if gauge_name == "NA":
    product_path = download_routing_product_lat_lon(lat,lon,product_name)
  else:
    product_path = download_routing_product_gauge(product_name,gauge_name)

  # extract drainage area
  if product_name == "NALRP":
    version_num = "v2-1"
  if product_name == "OLRP":
    version_num = "v1-0"
  if product_name == "OLRP_V2":
    version_num = "v2-0"

  extract_drainage_area(product_path,most_down_stream_subbasin_ids,
                            most_up_stream_subbasin_ids,temporary_dir,version_num,main_dir)

  # visualize final shapefile
  shp_file_path = os.path.join(main_dir, 'shapefile')
  visualize_shapefile(shp_file_path)

  # delete temorary folder
  remove_temp_data(main_dir, temporary_dir,product_path)
  """
  

