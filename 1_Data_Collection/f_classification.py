import sys
import os
import json
import pandas as pd

def landcover_classification(main_dir, landuse_ID, landuse_Classifications):
    """
    Classify landcover based on provided information.

    Parameters:
    - main_dir (str): Main directory path.
    - landuse_ID (str): Comma-separated string of landuse IDs.
    - landuse_Classifications (str): Comma-separated string of landuse classifications.
    """

    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Landcover classification')
    print('-----------------------------------------------------------------------------------------------')

    land_dir = os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'Landcover')

    if not os.path.exists(land_dir):
        os.makedirs(land_dir)

    if landuse_ID == "NA":
        print('Please upload landcover csv file to: ', land_dir)
    else:
        landID_val = list(landuse_ID.split(","))
        landclass_val = list(landuse_Classifications.split(","))

        land_class = {
            'Landuse_ID': landID_val,
            'LAND_USE_C': landclass_val
        }

        land_class_df = pd.DataFrame(land_class)
        land_class_df['Landuse_ID'] = land_class_df['Landuse_ID'].astype(str)

        print(land_class_df)

        land_class_df.to_csv(os.path.join(land_dir, 'landcover_info.csv'), index=False)

def vegetation_classification(main_dir, veg_ID, veg_Classifications):
    """
    Classify vegetation based on provided information.

    Parameters:
    - main_dir (str): Main directory path.
    - veg_ID (str): Comma-separated string of vegetation IDs.
    - veg_Classifications (str): Comma-separated string of vegetation classifications.
    """

    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Vegetation classification')
    print('-----------------------------------------------------------------------------------------------')

    veg_dir = os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'Landcover')

    if not os.path.exists(veg_dir):
        os.makedirs(veg_dir)

    if veg_ID == "NA":
        print('Please upload vegetation csv file to: ', veg_dir)
    else:
        vegID_val = list(veg_ID.split(","))
        vegClass_val = list(veg_Classifications.split(","))

        veg_class = {
            'Veg_ID':  vegID_val,
            'VEG_C': vegClass_val
        }

        veg_class_df = pd.DataFrame(veg_class)
        veg_class_df['Veg_ID'] = veg_class_df['Veg_ID'].astype(str)

        print(veg_class_df)

        veg_class_df.to_csv(os.path.join(veg_dir, 'veg_info.csv'), index=False)

def soil_classification(main_dir, soil_ID, soil_Classifications):
    """
    Classify soil based on provided information.

    Parameters:
    - main_dir (str): Main directory path.
    - soil_ID (str): Comma-separated string of soil IDs.
    - soil_Classifications (str): Comma-separated string of soil classifications.
    """

    print('\n-----------------------------------------------------------------------------------------------')
    print('( ) Soil classification')
    print('-----------------------------------------------------------------------------------------------')

    soil_dir = os.path.join(main_dir, 'workflow_outputs', '1_HRU_data', 'Soil')

    if not os.path.exists(soil_dir):
        os.makedirs(soil_dir)

    if soil_ID == "NA":
        print('Please upload soil csv file to: ', soil_dir)
    else:
        soilID_val = list(soil_ID.split(","))
        soilClass_val = list(soil_Classifications.split(","))

        soil_class = {
            'Soil_ID':  soilID_val,
            'SOIL_PROF': soilClass_val
        }

        soil_class_df = pd.DataFrame(soil_class)
        soil_class_df.to_csv(os.path.join(soil_dir, 'soil_info.csv'), index=False)

        print(soil_class_df)

"""
main_dir = sys.argv[1]

# read in information from configuration file
with open(os.path.join(main_dir,"configuration_file.json"), "r") as f:
    config_file_info = json.load(f)

# define land cover IDS
landuse_ids = config_file_info['landuse_ID']
# define landcover classification
landuse_classifications = config_file_info['landuse_Classifications']
# run function
landcover_classification(main_dir, landuse_ids, landuse_classifications)

# define veg IDS
veg_ids = config_file_info['veg_ID']
# define veg classification
veg_classifications = config_file_info['veg_Classifications']
# run function
vegetation_classification(main_dir, veg_ids, veg_classifications)

# define soil IDS
soil_ids = config_file_info['soil_ID']
# define soil classification
soil_classifications = config_file_info['soil_Classifications']
# run function
soil_classification(main_dir, soil_ids, soil_classifications)
"""




