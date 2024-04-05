# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 14:17:13 2024

@author: jpadil10
"""
import pandas as pd 
import chardet



# # file_2_merge = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Membrane/Membrane_SC_stat.csv'
# file_2_merge = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Membrane/Membrane_SC_calibrated.csv'
# weather_file = "C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/South Creek/MIMS/Surface Water/Atmosphere Data/weather_filtered.csv"
# sonde_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Manta/Manta/HAWKS_4_SC_modified.csv'
# ctd_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/Raw Logger data/CTD_Divers/House_A_Rickabys_SthCreek/CSV/Atmoshperic pressure_231113152722_V9769 - Copy.csv'


# =============================================================================
#  currently calibrate file works with timestamp as the time header and membrane file with time F
# =============================================================================


def merge_ancilliary(input_file, 
                     weather_file = None, weather_header = None,
                     sonde_file = None, sonde_header = None,
                     ctd_file = None, ctd_header = None, ctd_range = None, 
                     show_plots=False):
    
    merged_data = pd.read_csv(input_file)
    # merged_data.sort_values(['timestamp', 'mz'], ascending = True, inplace = True)
    dfs_to_concat = []
    
    if weather_file is not None:
        weather = read_weather(input_file, weather_file, weather_header)
        dfs_to_concat.append(weather)
    
    if sonde_file is not None:
        sonde = read_sonde(input_file, sonde_file)
        dfs_to_concat.append(sonde)
    
    if ctd_file is not None:
        ctd = read_CTD(input_file, ctd_file, ctd_header, ctd_range)
        dfs_to_concat.append(ctd)
    
    if len(dfs_to_concat) > 0:
        merged_data = pd.concat([merged_data] + dfs_to_concat, axis=1, join='outer')
        
    merged_data = merged_data.loc[:, ~merged_data.columns.duplicated()]
    
    print(f'\nSuccesfully merged: {len(dfs_to_concat)} ancilliary files\n')
    
    merged_data.to_csv(input_file.replace('stat.csv', 'calibrated_ancilliary.csv'))

    return merged_data


# =============================================================================
# Function to read, modify (if needed) and match weather data with mims time
# csv. Input file was pre modified to accomodate different dates, separate numbers from text
# and convert string to numbers in excel this should be done in python for future
# iterations

# Inputs:
#     - input_file: membrane stat file to attach all necessary ancilliary data for post-processing
#     - weather_file: raw csv containing atmospheric pressure
#     - time_header: name of time column (case sensitive)
#     - separation: in the case of file containing units attached to the value

# Output: 
#     - will merge weather/atmospheric pressure data to the membrane input file
# =============================================================================

def read_weather(input_file, weather_file, time_header, separation=False, input_header = 'timestamp',):
    
    print("\nCompiling Weather data..")

    data = pd.read_csv(input_file)
    data[input_header] = pd.to_datetime(data[input_header]) 
   
    weather_data = read_csv_with_encoding_and_retry(weather_file)
    weather_data[time_header] = pd.to_datetime(weather_data[time_header])

    # This section is to be used incase file is a copy paste from a source with strings attached like wunderground 
    if separation:
        weather_data['Atmospheric Pressure'] = weather_data['Atmospheric Pressure'].astype(str)
        weather_data[['Atmospheric Pressure','unit']] = weather_data['Atmospheric Pressure'].str.split(n=1, expand=True) # expand= True allocates the other split portion to another column 
        weather_data['Atmospheric Pressure'] = weather_data['Atmospheric Pressure'].str.replace(',','').fillna('')
        weather_data['Atmospheric Pressure'] = pd.to_numeric(weather_data['Atmospheric Pressure'],errors='coerce')
        weather_data['Atmospheric Pressure'].fillna(0, inplace= True)

        merged_bytime = pd.merge_asof(data, weather_data[['Time','Atmospheric Pressure']],
                                      left_on=input_header,right_on='Time', direction='nearest')
        merged_bytime = merged_bytime.drop(columns = 'Time')
   
    # merge_asof uses the nearest number based on the columns to match from the selected dataframes 
    # (like v lookup in excel)
    
    else:
        merged_bytime = pd.merge_asof(data, weather_data[['Time','Atmospheric Pressure']],
                                      left_on=input_header,right_on= time_header, direction='nearest')
        merged_bytime = merged_bytime.drop(columns = time_header)
        
    print("\nWeather data successfully merged")

    return merged_bytime



# =============================================================================
# Function to read, modify (if needed) and match sonde/logger data with mims time csv.
# Raw file was premodified to eliminate rows of documentation data.
# Below merges separate Date and Time columns but also works with stand alone 
# time_header parameter should the Date and Time columns not be present. 

# Inputs:
#     - input_file: membrane stat file to attach all necessary ancilliary data for post-processing
#     - sonde_file: raw csv containing physiochemical parameters
#     - time_header: name of time column (case sensitive)

# Output: 
#     - will merge physiochemical data to the membrane input file
# =============================================================================

def read_sonde(input_file, sonde_file, time_header = None):
    
    print("\nCompiling Sonde data..")

    data = pd.read_csv(input_file)
    data['timestamp'] = pd.to_datetime(data['timestamp'])    
    
    sonde_data = read_csv_with_encoding_and_retry(sonde_file)
    
    # if 'DATE' not in sonde_data.columns or 'TIME' not in sonde_data.columns:
    if time_header is not None:    
        sonde_data[time_header].fillna(0, inplace= True)
        sonde_data[time_header] = pd.to_datetime(sonde_data[time_header])
        sonde_data.drop(columns = ['Depth_m','Int_Batt_V' ], inplace= True)
        
        merged_bytime = pd.merge_asof(data, sonde_data, left_on='timestamp',right_on= time_header, direction='nearest')
        merged_bytime = merged_bytime.drop(columns = time_header)
        
    else:
        sonde_data['DATE'].fillna(0, inplace= True)
        sonde_data['TIME'].fillna(0, inplace= True)
        sonde_data['Timestamp'] = sonde_data['DATE'] + ' ' + sonde_data['TIME']
        sonde_data["Timestamp"] = pd.to_datetime(sonde_data['Timestamp'])
        sonde_data.drop(columns = ['DATE','TIME','Depth_m','Int_Batt_V' ], inplace= True) # inplace makes the drop action within the same target dataframe
        
        merged_bytime = pd.merge_asof(data, sonde_data, left_on='timestamp',right_on='Timestamp', direction='nearest')
        merged_bytime = merged_bytime.drop(columns = 'Timestamp')
    
    print("\nSonde data successfully merged")
    
    return merged_bytime

# =============================================================================
# Function to read, modify (if needed) and match ctd data with mims time csv.
# Raw file was premodified to this script to eliminate first rows of doc. data
# This function also calculates the depth. 
# It first averages pressure measurement of atmpospheric pressure and then 
# substracts that value from "submerged" values 
# to get the correct depth. Data range that represents atmopheric values must be defined

# Inputs:
#     - input_file: membrane stat file to attach all necessary ancilliary data for post-processing
#     - ctd_file: raw csv containing ctd data
#     - time_header: name of time column (case sensitive)
#     - atmos_p_range: list [a,b] of values to use for atmospheric average (based on index of the df)

# Output: 
#     - will merge depth (m) data to the membrane input file
# =============================================================================

def read_CTD(input_file, ctd_file, time_header, atmos_p_range):
    # time_header = 'Date/time'
    # atmos_p_range = [2484,2657]
    
    print("\nCompiling CTD data..")
    
    data = pd.read_csv(input_file)
    data['timestamp'] = pd.to_datetime(data['timestamp'])    
    
    ctd_data = read_csv_with_encoding_and_retry(ctd_file)
    ctd_data [time_header] = pd.to_datetime(ctd_data[time_header], dayfirst= True) # tells panda that the day (dd) is the first thing in the file format
        
    # Convert atmos_p_range to a slice object
    if isinstance(atmos_p_range, list):
        start_index = atmos_p_range[0]
        atmos_p_range = slice(*atmos_p_range) # * is used to "unpack'the elements within the variable

    elif not isinstance(atmos_p_range, slice):
        raise ValueError("atmos_p_range should be a list or slice object")
        
    # subset the data to take out only target portions of study interval
    ctd_data = ctd_data.iloc[start_index:]
    
    # subset the data again to take atmospheric pressure mean
    subset = ctd_data.loc[atmos_p_range]
    atmos_p_mean = subset['Pressure[cmH2O]'].mean()
    
    # calculate corrected depth 
    ctd_data['depth'] = (ctd_data['Pressure[cmH2O]'] - atmos_p_mean) / 100
        
    merged_bytime = pd.merge_asof(data, ctd_data[[time_header, 'depth']], left_on ='timestamp',right_on = time_header, 
                                  direction='nearest')
    merged_bytime = merged_bytime.drop(columns = time_header)
    
    print("\nCTD data successfully merged")

    return merged_bytime

# =============================================================================
# Function to read csv's. CTD data seems to be encoded using another codec that is not utf-8. 
# Here it prompts to read the csv with default parameters. If error it will then query what type 
# of encoding is present and output the string. If the codec matches "latin-1" (the one present for this
# type of file during testing) it will apply the read.csv using specific encoding 

# Inputs:
#     - input_file: membrane stat file to attach all necessary ancilliary data for post-processing

# Output: 
#     - will read the csv and store it in a variable
# =============================================================================   
    
def read_csv_with_encoding_and_retry(input_file):
    try:
        csv_data = pd.read_csv(input_file)
        
    except UnicodeDecodeError as pe:
        print(f"\nError reading the file {input_file} with encoding: {pe}") # With f-strings, you can embed expressions (such as variables or function calls) directly inside string literals, making string interpolation more concise and readable.
        
        with open(input_file, 'rb') as f:
            result = chardet.detect(f.read())
        print('\nFile is currently using encoding:', result['encoding'])
        
        if result['encoding'] == 'ISO-8859-1':
        
            csv_data = pd.read_csv(input_file, encoding= 'latin-1')
            
            print(f"\nSuccesful in reading file: {input_file} ")
            
        else:
        
            print("Other type of encoding is employed, switch manually in function file")
        
    return csv_data


# =============================================================================
# Testing section
# =============================================================================

# test = read_weather(file_2_merge , weather_file, time_header = 'Time')


# # # test2 = read_sonde(membrane_file, sonde_file)

# # # test3 = read_CTD(membrane_file, ctd_file, 'Date/time', [2484,2657])

# data = pd.read_csv(file_2_merge)
# data.sort_values(['timestamp', 'mz'], ascending = True, inplace = True)
# a = merge_ancilliary(file_2_merge,
#                           weather_file = weather_file, weather_header = 'Time')
#                            sonde_file = sonde_file,
#                            ctd_file = ctd_file, ctd_header = 'Date/time',
#                            ctd_range= [2484,2657])

