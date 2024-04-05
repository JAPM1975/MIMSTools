



import os
import pandas as pd
import datetime
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import sklearn
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score


from scipy.interpolate import UnivariateSpline


from src.postprocess.ProcessRaw import read_and_process_raw
from src.postprocess.Ancilliary_data_compiler import merge_ancilliary
from src.postprocess.Calibrate import calibrate


# =============================================================================
# Import Ancialliary
# =============================================================================

weather_file = "C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/South Creek/MIMS/Surface Water/Atmosphere Data/weather_filtered.csv"
sonde_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Manta/Manta/HAWKS_4_SC_modified.csv'
ctd_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/Raw Logger data/CTD_Divers/House_A_Rickabys_SthCreek/CSV/Atmoshperic pressure_231113152722_V9769 - Copy.csv'

# Process the raw files
# #######################

# 48H time series

mz_list = ['CH4', 'H2O', 'Ne 20', 'N2', 'O2', 'Ar', 'CO2', 'He', 'Ne 22', 'Kr', 'Xe']

composition_dry_air = {'N2': 0.780840, 
                       'O2': 0.209460,
                       'Ar': 0.009340,
                       'CO2': 0.000417,
                       'CH4': 0.00000187,
                       'Ne': 0.00001818,
                       'He': 0.00000524,
                       'Kr': 0.00000114,
                       'Xe': 0.000000087,
                       'Ne 20': 0.90483*0.00001818,
                       'Ne 22': 0.09253*0.00001818,
                       'H2O': 1/761}

atmospheric_pressure = 761 # mmHg
# atmospheric_pressure = 1013.25 # hPa


mz_upper = ['Ar', 'Ar']
mz_lower = ['Ar', 'Ar']
mz_upper = ['N2', 'O2']

# file_air = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Air/Air.csv'
# result_air = read_and_process_raw(file_air, mz_list)

# file_blank = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Blanks/Blank.csv'
# result_blank = read_and_process_raw(file_blank, mz_list)

# file_membrane = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Membrane/Membrane.csv'
# result_membrane = read_and_process_raw(file_membrane, mz_list)


# # Deep GW data

# file_air = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Deep GW/Deep GW Maroota/Valves/Air/AirAir.csv'
# result_air = read_and_process_raw(file_air, mz_list)

# file_blank = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Deep GW/Deep GW Maroota/Valves/Blanks/Blank.csv'
# result_blank = read_and_process_raw(file_blank, mz_list)

# file_membrane = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Deep GW/Deep GW Maroota/Valves/Membrane/Membrane.csv'
# result_membrane = read_and_process_raw(file_membrane, mz_list)

# =============================================================================
# # # Second SW time series
# =============================================================================

# file_air = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Air/Air_Nov_4.csv'
# result_air = read_and_process_raw(file_air, mz_list)

# file_blank = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Blanks/Blank_Nov_4.csv'
# result_blank = read_and_process_raw(file_blank, mz_list)

# file_membrane = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Membrane/Membrane_Nov_4.csv'
# result_membrane = read_and_process_raw(file_membrane, mz_list)

# =============================================================================
# Survey
# =============================================================================

# file_air = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Survey/Valves/Air/Air3.csv'
# result_air = read_and_process_raw(file_air, mz_list)

# file_blank = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Survey/Valves/Blanks/Blank3.csv'
# result_blank = read_and_process_raw(file_blank, mz_list)

# file_membrane = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Survey/Valves/Membrane/Membrane3.csv'
# result_membrane = read_and_process_raw(file_membrane, mz_list)

# # =============================================================================
# South Creek
# =============================================================================

file_air = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Air/Air_SC.csv'
result_air = read_and_process_raw(file_air, mz_list, refine= True)

file_blank = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Blanks/Blank_SC.csv'
result_blank = read_and_process_raw(file_blank, mz_list, refine= True)

file_membrane = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Membrane/Membrane_SC.csv'
result_membrane = read_and_process_raw(file_membrane, mz_list, refine= True)

# Do calibration
# #########################


blank_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Blanks/Blank_SC_stat.csv'
# refined_blank = refine_data(blank_file, mz_list)
# R_blank_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Blanks/Blank_SC_refined_stat.csv'


standard_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Air/Air_SC_stat.csv'
# refined_standard = refine_data(standard_file, mz_list)
# R_standard_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Air/Air_SC_refined_stat.csv'

membrane_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Membrane/Membrane_SC_stat.csv'



# refined_membrane = refine_data(membrane_file, mz_list)
# R_membrane_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Membrane/Membrane_SC_refined_stat.csv'

result_calibration = calibrate(blank_file, standard_file, membrane_file , mz_list, mz_upper, mz_lower, composition_dry_air, weather_file)

# calibrate(R_blank_file, R_standard_file, R_membrane_file, mz_list, mz_upper, mz_lower, composition_dry_air, atmospheric_pressure)
# membrane_ancilliary = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Membrane/Membrane_SC_stat_ancilliary.csv'
calibrated_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Membrane/Membrane_SC_calibrated.csv'

a = merge_ancilliary(calibrated_file,
                          weather_file = weather_file, weather_header = 'Time',
                          sonde_file = sonde_file,
                          ctd_file = ctd_file, ctd_header = 'Date/time',
                          ctd_range= [2484,2657])

# =============================================================================
# First House Time series
# =============================================================================

# blank_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Air/Air_stat.csv'
# standard_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Blanks/Blank_stat.csv'
# membrane_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Membrane/Membrane_stat.csv'
# calibrate(blank_file, standard_file, membrane_file, mz_list, mz_upper, mz_lower, composition_dry_air, atmospheric_pressure)

# # =============================================================================
# Deep GW Maroota
# =============================================================================

# blank_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Deep GW/Deep GW Maroota/Valves/Blanks/Blank_stat.csv'
# standard_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Deep GW/Deep GW Maroota/Valves/Air/AirAir_stat.csv'
# membrane_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Deep GW/Deep GW Maroota/Valves/Membrane/Membrane_stat.csv'
# calibrate(blank_file, standard_file, membrane_file, mz_list, mz_upper, mz_lower, composition_dry_air, atmospheric_pressure)

# =============================================================================
# Second House Time series
# =============================================================================
# blank_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Blanks/Blank_Nov_4_stat.csv'
# standard_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Air/Air_Nov_4_stat.csv'
# membrane_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Membrane/Membrane_Nov_4_stat.csv'
# calibrate(blank_file, standard_file, membrane_file, mz_list, mz_upper, mz_lower, composition_dry_air, atmospheric_pressure)


# =============================================================================
# Survey
# =============================================================================
# blank_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Survey/Valves/Blanks/Blank3_stat.csv'
# standard_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Survey/Valves/Air/Air3_stat.csv'
# membrane_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Survey/Valves/Membrane/Membrane3_stat.csv'
# calibrate(blank_file, standard_file, membrane_file, mz_list, mz_upper, mz_lower, composition_dry_air, atmospheric_pressure)




