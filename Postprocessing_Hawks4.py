



# import os
# import pandas as pd
# import datetime
# import seaborn as sns
# import matplotlib.pyplot as plt
# import matplotlib.dates as mdates
# import numpy as np
# import sklearn
# from sklearn.linear_model import LinearRegression
# from sklearn.metrics import r2_score


# from scipy.interpolate import UnivariateSpline



from src.postprocess.ProcessRaw import read_and_process_raw
from src.postprocess.Calibrate import calibrate


# Process the raw files
# #######################

# 48H time series

mz_list = ['CH4', 'H2O', 'Ne 20', 'N2', 'O2', 'Ar', 'CO2', 'He', 'Ne 22', 'Kr', 'Xe']

file_air = 'Data/Hawks4/SW1/Valves/Air/Air.csv'
result_air = read_and_process_raw(file_air, mz_list)

file_blank = 'Data/Hawks4/SW1/Valves/Blanks/Blank.csv'
result_blank = read_and_process_raw(file_blank, mz_list)

file_membrane = 'Data/Hawks4/SW1/Valves/Membrane/Membrane.csv'
result_membrane = read_and_process_raw(file_membrane, mz_list)


# Deep GW data

file_air = 'Data/Hawks4/Deep GW/Valves/Air/Air.csv'
result_air = read_and_process_raw(file_air, mz_list)

file_blank = 'Data/Hawks4/Deep GW/Valves/Blanks/Blank.csv'
result_blank = read_and_process_raw(file_blank, mz_list)

file_membrane = 'Data/Hawks4/Deep GW/Valves/Membrane/Membrane.csv'
result_membrane = read_and_process_raw(file_membrane, mz_list)

# Second SW time series

file_air = 'Data/Hawks4/SW2/Valves/Air/Air_Nov_4.csv'
result_air = read_and_process_raw(file_air, mz_list)

file_blank = 'Data/Hawks4/SW2/Valves/Blanks/Blank_Nov_4.csv'
result_blank = read_and_process_raw(file_blank, mz_list)

file_membrane = 'Data/Hawks4/SW2/Valves/Membrane/Membrane_Nov_4.csv'
result_membrane = read_and_process_raw(file_membrane, mz_list)

# Do calibration
# #########################

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


mz_list = ['N2', 'O2', 'Ar', 'CH4', 'CO2', 'Kr', 'He']
mz_upper = ['Ar', 'Ar']
mz_lower = ['N2', 'O2']

blank_file = 'Data/Hawks4/SW1/Valves/Blanks/Blank_stat.csv'
standard_file = 'Data/Hawks4/SW1/Valves/Air/Air_stat.csv'
membrane_file = 'Data/Hawks4/SW1/Valves/Membrane/Membrane_stat.csv'
calibrate(blank_file, standard_file, membrane_file, mz_list, mz_upper, mz_lower, composition_dry_air, atmospheric_pressure)

blank_file = 'Data/Hawks4/Deep GW/Valves/Blanks/Blank_stat.csv'
standard_file = 'Data/Hawks4/Deep GW/Valves/Air/Air_stat.csv'
membrane_file = 'Data/Hawks4/Deep GW/Valves/Membrane/Membrane_stat.csv'
calibrate(blank_file, standard_file, membrane_file, mz_list, mz_upper, mz_lower, composition_dry_air, atmospheric_pressure)

blank_file = 'Data/Hawks4/SW2/Valves/Blanks/Blank_Nov_4_stat.csv'
standard_file = 'Data/Hawks4/SW2/Valves/Air/Air_Nov_4_stat.csv'
membrane_file = 'Data/Hawks4/SW2/Valves/Membrane/Membrane_Nov_4_stat.csv'
calibrate(blank_file, standard_file, membrane_file, mz_list, mz_upper, mz_lower, composition_dry_air, atmospheric_pressure)




