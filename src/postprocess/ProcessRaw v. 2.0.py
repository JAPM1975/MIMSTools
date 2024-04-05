# Functions to read in raw data and calculate statistics for each cycle and mass
# 
# ############################################################################


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# input_file = 'C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Membrane/Membrane_SC.csv'
# mz_list = ['N2' ]

# ##################################
# This function reads in a raw file and processes it (i.e. calculating mean and 
# std) and stores it in an output file.
# 
# Inputs:
#     - input_file: string of the file to read in and process
#     - mz_list: list of masses to evaluate e.g. mz_list = ['CH4', 'H2O', 'Ne 20']
#     - show_plots: Should plots be shown (very slow!) or not?
#     - refine: if the designated number of starting and last rows of data be eliminated

#
# Output: A dataframe with the processed data. The output also gets stored to a file
# #################################

def read_and_process_raw(input_file, mz_list, show_plots = False, refine = False):
  
  # Read the file
  data = pd.read_csv(input_file)
  
  # Convert to datetime to make sure data is in correct data type for time values
  data['TimeF'] = pd.to_datetime(data['TimeF'])
  data['TimeEM'] = pd.to_datetime(data['TimeEM'])
  
  # Separate into cycles because they do not get written by field code. We use floor division 
  # by the known ammount of measurements per cycle. Number (10) here 
  # can be changed depending on how many measurements the MIMS takes when changing
  # the time it takes to complete one valve cycle
  data['cycle'] = np.ceil(data.index // 10)
  
  # reshape the data and summarise by cycle and gas
  data_t = data.melt(id_vars = ['TimeF', 'TimeEM', 'cycle'], 
                     value_vars = mz_list, 
                     value_name = 'signal', 
                     var_name = 'mz')
  
  # Groups each row by cycle # and gas type into micro dataframes for later use
  data_f = data_t.groupby(by = ['cycle', 'mz']) 
  
  # Applies specific functions to micro dataframes by column. This calculates the 
  # mean and standard deviation (1 degree of freedom) for each gas grouped by cycle
  # and stores this mean and std value of each gas in a row per cycle.
  data_f = data_f.agg(signal_mean=pd.NamedAgg(column="signal", aggfunc="mean"),
                      signal_std=pd.NamedAgg(column="signal", aggfunc="std"),
                      TimeF=pd.NamedAgg(column="TimeF", aggfunc="mean"),
                      TimeEM=pd.NamedAgg(column="TimeEM", aggfunc="mean"))
  
  data_f = data_f.reset_index()
  
  # to calculate the standard error of the means (estimate population size based on sample)?
  data_f['signal_std'] = data_f['signal_std']/np.sqrt(10-1) # we want the std of the mean
  
  if refine == True:
   
      refined_data = data_f
      # The below min and max can be adjusted to take out certain cycles (rows of 
      # data) to adjust the overall fit
      
      first_cycle = refined_data['cycle'].min()
      last_cycle = refined_data['cycle'].max() - 2
      
      # This selects the data available between the designated delimiters
      refined_data = refined_data[(refined_data['cycle'] > first_cycle) & (refined_data['cycle'] < last_cycle)]

      refined_data = refined_data.reset_index()
      refined_data = refined_data.drop(columns = ['index'])
      
      refined_data.to_csv(input_file.replace('.csv','_stat.csv'))
      
      data_f = refined_data
  else:
      
      data_f.to_csv(input_file.replace(".csv", "_stat.csv"))  
  
  # Creates a unique entry for each cycle to be used as a counter for the loop.
  cycles = np.unique(data_t['cycle'])
    
  # Plot each cycle/gas data if necessary
  for c in cycles:
    for mz in mz_list:
      
      if (show_plots):
        plot_raw(data_t[(data_t['cycle'] == c) & (data_t['mz'] == mz)],
                 data_f[(data_f['cycle'] == c) & (data_f['mz'] == mz)])

  

  return(data_f)




# #####################################
# Plot the raw data and the calculated mean and std for a single cycle/gas
#
# Inputs:
#     - data: the raw data points from the MS for one cycle/gas
#     - stats: the summary data for the raw data
#
# Outputs:
#     - will show a plot
# #######################################

def plot_raw(data, stats):
  plt.figure(figsize=(10,6))
  plt.scatter(data['TimeF'], data['signal'], label='data')
  plt.scatter(stats['TimeF'], stats['signal_mean'], label='mean')
  plt.plot(data['TimeF'], [stats['signal_mean']] * len(data), label='mean')
  plt.plot(data['TimeF'], [stats['signal_mean'] + stats['signal_std']] * len(data), label='std - upper')
  plt.plot(data['TimeF'], [stats['signal_mean'] - stats['signal_std']] * len(data), label='std - lower')
  plt.xlabel('Time')
  plt.ylabel('Signal')
  plt.legend()
  plt.title(label= str(data['cycle']) + str(data['mz']))
  plt.show()
  
  
# =============================================================================
# Testing section
# =============================================================================

# Raw = pd.read_csv('C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Air/Air_SC.csv')

# Raw['TimeF'] = pd.to_datetime(Raw['TimeF'])
# Raw['TimeEM'] = pd.to_datetime(Raw['TimeEM'])

# Raw['cycle']= np.ceil(Raw.index // 10)

# mz_list = ['CH4',  'N2']

# Raw_2 = Raw.melt(id_vars = ['TimeF', 'TimeEM', 'cycle'], 
#                    value_vars = mz_list, 
#                    value_name = 'signal', 
#                    var_name = 'mz')

# Raw_3 = Raw_2.groupby(by = ['cycle','mz'])
# Raw_3 = Raw_3.agg(signal_mean=pd.NamedAgg(column="signal", aggfunc="mean"),
#                     signal_std=pd.NamedAgg(column="signal", aggfunc="std"),
#                     TimeF=pd.NamedAgg(column="TimeF", aggfunc="mean"),
#                     TimeEM=pd.NamedAgg(column="TimeEM", aggfunc="mean"))

# Raw_3 = Raw_3.reset_index()

# Raw_3['signal_std'] = Raw_3['signal_std']/np.sqrt(10-1)

# cycles = np.unique(Raw_3['cycle'])

# for c in cycles:
#   for mz in mz_list:
    
#       plot_raw(Raw_2[(Raw_2['cycle'] == c) & (Raw_2['mz'] == mz)])
#                # data_f[(data_f['cycle'] == c) & (data_f['mz'] == mz)])

# plt.figure(figsize=(10,6))
# plt.scatter(Raw_2[(Raw_2['TimeF'] == 1)], Raw_2[(Raw_2['signal'] ==1)], label='data')
# plt.scatter(Raw_3['TimeF'], Raw_3['signal_mean'], label='mean')
# plt.plot(data['TimeF'], [stats['signal_mean']] * len(data), label='mean')
# plt.plot(data['TimeF'], [stats['signal_mean'] + stats['signal_std']] * len(data), label='std - upper')
# plt.plot(data['TimeF'], [stats['signal_mean'] - stats['signal_std']] * len(data), label='std - lower')
# plt.xlabel('Time')
# plt.ylabel('Signal')
# plt.legend()
# plt.show()
