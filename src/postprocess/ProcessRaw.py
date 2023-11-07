# Functions to read in raw data and calculate statistics for each cycle and mass
# 
# ############################################################################


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ##################################
# This function reads in a raw file and processes it (i.e. calculating mean and 
# std) and stores it in an output file.
# 
# Inputs:
#     - input_file: string of the file to read in and process
#     - mz_list: list of masses to evaluate e.g. mz_list = ['CH4', 'H2O', 'Ne 20']
#     - show_plots: Should plots be shown (very slow!) or not?
#
# Output: A dataframe with the processed data. The output also gets stored to a file
# #################################

def read_and_process_raw(input_file, mz_list, show_plots = False):
  
  # Read the file
  data = pd.read_csv(input_file)
  data['TimeF'] = pd.to_datetime(data['TimeF'])
  data['TimeEM'] = pd.to_datetime(data['TimeEM'])
  
  # Identify cycles because they do not get written by field code
  data['cycle'] = np.ceil(data.index // 10)
  
  # reshape the data and summarise by cycle and gas
  data_t = data.melt(id_vars = ['TimeF', 'TimeEM', 'cycle'], 
                     value_vars = mz_list, 
                     value_name = 'signal', 
                     var_name = 'mz')
  data_f = data_t.groupby(by = ['cycle', 'mz']) 
  data_f = data_f.agg(signal_mean=pd.NamedAgg(column="signal", aggfunc="mean"),
                      signal_std=pd.NamedAgg(column="signal", aggfunc="std"),
                      TimeF=pd.NamedAgg(column="TimeF", aggfunc="mean"),
                      TimeEM=pd.NamedAgg(column="TimeEM", aggfunc="mean"))
  data_f = data_f.reset_index()
  data_f['signal_std'] = data_f['signal_std']/np.sqrt(10-1) # we want the std of the mean
  
  
  # Plot each cycle/gas data if necessary
  cycles = np.unique(data_t['cycle'])
  for c in cycles:
    for mz in mz_list:
      
      if (show_plots):
        plot_raw(data_t[(data_t['cycle'] == c) & (data_t['mz'] == mz)],
                 data_f[(data_f['cycle'] == c) & (data_f['mz'] == mz)])
  
  # Save the data to the outputfile
  data_f.to_csv(input_file.replace(".csv", "_stat.csv"))    

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
  plt.show()
  
