# This file provides functions to calibrate membrane measurements using blanks
# and standards
#
# 


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline






def calibrate(blank_file, standard_file, membrane_file, mz_list, mz_upper, mz_lower, composition_dry_air, atmospheric_pressure):
  blank_data = pd.read_csv(blank_file)
  blank_data['TimeF'] = pd.to_datetime(blank_data['TimeF'])
  standard_data = pd.read_csv(standard_file)
  standard_data['TimeF'] = pd.to_datetime(standard_data['TimeF'])
  membrane_data = pd.read_csv(membrane_file)
  membrane_data['TimeF'] = pd.to_datetime(membrane_data['TimeF'])

  
  # Fit each gas
  result = pd.DataFrame()
  for mz in mz_list:
    blank_data_mz = blank_data[blank_data['mz'] == mz]
    standard_data_mz = standard_data[standard_data['mz'] == mz]
    membrane_data_mz = membrane_data[membrane_data['mz'] == mz]
    p_mz = composition_dry_air[mz]*atmospheric_pressure
    
    result_temp = fit_signal_vs_time(blank_data_mz, standard_data_mz, membrane_data_mz, mz, p_mz)
    result = pd.concat([result, result_temp])
    
  # Fit some ratios
  result_temp = fit_ratios(result, mz_upper, mz_lower, composition_dry_air)
  
  result = pd.concat([result, result_temp])
  
  result.to_csv(membrane_file.replace('stat.csv', 'calibrated.csv'))
    





def fit_signal_vs_time(blank_data, standard_data, membrane_data, mz, p_mz, p_mz_unc=None):
    
    print(f"Postprocessing {mz}...")
    
    m = len(blank_data)
    smoothing_factor = np.mean((m-np.sqrt(2*m),m+np.sqrt(2*m)))
    # smoothing_factor = 150
    if p_mz_unc is None:
      p_mz_unc = p_mz * 0.01
    
    # fit blank
    time_numeric = pd.to_numeric(blank_data['TimeF'])   
    
    cs_blk = UnivariateSpline(x = time_numeric, y = blank_data['signal_mean'], w = 1/blank_data['signal_std'])
    cs_blk.set_smoothing_factor(smoothing_factor)
    residual = cs_blk(time_numeric) - blank_data['signal_mean']
    cs_blk_res = UnivariateSpline(x = time_numeric, y = abs(residual), w = 1/blank_data['signal_std'])
    cs_blk_res.set_smoothing_factor(smoothing_factor)
                        
    
    plot_signal(title = f"{mz} - blank", 
                time = blank_data['TimeF'], 
                signal_mean = blank_data['signal_mean'], 
                signal_unc = blank_data['signal_std'], 
                signal_fit = cs_blk(time_numeric),
                signal_fit_unc = cs_blk_res(time_numeric))
                
    
    
     # fit standard
    time_numeric = pd.to_numeric(standard_data['TimeF'])  
    
    cs_std = UnivariateSpline(x = time_numeric, y = standard_data['signal_mean'], w = 1/standard_data['signal_std'])
    cs_std.set_smoothing_factor(smoothing_factor)
    residual = cs_std(time_numeric) - standard_data['signal_mean']
    cs_std_res = UnivariateSpline(x = time_numeric, y = abs(residual), w = 1/standard_data['signal_std'])
    cs_std_res.set_smoothing_factor(smoothing_factor)
                          
    plot_signal(title = f"{mz} - standard", 
                time = standard_data['TimeF'], 
                signal_mean = standard_data['signal_mean'], 
                signal_unc = standard_data['signal_std'], 
                signal_fit = cs_std(time_numeric), 
                signal_fit_unc = cs_std_res(time_numeric),
                other_fit = cs_blk(time_numeric),
                other_label = "blank")


    # fit blk corrected standard
    std_mean_blk_corr = standard_data['signal_mean'] - cs_blk(time_numeric)
    std_mean_blk_corr_unc = np.sqrt(standard_data['signal_std']**2 + cs_blk_res(time_numeric)**2)
    
    cs_std2 = UnivariateSpline(x = time_numeric, y = std_mean_blk_corr, w = 1/std_mean_blk_corr_unc)
    cs_std2.set_smoothing_factor(smoothing_factor)
    residual = cs_std2(time_numeric) - std_mean_blk_corr
    cs_std2_res = UnivariateSpline(x = time_numeric, y = abs(residual), w = 1/std_mean_blk_corr_unc)
    cs_std2_res.set_smoothing_factor(smoothing_factor)
    
    plot_signal(title = f"{mz} - standard blk corr", 
                time = standard_data['TimeF'], 
                signal_mean = std_mean_blk_corr, 
                signal_unc = std_mean_blk_corr_unc, 
                signal_fit = cs_std2(time_numeric), 
                signal_fit_unc = cs_std2_res(time_numeric))
    
    
    # fit efficiency
    efficiency = std_mean_blk_corr/p_mz
    efficiency_unc = efficiency * np.sqrt((std_mean_blk_corr_unc/std_mean_blk_corr)**2 + (p_mz_unc/p_mz)**2)
    
    cs_eff = UnivariateSpline(x = time_numeric, y = efficiency, w = 1/efficiency_unc)
    cs_eff.set_smoothing_factor(smoothing_factor)
    residual = cs_eff(time_numeric) - efficiency
    cs_eff_res = UnivariateSpline(x = time_numeric, y = abs(residual), w = 1/efficiency_unc)
    cs_eff_res.set_smoothing_factor(smoothing_factor)
    
    plot_signal(title = f"{mz} - efficiency", 
                time = standard_data['TimeF'], 
                signal_mean = efficiency, 
                signal_unc = efficiency_unc, 
                signal_fit = cs_eff(time_numeric), 
                signal_fit_unc = cs_eff_res(time_numeric))
                # other_fit = cs_blk(time_numeric),
                # other_label = "blank")


    # fit sample
    time_numeric = pd.to_numeric(membrane_data['TimeF']) 
    
    cs_memb = UnivariateSpline(x = time_numeric, y = membrane_data['signal_mean'], w = 1/membrane_data['signal_std'])
    cs_memb.set_smoothing_factor(smoothing_factor)
    residual = cs_memb(time_numeric) - membrane_data['signal_mean']
    cs_memb_res = UnivariateSpline(x = time_numeric, y = abs(residual), w = 1/membrane_data['signal_std'])
    cs_memb_res.set_smoothing_factor(smoothing_factor)

    plot_signal(title = f"{mz} - membrane", 
                time = membrane_data['TimeF'], 
                signal_mean = membrane_data['signal_mean'], 
                signal_unc = membrane_data['signal_std'],
                signal_fit = cs_memb(time_numeric), 
                signal_fit_unc = cs_memb_res(time_numeric),
                other_fit = cs_blk(time_numeric),
                other_label = "blank")


    # fit blk corrected sample
    signal_mean_blk_corr = membrane_data['signal_mean'] - cs_blk(time_numeric)
    signal_mean_blk_corr_unc = np.sqrt(membrane_data['signal_std']**2 + cs_blk_res(time_numeric)**2)
    
    cs_memb2 = UnivariateSpline(x = time_numeric, y = signal_mean_blk_corr, w = 1/signal_mean_blk_corr_unc)
    cs_memb2.set_smoothing_factor(smoothing_factor)
    residual = cs_memb2(time_numeric) - signal_mean_blk_corr
    cs_memb2_res = UnivariateSpline(x = time_numeric, y = abs(residual), w = 1/signal_mean_blk_corr_unc)
    cs_memb2_res.set_smoothing_factor(smoothing_factor)

    plot_signal(title = f"{mz} - membrane blk corr", 
                time = membrane_data['TimeF'], 
                signal_mean = signal_mean_blk_corr, 
                signal_unc = signal_mean_blk_corr_unc,
                signal_fit = cs_memb2(time_numeric),
                signal_fit_unc = cs_memb2_res(time_numeric),
                other_fit = cs_std2(time_numeric),
                other_label = "standard")
     
     # fit efficieny corrected sample => gives partial pressures
    signal_mean_std_corr = signal_mean_blk_corr/cs_eff(time_numeric)
    signal_mean_std_corr_unc_rel = np.sqrt((signal_mean_blk_corr_unc/signal_mean_blk_corr)**2 + (cs_eff_res(time_numeric)/cs_eff(time_numeric))**2)
    signal_mean_std_corr_unc = signal_mean_std_corr * signal_mean_std_corr_unc_rel
    
    cs_memb3 = UnivariateSpline(x = time_numeric, y = signal_mean_std_corr, w = 1/signal_mean_std_corr_unc)
    residual = cs_memb3(time_numeric) - signal_mean_std_corr
    cs_memb3_res = UnivariateSpline(x = time_numeric, y = abs(residual), w = 1/signal_mean_std_corr_unc)
    cs_memb3_res.set_smoothing_factor(smoothing_factor)
    
    
    plot_signal(title = f"{mz} - membrane eff corr", 
                time = membrane_data['TimeF'], 
                signal_mean = signal_mean_std_corr, 
                signal_unc = signal_mean_std_corr_unc, 
                signal_fit = cs_memb3(time_numeric), 
                signal_fit_unc = cs_memb3_res(time_numeric),
                # other_fit = cs_blk(time_numeric),
                other_label = "blank")
    
    
    result = pd.DataFrame(list(zip(membrane_data['TimeF'], signal_mean_std_corr, signal_mean_std_corr_unc)), 
                          columns=['timestamp', 'signal_calibrated', 'signal_calibrated_unc'])

    result['mz'] = mz
    
    return(result)



  
  
  
      

def fit_ratios(result, upper, lower, composition_dry_air):
  
  print(f"Postprocessing ratios...")
  
  n = len(upper)
  m = len(result)
  smoothing_factor = np.mean((m-np.sqrt(2*m),m+np.sqrt(2*m)))
  # smoothing_factor = 150

  results = pd.DataFrame()
  for i in range(0,n):
    calib_data_u = result[result['mz'] == upper[i]]
    calib_data_l = result[result['mz'] == lower[i]]
    
    # Calculate ratio
    
    time_numeric = pd.to_numeric(calib_data_u['timestamp']) 
    ratio = calib_data_u['signal_calibrated'] / calib_data_l['signal_calibrated']
    ratio_unc = ratio * np.sqrt((calib_data_u['signal_calibrated_unc']/calib_data_u['signal_calibrated'])**2 + (calib_data_l['signal_calibrated_unc']/calib_data_l['signal_calibrated'])**2)
    
    
    cs_ratio = UnivariateSpline(x = time_numeric, y = ratio, w = 1/ratio_unc, s = smoothing_factor)
    residual = cs_ratio(time_numeric) - ratio
    cs_ratio_unc = UnivariateSpline(x = time_numeric, y = abs(residual), w = 1/ratio_unc, s = smoothing_factor)
    
    plot_signal(title = f"{upper[i]}/{lower[i]} - sample ratio",
                time = calib_data_u['timestamp'],
                signal_mean = ratio,
                signal_unc = ratio_unc,
                signal_fit = cs_ratio(time_numeric),
                signal_fit_unc = cs_ratio_unc(time_numeric),
                other_fit = [composition_dry_air[upper[i]]/composition_dry_air[lower[i]]] * len(time_numeric),
                other_label = "air")
                # ylim = [composition_dry_air[upper[i]]/composition_dry_air[lower[i]]/2, composition_dry_air[upper[i]]/composition_dry_air[lower[i]]*2])

    result_temp = pd.DataFrame(list(zip(calib_data_u['timestamp'], ratio, ratio_unc)), 
                               columns=['timestamp', 'signal_calibrated', 'signal_calibrated_unc'])

    result_temp['mz'] = f"{upper[i]}/{lower[i]}"
    results = pd.concat([results, result_temp])
    
  return(results)






def plot_signal(title, time, signal_mean, signal_unc=None, signal_fit=None, signal_fit_unc=None, other_fit=None, other_label=None, ylim=None):
  
  # calculate y limits
  min_y = min(signal_mean)
  if signal_fit is not None:
    min_y = min([min_y, min(signal_fit)])
  if signal_fit_unc is not None:
    min_y = min([min_y, min(signal_fit - signal_fit_unc)])
  if other_fit is not None:
    min_y =min([min_y, min(other_fit)])
    
  max_y = max(signal_mean)
  if signal_fit is not None:
    max_y = max([max_y, max(signal_fit)])
  if signal_fit_unc is not None:
    max_y = max([max_y, max(signal_fit + signal_fit_unc)])
  if other_fit is not None:
    max_y = max([max_y, max(other_fit)])
  
  min_y = min_y - 0.05 * (max_y - min_y)
  max_y = max_y + 0.05 * (max_y - min_y)
  
  plt.figure(figsize=(10,6))
  plt.scatter(time, signal_mean, label='data', zorder = 10, s = 10)
  # min_y = min(signal_mean)
  # max_y = max(signal_mean)
  if signal_unc is not None:
    plt.errorbar(x = time, y = signal_mean, yerr = signal_unc, 
                 color = 'grey', elinewidth = 0.5, linewidth=0)
  if signal_fit is not None:
    # min_y = min([min_y, min(signal_fit)])
    # max_y = max([max_y, max(signal_fit)])
    plt.plot(time, signal_fit, label = "fit", color = "darkred", zorder = 10)
  if signal_fit_unc is not None:
    plt.fill_between(time, signal_fit + signal_fit_unc, signal_fit - signal_fit_unc, label = "uncertainty", color = "darkred", alpha = 0.3)
  if other_fit is not None:
    plt.plot(time, other_fit, label = other_label, color = "orange")
  plt.xlabel('Time')
  plt.ylabel('Signal')
  plt.legend()
  plt.title(title)
  if (not np.isnan(min_y)) and (not np.isinf(min_y)) and (not np.isnan(max_y)) and (not np.isinf(max_y)):
    plt.ylim([min_y, max_y])
  if ylim is not None:
    plt.ylim(ylim)
  # plt.xticks()
  # plt.Axes.xaxis_date()
  # plt.Axes.axes.s
  # plt.Axes.set_ymargin()
  # plt.Axes.
  plt.show()

