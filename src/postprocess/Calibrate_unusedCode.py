






# 
# 
# def fit_ratio_vs_time(blank_data_1, standard_data_1, membrane_data_1, blank_data_2, standard_data_2, membrane_data_2, mz1, mz2):
#     
#     m = len(blank_data)
#     smoothing_factor = np.mean((m-np.sqrt(2*m),m+np.sqrt(2*m)))
#     # smoothing_factor = 150
#     
#     # fit blank 1
#     time_numeric = pd.to_numeric(blank_data1['TimeF'])   
#     
#     cs_blk1 = UnivariateSpline(x = time_numeric, y = blank_data1['signal_mean'], w = 1/blank_data1['signal_std'], s = smoothing_factor)
#     residual = cs_blk1(time_numeric) - blank_data1['signal_mean']
#     cs_blk_res1 = UnivariateSpline(x = time_numeric, y = abs(residual), w = 1/blank_data1['signal_std'], s = smoothing_factor)
#                         
#     
#     plot_signal(title = f"{mz1} - blank", 
#                 time = blank_data1['TimeF'], 
#                 signal_mean = blank_data1['signal_mean'], 
#                 signal_unc = blank_data1['signal_std'], 
#                 signal_fit = cs_blk1(time_numeric),
#                 signal_fit_unc = cs_blk_res1(time_numeric))
#                 
#     # fit blank 2
#     time_numeric = pd.to_numeric(blank_data2['TimeF'])   
#     
#     cs_blk2 = UnivariateSpline(x = time_numeric, y = blank_data2['signal_mean'], w = 1/blank_data2['signal_std'], s = smoothing_factor)
#     residual = cs_blk1(time_numeric) - blank_data2['signal_mean']
#     cs_blk_res2 = UnivariateSpline(x = time_numeric, y = abs(residual), w = 1/blank_data2['signal_std'], s = smoothing_factor)
#                         
#     
#     plot_signal(title = f"{mz2} - blank", 
#                 time = blank_data2['TimeF'], 
#                 signal_mean = blank_data2['signal_mean'], 
#                 signal_unc = blank_data2['signal_std'], 
#                 signal_fit = cs_blk2(time_numeric),
#                 signal_fit_unc = cs_blk_res2(time_numeric))
#                 
#     
#     
#      
#     # fit blk corrected standard 1
#     time_numeric1 = pd.to_numeric(standard_data1['TimeF']) 
#     std_mean_blk_corr1 = standard_data1['signal_mean'] - cs_blk1(time_numeric1)
#     std_mean_blk_corr_unc1 = np.sqrt(standard_data1['signal_std']**2 + cs_blk_res1(time_numeric1)**2)
#     
#     cs_std1 = UnivariateSpline(x = time_numeric1, y = std_mean_blk_corr1, w = 1/std_mean_blk_corr_unc1, s = smoothing_factor)
#     residual = cs_std1(time_numeric1) - std_mean_blk_corr1
#     cs_std_res1 = UnivariateSpline(x = time_numeric1, y = abs(residual), w = 1/std_mean_blk_corr_unc1, s = smoothing_factor)
#     
#     plot_signal(title = f"{mz1} - standard blk corr", 
#                 time = standard_data1['TimeF'], 
#                 signal_mean = std_mean_blk_corr1, 
#                 signal_unc = std_mean_blk_corr_unc1, 
#                 signal_fit = cs_std1(time_numeric1), 
#                 signal_fit_unc = cs_std_res1(time_numeric1))
#                 # other_fit = cs_blk(time_numeric1),
#                 # other_label = "blank")
#                 
#     # fit blk corrected standard 2
#     time_numeric2 = pd.to_numeric(standard_data2['TimeF']) 
#     std_mean_blk_corr2 = standard_data2['signal_mean'] - cs_blk2(time_numeric2)
#     std_mean_blk_corr_unc2 = np.sqrt(standard_data2['signal_std']**2 + cs_blk_res2(time_numeric2)**2)
#     
#     cs_std2 = UnivariateSpline(x = time_numeric2, y = std_mean_blk_corr2, w = 1/std_mean_blk_corr_unc2, s = smoothing_factor)
#     residual = cs_std2(time_numeric2) - std_mean_blk_corr2
#     cs_std_res2 = UnivariateSpline(x = time_numeric2, y = abs(residual), w = 1/std_mean_blk_corr_unc2, s = smoothing_factor)
#     
#     plot_signal(title = f"{mz2} - standard blk corr", 
#                 time = standard_data2['TimeF'], 
#                 signal_mean = std_mean_blk_corr2, 
#                 signal_unc = std_mean_blk_corr_unc2, 
#                 signal_fit = cs_std2(time_numeric2), 
#                 signal_fit_unc = cs_std_res2(time_numeric2))
#                 # other_fit = cs_blk(time_numeric2),
#                 # other_label = "blank")
#                 
#     # Calculate ratio for standards
#     ratio1_std = std_mean_blk_corr1 / cs_std2(time_numeric1)
#     ratio2_std = cs_std1(time_numeric2) / std_mean_blk_corr2
#     ratio_std = pd.concat([ratio1_std, ratio2_std])
#     time_std = pd.concat([time_numeric1, time_numeric2])
#     ratio_std.sort_index(inplace=True)
#     time_std.sort_index(inplace=True)
#     ratio1_std_unc = ratio1_std * np.sqrt((std_mean_blk_corr_unc1/std_mean_blk_corr1)**2 + (cs_std_res2(time_numeric1)/cs_std2(time_numeric1))**2)
#     ratio2_std_unc = ratio2_std * np.sqrt((cs_std_res1(time_numeric2)/cs_std1(time_numeric2))**2 + (std_mean_blk_corr_unc2/std_mean_blk_corr2)**2)
#     ratio_std_unc = pd.concat([ratio1_std_unc, ratio2_std_unc])
#     ratio_std_unc.sort_index(inplace=True)
#     
#     cs_std_ratio = UnivariateSpline(x = time_std, y = ratio_std, w = 1/ratio_std_unc, s = smoothing_factor)
#     residual = cs_std_ratio(time_std) - ratio_std
#     cs_std_ratio_res = UnivariateSpline(x = time_std, y = abs(residual), w = 1/ratio_std_unc, s = smoothing_factor)
#     
#     plot_signal(title = f"{mz1}/{mz2} - standard ratio",
#                 time = time_std,
#                 signal_mean = ratio_std,
#                 signal_unc = ratio_std_unc,
#                 signal_fit = cs_std_ratio(time_std),
#                 signal_fit_unc = cs_std_ratio_res(time_std))
#     
# 
#     # fit blk corrected sample 1
#     time_numeric1 = pd.to_numeric(membrane_data1['TimeF']) 
#     
#     signal_mean_blk_corr1 = membrane_data1['signal_mean'] - cs_blk1(time_numeric1)
#     signal_mean_blk_corr1_unc = np.sqrt(membrane_data1['signal_std']**2 + cs_blk_res1(time_numeric1)**2)
#     
#     cs_memb1 = UnivariateSpline(x = time_numeric1, y = signal_mean_blk_corr1, s = smoothing_factor, w = 1/signal_mean_blk_corr1_unc)
#     residual = cs_memb1(time_numeric) - signal_mean_blk_corr1
#     cs_memb1_res = UnivariateSpline(x = time_numeric1, y = abs(residual), s = smoothing_factor, w = 1/signal_mean_blk_corr1_unc)
# 
#     plot_signal(title = f"{mz1} - membrane blk corr", 
#                 time = membrane_data1['TimeF'], 
#                 signal_mean = signal_mean_blk_corr1, 
#                 signal_unc = signal_mean_blk_corr1_unc,
#                 signal_fit = cs_memb1(time_numeric1),
#                 signal_fit_unc = cs_memb1_res(time_numeric),
#                 other_fit = cs_std1(time_numeric1),
#                 other_label = "standard")
#                 
#                 
#     # fit blk corrected sample 2
#     time_numeric2 = pd.to_numeric(membrane_data2['TimeF']) 
#     
#     signal_mean_blk_corr2 = membrane_data2['signal_mean'] - cs_blk2(time_numeric1)
#     signal_mean_blk_corr2_unc = np.sqrt(membrane_data2['signal_std']**2 + cs_blk_res2(time_numeric2)**2)
#     
#     cs_memb2 = UnivariateSpline(x = time_numeric2, y = signal_mean_blk_corr2, s = smoothing_factor, w = 1/signal_mean_blk_corr2_unc)
#     residual = cs_memb2(time_numeric) - signal_mean_blk_corr2
#     cs_memb2_res = UnivariateSpline(x = time_numeric2, y = abs(residual), s = smoothing_factor, w = 1/signal_mean_blk_corr2_unc)
# 
#     plot_signal(title = f"{mz2} - membrane blk corr", 
#                 time = membrane_data2['TimeF'], 
#                 signal_mean = signal_mean_blk_corr2, 
#                 signal_unc = signal_mean_blk_corr2_unc,
#                 signal_fit = cs_memb2(time_numeric2),
#                 signal_fit_unc = cs_memb2_res(time_numeric),
#                 other_fit = cs_std2(time_numeric2),
#                 other_label = "standard")
#                 
#                 
#     # Calculate ratio for samples
#     ratio1_memb = signal_mean_blk_corr1 / cs_memb2(time_numeric1)
#     ratio2_memb = cs_memb1(time_numeric2) / signal_mean_blk_corr2
#     ratio_memb = pd.concat([ratio1_memb, ratio2_memb])
#     time_memb = pd.concat([time_numeric1, time_numeric2])
#     ratio_memb.sort_index(inplace=True)
#     time_memb.sort_index(inplace=True)
#     ratio1_memb_unc = ratio1_memb * np.sqrt((signal_mean_blk_corr1_unc/signal_mean_blk_corr1)**2 + (cs_memb2_res(time_numeric1)/cs_memb2(time_numeric1))**2)
#     ratio2_memb_unc = ratio2_memb * np.sqrt((cs_memb1_res(time_numeric2)/cs_memb1(time_numeric2))**2 + (signal_mean_blk_corr2_unc/signal_mean_blk_corr2)**2)
#     ratio_memb_unc = pd.concat([ratio1_memb_unc, ratio2_memb_unc])
#     ratio_memb_unc.sort_index(inplace=True)
#     
#     cs_memb_ratio = UnivariateSpline(x = time_memb, y = ratio_memb, w = 1/ratio_memb_unc, s = smoothing_factor)
#     residual = cs_memb_ratio(time_memb) - ratio_memb
#     cs_memb_ratio_res = UnivariateSpline(x = time_memb, y = abs(residual), w = 1/ratio_memb_unc, s = smoothing_factor)
#     
#     plot_signal(title = f"{mz1}/{mz2} - sample ratio",
#                 time = time_memb,
#                 signal_mean = ratio_memb,
#                 signal_unc = ratio_memb_unc,
#                 signal_fit = cs_memb_ratio(time_memb),
#                 signal_fit_unc = cs_memb_ratio_res(time_memb))
#      
#     # fit std corrected sample 1
#     signal_mean_std_corr1 = signal_mean_blk_corr1/cs_std1(time_numeric1)
#     signal_mean_std_corr1_unc_rel = np.sqrt((signal_mean_blk_corr1_unc/signal_mean_blk_corr1)**2 + (cs_std_res1(time_numeric1)/cs_std1(time_numeric1))**2)
#     signal_mean_std_corr1_unc = signal_mean_std_corr1 * signal_mean_std_corr1_unc_rel
#     
#     cs_memb1b = UnivariateSpline(x = time_numeric1, y = signal_mean_std_corr1, w = 1/signal_mean_std_corr1_unc)
#     residual = cs_memb1b(time_numeric1) - signal_mean_std_corr1
#     cs_memb1b_res = UnivariateSpline(x = time_numeric1, y = abs(residual), w = 1/signal_mean_std_corr1_unc)
#     cs_memb1b_res.set_smoothing_factor(smoothing_factor)
#     
#     
#     plot_signal(title = f"{mz1} - membrane std corr", 
#                 time = membrane_data1['TimeF'], 
#                 signal_mean = signal_mean_std_corr1, 
#                 signal_unc = signal_mean_std_corr1_unc, 
#                 signal_fit = cs_memb1b(time_numeric1), 
#                 signal_fit_unc = cs_memb1b_res(time_numeric1),
#                 # other_fit = cs_blk(time_numeric),
#                 other_label = "blank")
#     
#     # fit std corrected sample 2
#     signal_mean_std_corr2 = signal_mean_blk_corr2/cs_std1(time_numeric2)
#     signal_mean_std_corr2_unc_rel = np.sqrt((signal_mean_blk_corr2_unc/signal_mean_blk_corr2)**2 + (cs_std_res2(time_numeric1)/cs_std2(time_numeric2))**2)
#     signal_mean_std_corr2_unc = signal_mean_std_corr2 * signal_mean_std_corr2_unc_rel
#     
#     cs_memb2b = UnivariateSpline(x = time_numeric2, y = signal_mean_std_corr2, w = 1/signal_mean_std_corr2_unc)
#     residual = cs_memb2b(time_numeric2) - signal_mean_std_corr2
#     cs_memb2b_res = UnivariateSpline(x = time_numeric2, y = abs(residual), w = 1/signal_mean_std_corr2_unc)
#     cs_memb2b_res.set_smoothing_factor(smoothing_factor)
#     
#     
#     plot_signal(title = f"{mz2} - membrane std corr", 
#                 time = membrane_data2['TimeF'], 
#                 signal_mean = signal_mean_std_corr2, 
#                 signal_unc = signal_mean_std_corr2_unc, 
#                 signal_fit = cs_memb2b(time_numeric2), 
#                 signal_fit_unc = cs_memb2b_res(time_numeric2),
#                 # other_fit = cs_blk(time_numeric),
#                 other_label = "blank")
#     
#     
#     # Calculate ratio for samples after std correction
#     ratio1_memb = signal_mean_std_corr1 / cs_memb2b(time_numeric1)
#     ratio2_memb = cs_memb1b(time_numeric2) / signal_mean_std_corr2
#     ratio_memb_std_corr = pd.concat([ratio1_memb, ratio2_memb])
#     time_memb_std_corr = pd.concat([time_numeric1, time_numeric2])
#     ratio_memb_std_corr.sort_index(inplace=True)
#     time_memb_std_corr.sort_index(inplace=True)
#     ratio1_memb_unc = ratio1_memb * np.sqrt((signal_mean_std_corr1_unc/signal_mean_std_corr1)**2 + (cs_memb2b_res(time_numeric1)/cs_memb2b(time_numeric1))**2)
#     ratio2_memb_unc = ratio2_memb * np.sqrt((cs_memb1b_res(time_numeric2)/cs_memb1b(time_numeric2))**2 + (signal_mean_std_corr2_unc/signal_mean_std_corr2)**2)
#     ratio_memb_std_corr_unc = pd.concat([ratio1_memb_unc, ratio2_memb_unc])
#     ratio_memb_std_corr_unc.sort_index(inplace=True)
#     
#     cs_memb_ratio_std_corr = UnivariateSpline(x = time_memb_std_corr, y = ratio_memb_std_corr, w = 1/ratio_memb_std_corr_unc, s = smoothing_factor)
#     residual = cs_memb_ratio_std_corr(time_memb_std_corr) - ratio_memb_std_corr
#     cs_memb_ratio_std_corr_res = UnivariateSpline(x = time_memb_std_corr, y = abs(residual), w = 1/ratio_memb_std_corr_unc, s = smoothing_factor)
#     
#     plot_signal(title = f"{mz1}/{mz2} - sample ratio",
#                 time = time_memb_std_corr,
#                 signal_mean = ratio_memb_std_corr,
#                 signal_unc = ratio_memb_std_corr_unc,
#                 signal_fit = cs_memb_ratio_std_corr(time_memb_std_corr),
#                 signal_fit_unc = cs_memb_ratio_std_corr_res(time_memb_std_corr),
#                 ylim = [0.002, 0.016])
#   



  
  # blank_fit = fit_signal_vs_time(blank_data, mz_list, type = 'blank')
  # 
  # standard_fit = fit_signal_vs_time(standard_data, mz_list, type = 'standard')
  


    
  # blank_data_Ar = blank_data[blank_data['mz'] == 'Ar']
  # standard_data_Ar = standard_data[standard_data['mz'] == 'Ar']
  # membrane_data_Ar = membrane_data[membrane_data['mz'] == 'Ar']
  # blank_data_N2 = blank_data[blank_data['mz'] == 'N2']
  # standard_data_N2 = standard_data[standard_data['mz'] == 'N2']
  # membrane_data_N2 = membrane_data[membrane_data['mz'] == 'N2']
  # 
  # fit_ratio_vs_time(blank_data_Ar, standard_data_Ar, membrane_data_Ar, 
  #                   blank_data_N2, standard_data_N2, membrane_data_N2,  'Ar', 'N2')
  



  
# def fit_signal_vs_time(data, mz_list, type):
#   
#   for mz in mz_list:
#     data_mz = data[data['mz'] == mz]
#     
#     
#     degree = 1
#     time_numeric = pd.to_numeric(data_mz['TimeF'])
#     poly_fit,poly_cov = np.polyfit(x = time_numeric, y = data_mz['signal_mean'], deg = degree, 
#                           rcond=None, full=False, w=1/data_mz['signal_std'], cov=True)
#     fit_y = np.polyval(poly_fit,time_numeric)
#     TT = np.vstack([time_numeric**(degree-i) for i in range(degree+1)]).T
#     C_y = np.dot(TT, np.dot(poly_cov, TT.T)) # C_y = TT*C_z*TT.T
#     sig_y = np.sqrt(np.diag(C_y))  # Standard deviations are sqrt of diagonal     
#     
#     cs = UnivariateSpline(x = time_numeric, y = data_mz['signal_mean'], w = 1/data_mz['signal_std'])
#                           
#   
#     plt.figure(figsize=(10,6))
#     plt.scatter(data_mz['TimeF'], data_mz['signal_mean'], label='data', zorder = 10, s = 10)
#     plt.errorbar(data = data_mz, x = 'TimeF', y = 'signal_mean', yerr = 'signal_std', label = "uncertainty",
#                  color = 'grey', elinewidth = 0.5, linewidth=0)
#     plt.plot(data_mz['TimeF'],fit_y,'r')
#     plt.fill_between(data_mz['TimeF'], fit_y+sig_y*3, fit_y-sig_y*3, alpha=.25, color = 'r')
#     plt.plot(data_mz['TimeF'], cs(time_numeric))
#     # plt.scatter(stats['TimeF'], stats['signal_mean'], label='mean')
#     # plt.plot(data['TimeF'], [stats['signal_mean']] * len(data), label='mean')
#     # plt.plot(data['TimeF'], [stats['signal_mean'] + stats['signal_std']] * len(data), label='std - upper')
#     # plt.plot(data['TimeF'], [stats['signal_mean'] - stats['signal_std']] * len(data), label='std - lower')
#     plt.xlabel('Time')
#     plt.ylabel('Signal')
#     plt.legend()
#     plt.title(f"{mz} - {type}")
#     plt.show()
#   
