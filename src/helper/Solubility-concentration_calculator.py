# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 11:03:04 2023

@author: jpadil10
"""

# =============================================================================
"                          Soulubility - Henry Calculator                     "
# =============================================================================

import os
import pandas as pd
import datetime
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np

Membrane = pd.read_csv("Data/SCU POND/Membrane_sub2.csv")                                  # Using only the second Half of the time series 

temp = Membrane['water']
  
# temperature_range = range(0, 26)
# Temp_df= pd.DataFrame(temperature_range, columns=['temp'])        
# temp = Temp_df['temp']  
                                                        # Temperature in Celcius
patm = Membrane['atmos']                                                        # Pressure in hPa
pATM = patm * 0.000986923267                                                    # To convert pressure into atm

# =============================================================================
# Solubility equations used in (Maier et al., 2022) || Water vapor correction and Henry coefficient equations used in 'Noblefit'/'Ruedipy'
# =============================================================================
# In order to get exact values as Maier Substitute PA in H calculation for  1013.25 (except in CO2 and CH4)

def O2_calculation():                                                           # Source: Garcia and Gordon, 1992 https://doi.org/10.4319/lo.1992.37.6.1307 ||  Units:(micromol/kg)

    vol = 0.2095
    O2 = {'O2': [5.80818, 3.20684, 4.11890, 4.93845, 1.01567, 1.411575, -7.01211e-3, -7.25958e-3, -7.93334e-3, -5.4491e-3, -1.32412e-7]}
    O2_df = pd.DataFrame(O2)
    O2_df.index = ['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'B0', 'B1', 'B2', 'B3', 'C0']
    Concentration_df = pd.DataFrame(columns=['Temperature (°C)'])
    
    for t, pm, PA in zip(temp,pATM,patm):
        
        pw_hPa = 10 ** ((0.7859 + 0.03477 * t) / (1 + 0.00412 * t))             # Water vapor pressure in hPa            
        pw_atm = 10 ** ((0.7859 + 0.03477 * t) / (1 + 0.00412 * t)) / 1013.25   # Water vapor pressure converted to atm
        p = pm - pw_atm                                                         # PRESSURE DRY AIR
        S = 0                                                                   # Salinity in per mille (g/kg)
        Ts = (298.15 - t) / (273.15 + t)                                        
        Ts = np.log(Ts) 
        
        concentrations = {}
        for species in O2_df.columns:
            coef_series = O2_df[species]
            
            # Brenwald divides at the end by 1000 to convert from kg to g. || Also his sources give solubility in ml/kg, hence 1 ml = 1 cc 
            
            lnC = coef_series['A0'] + (coef_series['A1'] * Ts) + (coef_series['A2'] * Ts**2) + (coef_series['A3'] * Ts**3) + (coef_series['A4'] * Ts**4) + (coef_series['A5']* Ts**5) + S * ((coef_series['B0']) + (coef_series['B1'] * Ts) + (coef_series['B2'] * Ts**2) + (coef_series['B3'] * Ts**3)) + (coef_series['C0'] * S**2)       
            c = (np.exp(lnC)) * (p / (1 - pw_atm))                                  # original (p / (1 - pw_atm) / 1000)
            concentrations[species] = c

        #Henry coefficient  calculation

        print ('Pa is ', PA)
        Hcp = (concentrations['O2'] / ((PA - pw_hPa) * vol))                    # PA and pw_hPa are in hPa, therefore this is micromol/(kg*hPa)
        concentrations['Henry Solubility Constant [O2]'] = Hcp
        
        row_data = {'Temperature (°C)': t,'Patmos (atm)': pm, 'Patmos (hPa)' : PA , **concentrations}
        Concentration_df = Concentration_df.append(row_data, ignore_index=True)

    return Concentration_df

def N2_calculation():                                                           # Source: Hamme and Emerson, 2004 https://doi.org/10.1016/j.dsr.2004.06.009 ||  Units:(micromol/kg * atm)
    
    
    vol = 0.781
    N2 = {'N2': [6.42931, 2.92704, 4.32531, 4.69149, -7.44129e-3, -8.02566e-3, -1.46775e-2] }
    N2_df = pd.DataFrame(N2)
    N2_df.index = ['A0', 'A1', 'A2', 'A3', 'B0', 'B1', 'B2']
    Concentration_df = pd.DataFrame(columns=['Temperature (°C)'])
    
    for t, pm, PA in zip(temp,pATM,patm):
        
        pw_hPa = 10 ** ((0.7859 + 0.03477 * t) / (1 + 0.00412 * t))             # Water vapor pressure in hPa                             
        pw_atm = 10 ** ((0.7859 + 0.03477 * t) / (1 + 0.00412 * t)) / 1013.25   # Water vapor pressure converted to atm
        p = pm - pw_atm                                                         # PRESSURE DRY AIR
        S = 0                                                                   # Assuming a constant salinity value, This should be 0 in freshwater
        Ts = (298.15 - t) / (273.15 + t)
        Ts = np.log(Ts)
        
        concentrations = {}
        for species in N2_df.columns:
            coef_series = N2_df[species]
            lnC = coef_series['A0'] + (coef_series['A1'] * Ts) + (coef_series['A2'] * Ts**2) + (coef_series['A3'] * Ts**3) + S * ((coef_series['B0']) + (coef_series['B1'] * Ts) + (coef_series['B2'] * Ts**2))
            c = (np.exp(lnC)) * (p / (1 - pw_atm))                                  
            concentrations[species] = c
            
        #Henry coefficient  calculation

        Hcp = concentrations['N2'] / ((PA - pw_hPa) * vol) 
        concentrations['Henry Solubility Constant [N2]'] = Hcp
        
        row_data = {'Temperature (°C)': t,'Patmos (atm)': pm, **concentrations}
        Concentration_df = Concentration_df.append(row_data, ignore_index=True)
        
    return Concentration_df

def Ar_calculation():                                                           # Source: Hamme and Emerson, 2004 https://doi.org/10.1016/j.dsr.2004.06.009 ||  Units:(micromol/kg * atm)
       
    vol = 0.0094
    Ar = {'Ar': [2.7915, 3.17609, 4.13116, 4.90379, -6.96233e-3, -7.6667e-3, -1.16888e-2] }
    Ar_df = pd.DataFrame(Ar)
    Ar_df.index = ['A0', 'A1', 'A2', 'A3', 'B0', 'B1', 'B2']
    Concentration_df = pd.DataFrame(columns=['Temperature (°C)'])
    
    for t, pm, PA in zip(temp,pATM,patm):
        
        pw_hPa = 10 ** ((0.7859 + 0.03477 * t) / (1 + 0.00412 * t))             # Water vapor pressure in hPa 
        pw_atm = 10 ** ((0.7859 + 0.03477 * t) / (1 + 0.00412 * t)) / 1013.25   # Water vapor pressure converted to atm
        p = pm - pw_atm                                                         # PRESSURE DRY AIR
        S = 0                                                                   # Assuming a constant salinity value, This should be 0 in freshwater
        Ts = (298.15 - t) / (273.15 + t)
        Ts = np.log(Ts)
        
        concentrations = {}
        for species in Ar_df.columns:
            coef_series = Ar_df[species]
            lnC = coef_series['A0'] + (coef_series['A1'] * Ts) + (coef_series['A2'] * Ts**2) + (coef_series['A3'] * Ts**3) + S * ((coef_series['B0']) + (coef_series['B1'] * Ts) + (coef_series['B2'] * Ts**2))
            c = (np.exp(lnC)) * (p / (1 - pw_atm))                                  
            concentrations[species] = c
            
        #Henry coefficient  calculation

        Hcp = concentrations['Ar'] / ((PA - pw_hPa) * vol) 
        concentrations['Henry Solubility Constant [Ar]'] = Hcp
        
        row_data = {'Temperature (°C)': t,'Patmos (atm)': pm, **concentrations}
        Concentration_df = Concentration_df.append(row_data, ignore_index=True)
        
    return Concentration_df

def CO2_calculation():                                                          # Source Weiss, 1974 https://doi.org/10.1016/0304-4203(74)90015-2           ||  Units:(mol/kg * atm)
    
    # % Inverse Henry's Law constant, using formulae and data from review/compilation by Rolf Sander (1999), with T in Kelvin:
    # %   k°H = Inverse Henry's law constant for solubility in water at 298.15 K, in (mol/L)/atm
    # %   d(ln(kH))/d(1/T) = Temperature dependence constant, in K 
    # %   kH(T) = k°H * exp( d(ln(kH))/d(1/T) * (1/T - 1/298.15 K) )  Inverse Henry's Law coefficient, in (mol/L)/atm
    
    kH0 = np.mean([3.3e-2, 3.3e-2, 3.3e-2, 3.4e-2, 3.4e-2, 3.4e-2, 3.4e-2, 3.4e-2, 3.3e-2])
    uT  = np.mean([ 2400, 2400, 2300, 2300, 2400, 2300, 2400, 2400, 2400 ])
    vol = 0.0004
    
    CO2 = {'CO2': [-60.2409, 93.4517, 23.3585, 0.023517, -0.023656, 0.0047036]}
    CO2_df = pd.DataFrame(CO2)
    CO2_df.index = ['A1', 'A2', 'A3', 'B1', 'B2', 'B3']
    Concentration_df = pd.DataFrame(columns=['Temperature (°C)'])

    for t, pm, PA in zip(temp,pATM,patm):
        
        pw_hPa = 10 ** ((0.7859 + 0.03477 * t) / (1 + 0.00412 * t))             # Water vapor pressure in hPa                                   
        pw_atm = 10 ** ((0.7859 + 0.03477 * t) / (1 + 0.00412 * t)) / 1013.25   # Water vapor pressure in hPa converted to atm
        p = pm - pw_atm                                                         # PRESSURE DRY AIR
        S = 0                                                                   # salinity in per mille (g/kg)
        T = t + 273.15

        concentrations = {}
        for species in CO2_df.columns:
            coef_series = CO2_df[species]
            lnKo = coef_series['A1'] + (coef_series['A2'] *(100/T)) + (coef_series['A3'] * np.log(T/100)) + S * ((coef_series['B1']) + (coef_series['B2'] * (T/100)) + (coef_series['B3'] * (T/100)**2))
            C = (np.exp(lnKo)) * (p / (1 - pw_atm))  * 1000000                  # Conversion to micromol  
            concentrations[species] = C
        
        #Henry coefficient  calculation (Taken directly from Brenwalds code)
# =============================================================================
#         try Inserting other equation and compare differences. 
# =============================================================================
        
        kH  = (kH0 * np.exp( uT * ((1/288.15) - (1/298.15))))* 1013.25
        concentrations['Henry Solubility Constant [CO2]'] = kH
        
        row_data = {'Temperature (°C)': t,'Patmos (atm)': pm, **concentrations}
        Concentration_df = Concentration_df.append(row_data, ignore_index=True)

    return Concentration_df

def CH4_calculation():                                                          # Source Wiesenburg and Guinasso Jr, 1979 https://doi.org/10.1021/je60083a006 ||  Units:(nanomol/kg * atm)
   
    # % Inverse Henry's Law constant, using formulae and data from review/compilation by Rolf Sander (1999), with T in Kelvin:
    # %   k°H = Inverse Henry's law constant for solubility in water at 298.15 K, in (mol/L)/atm
    # %   d(ln(kH))/d(1/T) = Temperature dependence constant, in K 
    # %   kH(T) = k°H * exp( d(ln(kH))/d(1/T) * (1/T - 1/298.15 K) )  Inverse Henry's Law coefficient, in (mol/L)/atm
    
    kH0 = np.mean([ 1.4e-3, 1.4e-3, 1.4e-3, 1.4e-3, 1.5e-3, 1.4e-3])            # Values taken from Sander (version 4)
    uT  = np.mean([ 1900, 1600, 1600, 1500, 1600, 1700 ]) 
    vol = 0.00000184
    
    CH4 = {'CH4': [-417.5053, 599.8626, 380.3636, -62.0764, -0.064236, 0.034980, -0.0052732] }
    CH4_df = pd.DataFrame(CH4)
    CH4_df.index = ['A1', 'A2', 'A3', 'A4', 'B1', 'B2', 'B3']
    Concentration_df = pd.DataFrame(columns=['Temperature (°C)'])
    
    for t in temp:
        S = 0  
        fg = 1.41e-6  # fg = atmospheric concentration 
        T = t + 273.15                                                          # temperature (T) in Kelvin
        
        concentrations = {}
        for species in CH4_df.columns:
            coef_series = CH4_df[species]
            lnC = np.log(fg) + coef_series['A1'] + (coef_series['A2'] *(100/T)) + (coef_series['A3'] * np.log(T/100)) + (coef_series['A4'] * (T/100)) + S * ((coef_series['B1']) + (coef_series['B2'] * (T/100)) + (coef_series['B3'] * (T/100)**2))
            C = (np.exp(lnC))/ 1000                                             #conversion to micromol
            concentrations[species] = C
        
        
        #Henry coefficient  calculation (Taken directly from Brenwalds code)

        kH  = (kH0 * np.exp( uT * ((1/288.15) - (1/298.15))))* 1013.25
        concentrations['Henry Solubility Constant [CH4]'] = kH
                
        row_data = {'Temperature (°C)': t, **concentrations}
        Concentration_df = Concentration_df.append(row_data, ignore_index=True)
        Concentration_df['CH4'] = Concentration_df['CH4'] 

    return Concentration_df


def He_calculation():                                                           # Source Weiss, 1971 https://doi.org/10.1021/je60049a019                    ||  Units:(mL/kg * atm)
    
    vol = 0.00000554
    He = {'He': [-152.9405, 196.8840, 126.8015, -20.6767, -0.040543, 0.021315, -0.0030732] }
    He_df = pd.DataFrame(He)
    He_df.index = ['A1', 'A2', 'A3', 'A4', 'B1', 'B2', 'B3']
    Concentration_df = pd.DataFrame(columns=['Temperature (°C)'])
    
    for t, pm, PA in zip(temp,pATM,patm):
        
        pw_hPa = 10 ** ((0.7859 + 0.03477 * t) / (1 + 0.00412 * t))             # Water vapor pressure in hPa                              
        pw_atm = 10 ** ((0.7859 + 0.03477 * t) / (1 + 0.00412 * t)) / 1013.25   # Water vapor pressure converted to atm
        p = pm - pw_atm                                                         # PRESSURE DRY AIR
        S = 0                                                                   # Salinity in per mille (g/kg)
        T = t + 273.15 
        
        concentrations = {}
        for species in He_df.columns:
            coef_series = He_df[species]
            lnC = coef_series['A1'] + (coef_series['A2'] *(100/T)) + (coef_series['A3'] * np.log(T/100)) + (coef_series['A4'] * (T/100)) + S * ((coef_series['B1']) + (coef_series['B2'] * (T/100)) + (coef_series['B3'] * (T/100)**2))
            c = (np.exp(lnC)) * (p / (1 - pw_atm))
            c = (c/22414.1) * 1000000                                           #conversion to mol (by division) and then to micromol 
            concentrations[species] = c
            
        #Henry coefficient  calculation

        Hcp = concentrations['He'] / ((PA - pw_hPa) * vol)
        concentrations['Henry Solubility Constant [He]'] = Hcp
        
        row_data = {'Temperature (°C)': t,'Patmos (atm)': pm, **concentrations}
        Concentration_df = Concentration_df.append(row_data, ignore_index=True)

    return Concentration_df

def Ne_calculation():                                                           # Source: Hamme and Emerson, 2004 https://doi.org/10.1016/j.dsr.2004.06.009 ||  Units:(micromol/kg * atm)
   
    vol = 0.00001818
    Ne = {'Ne': [2.18156, 1.29108, 2.12504, 0, -5.94737e-3, -5.13896e-3, 0]}
    Ne_df = pd.DataFrame(Ne)
    Ne_df.index = ['A0', 'A1', 'A2', 'A3', 'B0', 'B1', 'B2']
    Concentration_df = pd.DataFrame(columns=['Temperature (°C)'])
    
    # Loop through the temperature range and calculate concentrations
    for t, pm, PA in zip(temp,pATM,patm):
        pw_hPa = 10 ** ((0.7859 + 0.03477 * t) / (1 + 0.00412 * t))             # Water vapor pressure in hPa                                   
        pw_atm = 10 ** ((0.7859 + 0.03477 * t) / (1 + 0.00412 * t)) / 1013.25   # Water vapor pressure converted to atm
        p = pm - pw_atm                                                         # PRESSURE DRY AIR
        S = 0                                                                   # Assuming a constant salinity value, This should be 0 in freshwater
        Ts = (298.15 - t) / (273.15 + t)
        Ts = np.log(Ts)
        
        concentrations = {}
        for species in Ne_df.columns:
            coef_series = Ne_df[species]
            lnC = coef_series['A0'] + (coef_series['A1'] * Ts) + (coef_series['A2'] * Ts**2) + (coef_series['A3'] * Ts**3) + S * ((coef_series['B0']) + (coef_series['B1'] * Ts) + (coef_series['B2'] * Ts**2))
            c = (np.exp(lnC)) * (p / (1 - pw_atm))
            concentrations[species] = c
            
        #Henry coefficient  calculation

        Hcp = concentrations['Ne'] / ((PA - pw_hPa) * vol)
        concentrations['Henry Solubility Constant [Ne]'] = Hcp
        
        row_data = {'Temperature (°C)': t,'Patmos (atm)': pm, **concentrations}
        Concentration_df = Concentration_df.append(row_data, ignore_index=True)
        
    return Concentration_df
# =============================================================================
# For storing solubilities/concentration results
# =============================================================================

def Solubility_merging():
    
    All_gases_df = pd.DataFrame(columns=['Temperature (°C)'])  
    All_gases_df = pd.merge(All_gases_df, O2_calculation(), how='right')
    All_gases_df = pd.merge(All_gases_df, N2_calculation(), how='right')
    All_gases_df = pd.merge(All_gases_df, Ar_calculation(), how='right')
    All_gases_df = pd.merge(All_gases_df, CO2_calculation(), how='right')
    All_gases_df = pd.merge(All_gases_df, CH4_calculation(), how='right')
    All_gases_df = pd.merge(All_gases_df, He_calculation(), how='right')  
    All_gases_df = pd.merge(All_gases_df, Ne_calculation(), how='right')
    All_gases_df.insert(0, "Time", Membrane["TimeF"])
    All_gases_df.set_index('Time', inplace=True)
    
    return All_gases_df

All_gases_df = Solubility_merging()

# =============================================================================
# Converting Partial pressures to hPa
# =============================================================================

def convert2hPa():
    
    gases = ['CH4', 'H2O', 'Neon 20', 'N2', 'O2', 'Ar', 'CO2','He', 'Neon 22', 'Kr', 'Xe']
    converted_gases = pd.DataFrame()
    
    for gas in gases:
        converted_values = Membrane[gas] * 1.3332236842                         # Convert from Torr to hPa
        converted_gases[gas] = converted_values
        
    converted_gases = converted_gases.set_index(All_gases_df.index)
    return converted_gases

converted_gases_df = convert2hPa()

def Concentration_calculation():
    
    result_df = pd.DataFrame()
    column_mapping = {
    
    'O2': 'Henry Solubility Constant [O2]',
    'N2': 'Henry Solubility Constant [N2]',
    'Ar': 'Henry Solubility Constant [Ar]',
    'CH4': 'Henry Solubility Constant [CH4]',
    'CO2': 'Henry Solubility Constant [CO2]',
    'He': 'Henry Solubility Constant [He]',
      }
   
    for gas_column, henry_column in column_mapping.items():
    
        result = converted_gases_df[gas_column] * All_gases_df[henry_column]
        result_df[f'{gas_column} Concentration'] = result
        
    return result_df
    
Concentrations = Concentration_calculation()
Concentrations = Concentrations.reset_index()

x = Concentrations['Time']
y_columns = ['O2 Concentration']
y = Concentrations[y_columns]
plt.figure(figsize=(10,6))
for column in y_columns:
    plt.scatter(x, y[column], label=column)
plt.xlabel('Time')
plt.legend()
plt.grid(True)
plt.yscale('log')
plt.ylim(1e-10, 1e-5)

plt.show()

    
    
    
    
    
    
