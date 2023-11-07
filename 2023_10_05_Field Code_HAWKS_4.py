# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 10:46:05 2023

@author: jpadil10
"""
# =============================================================================
#                               Main Comments for future Work
# =============================================================================

import serial
import logging
import matplotlib.pyplot as plt
import pyrga
import pandas as pd
import time
import timer
from ftdi_serial import Serial
from vicivalve import VICI
import datetime
import termcolor

Serial.list_device_serials()


# =============================================================================
# Initializing RGA connection
# =============================================================================

if __name__ == "__main__":
    # initialize client with non-default noise floor setting
    RGA = pyrga.RGAClient("COM3" , noise_floor=7)                              # will probably need to change this to 6 or five, test in lab
    # check filament status and turn it on if necessary
    if not RGA.get_filament_status():
        RGA.turn_on_filament()
    # read partial pressures of air constituent

# =============================================================================
# Intializing Selector Valve connection
# =============================================================================

serial = Serial(device_serial='AI03K9K1', baudrate=9600)
valve = VICI(serial=serial, positions= 6)
print('Connection has been established with', serial)
print()
print()

# =============================================================================
# Initializing masses to be measured (with the pyrga library) 
# =============================================================================

FMASSES = {
        16: "CH4",
        18: "H2O",
        20: 'Ne 20',
        28: "N2",
        32: "O2",
        40: "Ar",
        44: 'CO2',
}

CDEMMASSES = {
         4: 'He',
        22: 'Ne 22',
        84: 'Kr', 
        132:'Xe', 
}


# =============================================================================
# Create Empty Dataframe for valve switch time results
# =============================================================================

Time_Series_Start= time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())

df = pd.DataFrame([Time_Series_Start])
df.columns = (['Timestamp'])
df.index = ['Time_Series_start_time']

# =============================================================================
# Valve switching functions
# =============================================================================

def switch_to_valve2():
    valve.switch_valve(2)
    t1 = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    time.sleep(5)
    print('\n','Position is in valve', valve._send_and_receive('CP'), '\n')
    return t1
    
def switch_to_valve1():
    valve.switch_valve(1)
    t2 = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    time.sleep(5)
    print('Position is in valve', valve._send_and_receive('CP'), '\n')
    return t2
   
def switch_to_valve4():
    valve.switch_valve(4)
    t3 = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    time.sleep(5)
    print('Position is in valve', valve._send_and_receive('CP'), 'the Air standard','\n')
    return t3
    
def switch_to_valve2B():
    valve.switch_valve(2)
    t4 = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    time.sleep(5)
    print('Position is in valve', valve._send_and_receive('CP'), 'B', '\n')
    return t4

def switch_to_valve6B():
    valve.switch_valve(6)
    t5 = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    time.sleep(5)
    print('Position is in valve', valve._send_and_receive('CP'), 'the Gas Standard','\n')
    return t5
    
def switch_to_valve2C():
    valve.switch_valve(2)
    t6 = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    time.sleep(5)
    print('Position is in valve', valve._send_and_receive('CP'), 'C','\n')
    return t6
    
def switch_to_valve6():
    valve.switch_valve(6)
    t7 = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    time.sleep(5)
    print('Position is in valve', valve._send_and_receive('CP'),'the Membrane', '\n')
    return t7

# =============================================================================
# Timed Interval Functions
# =============================================================================

def execute_pumpback_switchA():
    start_time = time.time()
    fpbackA = pd.DataFrame()
    
    while True:
        FpbackPartial_pressures = [] 
        CpbackPartial_pressures = []
                
        for m, i in FMASSES.items():
            
            ft0a = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            fpressrue= RGA.read_mass(m)
            print("partial pressure of element {} is {} Torr".format(i, fpressrue))
            print()
            FpbackPartial_pressures.append(fpressrue)
            
        for m, i in CDEMMASSES.items():
            
            ft0b = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            RGA._send_command("HV", 1400) 
            cpressrue= RGA.read_mass(m)/1000
            print("partial pressure of element {} is {} Torr".format(i, cpressrue)) 
            print()
            CpbackPartial_pressures.append(cpressrue)
            
        RGA._send_command("HV", 0)  
        ftdf = pd.DataFrame([ft0a])  
        ftdf.index = ['Valve2A']
        ftdf.columns = ['TimeF']
        
        funcpdf = pd.DataFrame([FpbackPartial_pressures], columns= list(FMASSES.values()))
        funcpdf.index = ['Valve2A']
        
        fjoin = ftdf.join(funcpdf)  
        
        ftdf2 = pd.DataFrame([ft0b])
        ftdf2.index = ['Valve2A']
        ftdf2.columns = ['TimeEM'] 
        
        funcempdf = pd.DataFrame([CpbackPartial_pressures], columns= list(CDEMMASSES.values()))
        funcempdf.index = ['Valve2A']
        
        fjoin2 = ftdf2.join(funcempdf)
        fjoin3 = fjoin.join(fjoin2)      
               
        fpbackA = fpbackA.append(fjoin3)
    
        current_time = time.time()
        elapsed_time = current_time - start_time
        
        print()
        print('Position is in Valve 2A','\n')
        
        if elapsed_time >= 180:
            break
        
    return  fpbackA

# blank can do 4.6 (5) scans in the 300 seconds
def execute_blank_scan_seconds():
   start_time = time.time()
   fblankdf = pd.DataFrame()
   
   while True:
       FblankPartial_pressures = [] 
       CblankPartial_pressures = []
              
       for m, i in FMASSES.items():
           
           ft1a = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
           fpressrue= RGA.read_mass(m)
           print("partial pressure of element {} is {} Torr".format(i, fpressrue))
           print()
           FblankPartial_pressures.append(fpressrue)
           
       for m, i in CDEMMASSES.items():
           
           ft1b = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
           RGA._send_command("HV", 1400) 
           cpressrue= RGA.read_mass(m)/1000
           print("partial pressure of element {} is {} Torr".format(i, cpressrue)) 
           print()
           CblankPartial_pressures.append(cpressrue)
                      
       RGA._send_command("HV", 0)
       ftdf = pd.DataFrame([ft1a])
       ftdf.index = ['Valve1']
       ftdf.columns = ['TimeF']
        
       funcbdf = pd.DataFrame([FblankPartial_pressures], columns= list(FMASSES.values()))
       funcbdf.index = ['Valve1']
       
       fjoin = ftdf.join(funcbdf)
       
       ftdf2 = pd.DataFrame([ft1b])
       ftdf2.index = ['Valve1']
       ftdf2.columns = ['TimeEM'] 
       
       funcembdf = pd.DataFrame([CblankPartial_pressures], columns= list(CDEMMASSES.values()))
       funcembdf.index = ['Valve1']
       
       fjoin2 = ftdf2.join(funcembdf)
       fjoin3 = fjoin.join(fjoin2)        
              
       fblankdf = fblankdf.append(fjoin3)
      
       current_time = time.time()
       elapsed_time = current_time - start_time
       
       print()
       print('Position is in Valve Blank','\n')
       
       if elapsed_time >= 300:
           break
       
   return  fblankdf

# Air can do 4.6 (5) scans in the 300 seconds
def execute_air_scan_seconds():
    start_time = time.time()
    fairdf = pd.DataFrame()
    
    while True:
        FairPartial_pressures = []
        CairPartial_pressures = []
        
        for m, i in FMASSES.items():
            
            ft2a = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            fpressrue= RGA.read_mass(m)
            print("partial pressure of element {} is {} Torr".format(i, fpressrue))
            print()
            FairPartial_pressures.append(fpressrue)
            
        for m, i in CDEMMASSES.items():
            
            ft2b = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            RGA._send_command("HV", 1400)
            cpressrue= RGA.read_mass(m)/1000
            print("partial pressure of element {} is {} Torr".format(i, cpressrue)) 
            print()
            CairPartial_pressures.append(cpressrue)
       
        RGA._send_command("HV", 0)
        ftdf = pd.DataFrame([ft2a])
        ftdf.index = ['Valve4']
        ftdf.columns = ['TimeF']
         
        funcadf = pd.DataFrame([FairPartial_pressures], columns= list(FMASSES.values()))
        funcadf.index = ['Valve4']
        
        fjoin = ftdf.join(funcadf)
        
        ftdf2 = pd.DataFrame([ft2b])
        ftdf2.index = ['Valve4']
        ftdf2.columns = ['TimeEM'] 
        
        funcemairdf = pd.DataFrame([CairPartial_pressures], columns= list(CDEMMASSES.values()))
        funcemairdf.index = ['Valve4']
        
        fjoin2 = ftdf2.join(funcemairdf)
        fjoin3 = fjoin.join(fjoin2)
                
        fairdf = fairdf.append(fjoin3)
        
        current_time = time.time()
        elapsed_time = current_time - start_time
        
        print()
        print('Position is in Valve Air','\n')
        
        if elapsed_time >= 300:
            break
        
    return  fairdf

def execute_pumpback_switchB():
    start_time = time.time()
    fpbackB = pd.DataFrame()
    
    while True:
        FpbackPartial_pressures = [] 
        CpbackPartial_pressures = []
                
        for m, i in FMASSES.items():
            
            ft3a = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            fpressrue= RGA.read_mass(m)
            print("partial pressure of element {} is {} Torr".format(i, fpressrue))
            print()
            FpbackPartial_pressures.append(fpressrue)
            
        for m, i in CDEMMASSES.items():
            
            ft3b = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            RGA._send_command("HV", 1400) 
            cpressrue= RGA.read_mass(m)/1000
            print("partial pressure of element {} is {} Torr".format(i, cpressrue)) 
            print()
            CpbackPartial_pressures.append(cpressrue)
            
        RGA._send_command("HV", 0)
        ftdf = pd.DataFrame([ft3a])  
        ftdf.index = ['Valve2B']
        ftdf.columns = ['TimeF']
        
        funcpdf = pd.DataFrame([FpbackPartial_pressures], columns= list(FMASSES.values()))
        funcpdf.index = ['Valve2B']
        
        fjoin = ftdf.join(funcpdf) 
        
        
        ftdf2 = pd.DataFrame([ft3b])
        ftdf2.index = ['Valve2B']
        ftdf2.columns =['TimeEM'] 
        
        funcempdf = pd.DataFrame([CpbackPartial_pressures], columns= list(CDEMMASSES.values()))
        funcempdf.index = ['Valve2B']
        
        fjoin2 = ftdf2.join(funcempdf)
        fjoin3 = fjoin.join(fjoin2)       
               
        fpbackB = fpbackB.append(fjoin3)
    
        current_time = time.time()
        elapsed_time = current_time - start_time
        
        print()
        print('Position is in Valve 2B','\n')
        
        if elapsed_time >= 180:
            break
        
    return  fpbackB

def execute_gas_standard_scan_seconds():
    start_time = time.time()
    fairdf = pd.DataFrame()
    
    while True:
        FairPartial_pressures = []
        CairPartial_pressures = []
        
        for m, i in FMASSES.items():
            
            ft2a = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            fpressrue= RGA.read_mass(m)
            print("partial pressure of element {} is {} Torr".format(i, fpressrue))
            print()
            FairPartial_pressures.append(fpressrue)
            
        for m, i in CDEMMASSES.items():
            
            ft2b = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            RGA._send_command("HV", 1400)
            cpressrue= RGA.read_mass(m)/1000
            print("partial pressure of element {} is {} Torr".format(i, cpressrue)) 
            print()
            CairPartial_pressures.append(cpressrue)
       
        RGA._send_command("HV", 0)
        ftdf = pd.DataFrame([ft2a])
        ftdf.index = ['Valve6B']
        ftdf.columns = ['TimeF']
         
        funcadf = pd.DataFrame([FairPartial_pressures], columns= list(FMASSES.values()))
        funcadf.index = ['Valve6B']
        
        fjoin = ftdf.join(funcadf)
        
        ftdf2 = pd.DataFrame([ft2b])
        ftdf2.index = ['Valve6B']
        ftdf2.columns = ['TimeEM'] 
        
        funcemairdf = pd.DataFrame([CairPartial_pressures], columns= list(CDEMMASSES.values()))
        funcemairdf.index = ['Valve6B']
        
        fjoin2 = ftdf2.join(funcemairdf)
        fjoin3 = fjoin.join(fjoin2)
                
        fstddf = fairdf.append(fjoin3)
        
        current_time = time.time()
        elapsed_time = current_time - start_time
        
        print()
        print('Position is in Valve Gas Standard','\n')
        
        if elapsed_time >= 300:
            break
        
    return  fstddf

def execute_pumpback_switchC():
    start_time = time.time()
    fpbackC = pd.DataFrame()
    
    while True:
        FpbackPartial_pressures = [] 
        CpbackPartial_pressures = []
                
        for m, i in FMASSES.items():
            
            ft3a = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            fpressrue= RGA.read_mass(m)
            print("partial pressure of element {} is {} Torr".format(i, fpressrue))
            print()
            FpbackPartial_pressures.append(fpressrue)
            
        for m, i in CDEMMASSES.items():
            
            ft3b = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            RGA._send_command("HV", 1400) 
            cpressrue= RGA.read_mass(m)/1000
            print("partial pressure of element {} is {} Torr".format(i, cpressrue)) 
            print()
            CpbackPartial_pressures.append(cpressrue)
            
        RGA._send_command("HV", 0)
        ftdf = pd.DataFrame([ft3a])  
        ftdf.index = ['Valve2C']
        ftdf.columns = ['TimeF']
        
        funcpdf = pd.DataFrame([FpbackPartial_pressures], columns= list(FMASSES.values()))
        funcpdf.index = ['Valve2C']
        
        fjoin = ftdf.join(funcpdf) 
        
        
        ftdf2 = pd.DataFrame([ft3b])
        ftdf2.index = ['Valve2C']
        ftdf2.columns =['TimeEM'] 
        
        funcempdf = pd.DataFrame([CpbackPartial_pressures], columns= list(CDEMMASSES.values()))
        funcempdf.index = ['Valve2C']
        
        fjoin2 = ftdf2.join(funcempdf)
        fjoin3 = fjoin.join(fjoin2)       
               
        fpbackC = fpbackC.append(fjoin3)
    
        current_time = time.time()
        elapsed_time = current_time - start_time
        
        print()
        print('Position is in Valve 2C','\n')
        
        if elapsed_time >= 180:
            break
        
    return  fpbackC
        
# blank can do 9.2 (10) scans in the 600 seconds
def execute_membrane_scan_seconds():
    start_time = time.time()
    fmemdf = pd.DataFrame()

    while True:
        FmemPartial_pressures = []
        CmemPartial_pressures = []
        
        for m, i in FMASSES.items():
            
            ft3a = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            fpressrue= RGA.read_mass(m)
            print("partial pressure of element {} is {} Torr".format(i, fpressrue))
            print()
            FmemPartial_pressures.append(fpressrue)

        for m, i in CDEMMASSES.items():
            
            ft3b = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            RGA._send_command("HV", 1400) 
            cpressrue= RGA.read_mass(m)/1000
            print("partial pressure of element {} is {} Torr".format(i, cpressrue)) 
            print()
            CmemPartial_pressures.append(cpressrue)
        
        RGA._send_command("HV", 0)
        ftdf = pd.DataFrame([ft3a])
        ftdf.index = ['Valve6']
        ftdf.columns = ['TimeF']
         
        funcmdf = pd.DataFrame([FmemPartial_pressures], columns= list(FMASSES.values()))
        funcmdf.index = ['Valve6']
        
        fjoin = ftdf.join(funcmdf)
        
        ftdf2 = pd.DataFrame([ft3b])
        ftdf2.index = ['Valve6']
        ftdf2.columns = ['TimeEM'] 
        
        funcemMEMdf = pd.DataFrame([CmemPartial_pressures], columns= list(CDEMMASSES.values()))
        funcemMEMdf.index = ['Valve6']
        
        fjoin2 = ftdf2.join(funcemMEMdf)
        fjoin3 = fjoin.join(fjoin2)
        
        fmemdf = fmemdf.append(fjoin3)
        
        current_time = time.time()
        elapsed_time = current_time - start_time
        
        print()
        print('Position is in Valve Membrane','\n')
        
        if elapsed_time >= 300:
            break
        
    return  fmemdf

# =============================================================================
# Create Empty Dataframe for each valve
# =============================================================================

fpback1 = pd.DataFrame()

fblank = pd.DataFrame()

fair = pd.DataFrame()

fpback2 = pd.DataFrame()

fstd = pd.DataFrame()

fpback3 = pd.DataFrame()

fmem = pd.DataFrame()

fmaster = pd.DataFrame()
    
# =============================================================================
# Main Loop 
# =============================================================================

print("Sequence will now commence")
print()
print("Warming up filament for 10 seconds")
print()
time.sleep(10)


# Two scans takes 22 minutes
# 134 cycles are needed to cover a 48 hour period
cycle_counter = 0 
while cycle_counter < 134:
    
    print('Currently in Cycle', cycle_counter)
          
    t1 = switch_to_valve2()
    Pback1 = execute_pumpback_switchA()    
        
    t2 = switch_to_valve1()
    Blank = execute_blank_scan_seconds()
        
    t3 = switch_to_valve4()
    Air = execute_air_scan_seconds()
    
    t4 = switch_to_valve2B()
    Pback2 = execute_pumpback_switchB()
    
    if cycle_counter == 0:
        t5 = switch_to_valve6B()
        Std = execute_gas_standard_scan_seconds()
        t6 = switch_to_valve2C()
        Pback3 = execute_pumpback_switchC()
    
    t7 = switch_to_valve6()
    Membrane = execute_membrane_scan_seconds()
    

# =============================================================================
# Create Self appending Dataframe for Valve switch times
# =============================================================================
              
    df2 = pd.DataFrame([t1,t2,t3,t4,t5,t6,t7])
    df2.index = ['Valve2', 'Valve1', 'Valve4', 'Valve2B', 'Valve4B', 'Valve2C','Valve6']
    df2.columns = (['Timestamp'])
    df = df.append(df2)
    df.to_csv("C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Valve time/valve_time.csv")

# =============================================================================
# Create Self appending Dataframe for Pumpback inlet partial pressure results
# =============================================================================
  
    fpback1 = fpback1.append(Pback1)
    fpback1['TimeF'] = pd.to_datetime(fpback1['TimeF'])
    fpback1['TimeEM'] = pd.to_datetime(fpback1['TimeEM'])
    fpback1.to_csv("C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Pump_back/Pump_back_A.csv")
    
    fpback2 = fpback2.append(Pback2)
    fpback2['TimeF'] = pd.to_datetime(fpback2['TimeF'])
    fpback2['TimeEM'] = pd.to_datetime(fpback2['TimeEM'])
    fpback2.to_csv("C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Pump_back/Pump_back_B.csv")
    
# =============================================================================
# Create Self appending Dataframe for Blank inlet partial pressure results
# =============================================================================
    
    fblank = fblank.append(Blank)
    fblank['TimeF'] = pd.to_datetime(fblank['TimeF'])
    fblank['TimeEM'] = pd.to_datetime(fblank['TimeEM'])
    fblank.to_csv("C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Blanks/Blank.csv")
    
# =============================================================================
# Create Self appending Dataframe for Air inlet partial pressure results
# =============================================================================
        
    fair = fair.append(Air)
    fair['TimeF'] = pd.to_datetime(fair['TimeF'])
    fair['TimeEM'] = pd.to_datetime(fair['TimeEM'])
    fair.to_csv("C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Air/Air.csv")

# =============================================================================
# Create Self appending Dataframe for Gas Standard inlet partial pressure results
# =============================================================================
   
    if cycle_counter == 0:
        fstd = fstd.append(Std)
        fstd['TimeF'] = pd.to_datetime(fstd['TimeF'])
        fstd['TimeEM'] = pd.to_datetime(fstd['TimeEM'])
        fstd.to_csv("C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Std/Std.csv")
        fpback3 = fpback3.append(Pback3)
        fpback3['TimeF'] = pd.to_datetime(fpback3['TimeF'])
        fpback3['TimeEM'] = pd.to_datetime(fpback3['TimeEM'])
        fpback3.to_csv("C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Pump_back/Pump_back_C.csv")
        
# =============================================================================
# Create Self appending Dataframe for Membrane inlet partial pressure results
# =============================================================================
    
    fmem = fmem.append(Membrane)
    fmem['TimeF'] = pd.to_datetime(fmem['TimeF'])
    fmem['TimeEM'] = pd.to_datetime(fmem['TimeEM'])
    fmem.to_csv("C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Membrane/Membrane.csv")


# =============================================================================
# Create Self appending Dataframe for Membrane inlet partial pressure results
# =============================================================================
   
    fmaster = fmaster.append(Pback1)
    fmaster = fmaster.append(Blank)
    fmaster = fmaster.append(Air)
    fmaster = fmaster.append(Pback2)
    if cycle_counter == 0:
        fmaster = fmaster.append(Std)
        fmaster = fmaster.append(Pback3)
    fmaster = fmaster.append(Membrane)
    fmaster['TimeF'] = pd.to_datetime(fmaster['TimeF'])
    fmaster['TimeEM'] = pd.to_datetime(fmaster['TimeEM'])
    fmaster.to_csv("C:/Users/jpadil10/OneDrive - Southern Cross University/Documents/AUSTRALIA/PHD/Data/Hawkesburry 4/MIMS/Surface Water/Valves/Master/master_dataframe.csv")
  
# =============================================================================
# Create Plot for visualization
# =============================================================================
  
    x = fmaster['TimeF']
    y_columns = ['CH4','H2O','N2','O2','Ar','CO2']
    y = fmaster[y_columns]

    plt.figure(figsize=(10,6))
    for column in y_columns:
        plt.plot(x, y[column], label=column)
    plt.xlabel('Time')
    plt.legend()
    plt.grid(True)
    plt.yscale('log')
    plt.ylim(1e-10, 1e-5)

    plt.show()
    
    xem = fmaster['TimeEM']
    yem_columns = ['He','Kr','Xe']
    yem = fmaster[yem_columns]
    plt.figure(figsize=(10,6))
    for column in yem_columns:
        plt.plot(xem, yem[column], label=column)
    plt.xlabel('Time')
    plt.legend()
    plt.grid(True)
    plt.yscale('log')
    plt.ylim(1e-11, 1e-5)

    plt.show()
# =============================================================================
# =============================================================================
    
    
    print()
    print('Cycle is:', cycle_counter)
    print()
    print('Pumpback A values:', fpback1)
    print()
    print('Blank values:', fblank)
    print()
    print('Air values:', fair)
    print()
    print('Pumpback B values:', fpback2)
    print()
    print('Membrane values:', fmem)
    print()

    print('Cycle:', cycle_counter, 'complete')
    cycle_counter = cycle_counter + 1
    print()
    
    
# =============================================================================
#  Final wrap ups  
# =============================================================================

RGA._send_command("HV",0)
# RGA.turn_off_filament()
# serial.disconnect()


