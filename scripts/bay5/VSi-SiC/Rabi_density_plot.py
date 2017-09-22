#!/usr/bin/python
# -*- coding: utf-8 -*-
import cPickle, pylab, os
import matplotlib.pyplot as plt
import numpy as np
import re
import matplotlib, pandas, time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.pylab import *
import matplotlib.pyplot as plt
from tools import save_toolbox
from stat import S_ISREG, ST_CTIME, ST_MODE,ST_MTIME
from matplotlib.backends.backend_pdf import PdfPages




import threading, time
from tools import save_toolbox
from tools.emod import FreeJob
from traits.api import HasTraits, Trait, Instance, Property, Float,Int, Array,Bool, Array, String, Str, Enum, Button
from datetime import datetime, timedelta
from tools.utility import timestamp
import logging
import cPickle


"""HELPFULL METHODS"""
def get_file_list(some_path):
    orgin_path=os.getcwd()
    os.chdir(some_path)
    path = os.getcwd()
    file_list = [f for f in os.listdir('.') if os.path.isfile(f)]
    os.chdir(orgin_path)
    return file_list

def find_string_in_array(f,seq):
  """Return items containing f in sequence """
  regex=re.compile(".*"+f+".*")#using regulary expressions for finding part of string in a list
  item=[m.group(0) for l in seq for m in [regex.search(l)] if m]
  return item[0] #item itself is a 1d-array so [0] has to be added

def get_file_list_sorted(some_path):
    """sorts a list as a function of data change"""
    #NOTE: use `ST_MTIME` to sort by a modification date
    #NOTE: use `ST_CTIME` to sort by a creation date
    # path to the directory (relative or absolute)
    dirpath = some_path
    orgin_path=os.getcwd()
    date_list=[]
    new_file_list=[]
    
    # get all entries in the directory w/ stats
    entries = (os.path.join(dirpath, fn) for fn in os.listdir(dirpath))
    entries = ((os.stat(path), path) for path in entries)
    
    # leave only regular files, insert creation date
    entries = ((stat[ST_MTIME], path)
               for stat, path in entries if S_ISREG(stat[ST_MODE]))
    #NOTE: on Windows `ST_CTIME` is a creation date 
    #  but on Unix it could be something else
    #NOTE: use `ST_MTIME` to sort by a modification date
    for mdate, path in sorted(entries):
         if path.endswith('.pys'):
            date_list.append(time.ctime(mdate))
            new_file_list.append(os.path.basename(path))
    os.chdir(orgin_path)
    return new_file_list

def return_value(k, string):
    """Returns a wanted number in front of a string as an integer found in a given string k, output limited to 10 digits"""
    v_idx = k.find(string)
    #print 'string ',string,' found at position: ',v_idx
    returnvalue = 0
    for i in range(1,10):        
        try:
            returnvalue = float(k[v_idx-i:v_idx])
        except:
            pass
    return returnvalue

def dictToCsv(dictionary, filename):
    df=pandas.DataFrame(dict([ (k,pandas.Series(v)) for k,v in dictionary.iteritems() ]))        
    df.to_excel(filename,sheet_name='book1')

def getRabiData(filename):
    """returns wavelength and intensity, must be in the folder"""
    with open(filename,'rb') as f:
        d=cPickle.load(f)
        return d['measurement']['tau'], d['measurement']['spin_state']  

def getCalibration(folder):
    laserlist=get_file_list_sorted(folder)[::-1]
    i=len(laserlist)-1
    laser_array=[]
    angle_array=[]

    while i >= 0:
        filename=laserlist[i]
        os.chdir(folder)    
        wavelength, intensity = getRabiData(filename)
        #wavelength = wavelength-np.mean(wavelength[500:])#substract BG
        intensity = intensity - np.mean(intensity[500:])#substract BG
        mymax = np.max(intensity[0:200])
        angle = return_value(filename,'_degree_')

        if i==len(laserlist)-1:#do this for first dataset to make it stackable
            laser_array = mymax
            angle_array = angle
        else:
            laser_array = np.append(laser_array,mymax)
            angle_array = np.append(angle_array,angle)
        i -=1
    return laser_array,angle_array



#----- Folder Specs ------------------------

startfolder = os.getcwd()

#data_path = 'F:/Messdaten/2016-07-06/CREE/Spec_Polar_P70K_Emission_Quater_Waveplate/'
base= 'D:\\data\\2016\\'
date= '2016-10-17\\'
sample='CREE5E17\\'
meas = 'Rabi_160_190MHz_-22dBm\\'
#meas= 'Rabi_160_190MHz_-16dBm\\'
off = 'off\\'
data_path=base+date+sample+meas
data_path = 'D:/Research/WorkData/VSi_V1_rabi/Rabi_160_190MHz_-13dBm/'
#Off_resonance_path=base+date+sample+meas+off

#file_list = get_file_list_sorted(data_path)[::-1]
file_list = get_file_list_sorted(data_path)[::-1]
#file_list = get_file_list(data_path)[::-1]
#file_list_off_resonance = get_file_list_sorted(Off_resonance_path)[::-1]
#cal='laser_calibration 5.5K excitation_HWP 10d_degree BP_Filter'
#calibration_path= base+date_laser+sample+cal+'/First_Measurement/pys'
#laser,angle_laser = getCalibration(calibration_path)# our laser polarization
os.chdir(data_path)
#tau,spin_state = getRabiData(file_list[1])
#mydict={'wavelength':wavelength}

#get some basic measurement parameters, as temperature etc
#RF_Frequency = return_value(file_list[1],'MHz_')
MW_Power = return_value(file_list[1],'dBm')

#create empty array first
spin_state_array=[]
tau_array=[]
RF_Frequency_array=[]

#off Resonance Data
#os.chdir(Off_resonance_path)
#off=file_list_off_resonance[0]
#tau_off, spin_state_off = getRabiData(off) 
#--------Loop over all files-----------
i=len(file_list)-1
while i >= 0:
    os.chdir(data_path)
    g=file_list[i]
    tau, spin_state = getRabiData(g) 
    RF_Frequency = return_value(g,'MHz')

    #-----Data stacking
    if i==len(file_list)-1:#do this for first dataset to make it stackable
        #spin_state=spin_state-spin_state_off
        spin_state_array = spin_state
        tau_array=tau
        RF_Frequency_array=RF_Frequency

    else:
        #spin_state_array = np.vstack((spin_state_array,spin_state-spin_state_off))
        spin_state_array = np.vstack((spin_state_array,spin_state))
        tau_array=np.append(tau_array,tau)
        RF_Frequency_array=np.append(RF_Frequency_array,RF_Frequency)
    #mydict[angle] = intensity#you can also save the dictonary right away
    i -=1    


plt.legend(bbox_to_anchor=(.05, -0.2), loc=2, borderaxespad=0)
fig, ax = plt.subplots(nrows=1, figsize=(6,10))
#spin_state_array_tr=np.transpose(np.transpose(spin_state_array))

ax.imshow(spin_state_array, extent=[tau_array[0],tau_array[5062],RF_Frequency_array[60], RF_Frequency_array[0]],aspect='auto',cmap='viridis',interpolation='none', vmin=4e5, vmax=4.5e5)
ax.grid(True)

ax.set_ylabel('RF Frequency [MHz]')
ax.set_xlabel('MW Pulse Length [ns]')
ax.set_title('Rabi Oscillation Plot')
plt.colorbar(ax.imshow(spin_state_array, extent=[tau_array[0],tau_array[5062],RF_Frequency_array[60], RF_Frequency_array[0]],aspect='auto',cmap='viridis',interpolation='none',vmin=4e5, vmax=4.5e5))#cmap='viridis',clim=(0, 1e5)), interpolation='none',
#plt.show()
plt.savefig('Rabi_Oscilation_MHz_-22dBm.png')
plt.clf()




#DFT of the data
i=len(file_list)-1
while i >= 0:
    plt.clf()
    os.chdir(data_path)
    g=file_list[i]
    #wavelength, intensity = getRabiData(g)
    tau, spin_state = getRabiData(g) 
    RF_Frequency = return_value(g,'MHz')
    z=np.ones(10000)
    #spin_state= spin_state-np.mean(spin_state)
    spin_state=spin_state/np.median(spin_state)
    z[0:spin_state.size]=spin_state
    #z[0:tau.size]=tau
    #z=spin_state
    fft=np.fft.rfft(z,)#-np.median(spin_state[0:50])
    #fft=np.fft.rfft(spin_state)
    #freq=np.fft.rfftfreq(np.size(3),d=1./(1./12e-9))
    freq=np.fft.rfftfreq(np.size(z),d=1./(1./12e-9))
    plt.plot(freq[:fft[7:].size],np.abs(fft[7:]))
    plt.savefig('FFT_'+str(RF_Frequency)+'_MHz.png')
    print 'FFT Plot saved'+str(RF_Frequency)+'_MHz'

    #-----Data stacking
    if i==len(file_list)-1:#do this for first dataset to make it stackable
        fft_array = fft
        freq_array= freq
        #RF_Frequency_array=RF_Frequency

    else:
        fft_array = np.vstack((fft_array,fft))
        freq_array=np.append(freq_array,freq)
        #RF_Frequency_array=np.append(RF_Frequency_array,RF_Frequency)

    
    i -=1

#fft_array=fft_array[:,10:]
freq_array=freq_array/1e6

plt.legend(bbox_to_anchor=(.05, -0.2), loc=2, borderaxespad=0)
fig2, ax2 = plt.subplots(nrows=1, figsize=(6,10))
#spin_state_array_tr=np.transpose(np.transpose(spin_state_array))
#ax2.imshow(np.abs(fft_array)**2,extent=[freq_array[0],freq_array[305060],RF_Frequency_array[60],RF_Frequency_array[0]], aspect='auto', interpolation='none',cmap='viridis', clim=(0, 0.015e1), vmin=0, vmax=0.02e1)#cmap='viridis',clim=(0, 1e5), interpolation='none',
ax2.imshow(np.abs(fft_array),extent=[freq_array[0],freq_array[305060],RF_Frequency_array[60],RF_Frequency_array[0]], aspect='auto', interpolation='none',cmap='viridis', clim=(0, 0.015e1), vmin=0, vmax=0.042e1)#cmap='viridis',clim=(0, 1e5), interpolation='none',

ax2.grid(True)
#extent=[freq_array[0],freq_array[4115],RF_Frequency_array[48],RF_Frequency_array[0]]
ax2.set_ylabel('RF Frequency [MHz]')
ax2.set_xlabel('Frequency [MHz]')
ax2.set_title('Rabi Oscillation Plot (FFT)')


# freq_array, RF_Frequency_array = np.mgrid[0:59,0:59]
# mesh = ax2.pcolormesh(freq_array, RF_Frequency_array, np.abs(fft_array))
#plt.colorbar(ax2.imshow(np.abs(fft_array)**2,extent=[freq_array[0],freq_array[305060],RF_Frequency_array[60],RF_Frequency_array[0]], aspect='auto',interpolation='none', cmap='viridis', clim=(0, 0.015e1), vmin=0, vmax=0.02e1))#cmap='viridis',clim=(0, 1e5)), interpolation='none',
plt.colorbar(ax2.imshow(np.abs(fft_array),extent=[freq_array[0],freq_array[305060],RF_Frequency_array[60],RF_Frequency_array[0]], aspect='auto',interpolation='none', cmap='viridis', clim=(0, 0.015e1), vmin=0, vmax=0.042e1))#cmap='viridis',clim=(0, 1e5)), interpolation='none',
#plt.show()
plt.savefig('Rabi_Oscilation_FFT_-22dBm.png')
plt.clf()



drifts = np.array([[fft_array], [freq_array], [RF_Frequency_array]]) 
np.savetxt('./fft_array.csv', fft_array, fmt='%.5f' )
np.savetxt('./freq_array.csv', freq_array, fmt='%.5f' )
np.savetxt('./RF_Frequency_array.csv', RF_Frequency_array, fmt='%.5f' )


os.chdir(startfolder)
