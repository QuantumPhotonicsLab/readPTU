import cPickle
import os
import numpy as np
import pylab as plt
from analysis_simulations.libs.tools import toolbox
from matplotlib.colors import ListedColormap

reload (toolbox)
_viridis_data = np.load ('D:/Research/bonato-lab/analysis_simulations/viridis_cmp.npy')
viridis = ListedColormap(_viridis_data, name='viridis')
plt.register_cmap(name='viridis', cmap=viridis)

class rabi_analysis_VSi_SiC ():

	def __init__(self):
		self.folder = 'D:/Research/WorkData/VSi_V1_rabi'

	def getRabiData(self, filename):
	    with open(filename,'rb') as f:
	        d=cPickle.load(f)
	        return d['measurement']['tau'], d['measurement']['spin_state']  

	def load_data (self, dBm):

		data_path=self.folder + '/Rabi_160_190MHz_-'+str(45-dBm)+'dBm/'
		file_list = toolbox.get_file_list(data_path)[::-1]

		os.chdir(data_path)

		self.data_dict = {}
		self.RF_frq_array = np.zeros (len(file_list))

		i = 0
		for g in file_list:
		    tau, spin_state = self.getRabiData(g) 
		    RF_frq = toolbox.return_value(g,'MHz')
		    extr_dbm = toolbox.return_value(g,'dBm')
		    self.RF_frq_array [i] = RF_frq
		    self.data_dict[str(RF_frq)+'MHz'] = {}
		    self.data_dict[str(RF_frq)+'MHz']['tau'] = tau
		    self.data_dict[str(RF_frq)+'MHz']['PL'] = spin_state
		    i +=1

		self.RF_frq_array = np.sort(self.RF_frq_array)
		self.data_dict ['RF_frq'] = self.RF_frq_array

	def plot_rabi (self, frq='all', do_plot=True):
		if frq=='all':
			i = 0
			for k in self.data_dict['RF_frq']:
				if i==0:
					tau = self.data_dict[str(k)+'MHz']['tau']
					all_PL = np.zeros((len(self.data_dict['RF_frq']), len(tau)))
				all_PL [i, :] = self.data_dict[str(k)+'MHz']['PL']
				i += 1
			if do_plot:
				plt.figure(figsize = (10,5))
				X, Y = plt.meshgrid (tau*1e-3, self.data_dict['RF_frq'])
				plt.pcolor (X, Y, all_PL, cmap= 'viridis')
				plt.xlabel ('time [us]', fontsize=16)
				plt.ylabel ('driving frq [MHz]', fontsize=16)
				plt.show()
		else:		
			t_us = self.data_dict[str(frq)+'MHz']['tau']*1e-3
			PL = self.data_dict [str(frq)+'MHz']['PL']*1e-3
			PL = PL - PL[0]
			PL = PL/np.median(PL)
			plt.figure(figsize = (15,5))
			plt.plot (t_us, PL, 'crimson')
			plt.plot (t_us, PL, 'o', color='crimson')
			plt.xlabel ('time [us]', fontsize=15)
			plt.ylabel ('PL [kcounts]', fontsize = 15)
			plt.show()
			return t_us, PL

	def return_rabi (self, frq):
		t = self.data_dict[str(frq)+'MHz']['tau']
		PL = self.data_dict [str(frq)+'MHz']['PL']*1e-3
		return t, PL		

	def plot_rabi_fft (self, cristian = True):
		i = 0
		for k in self.data_dict['RF_frq']:
			if i==0:
				tau = self.data_dict[str(k)+'MHz']['tau']
				N = len(tau)
				if cristian:
					PL_fft = np.zeros((len(self.data_dict['RF_frq']), N/2+1-3))
				else:
					PL_fft = []
			PL = self.data_dict[str(k)+'MHz']['PL']
			#PL = PL-PL[0]
			#PL = PL/np.median(PL)
			if cristian:
				fff = np.fft.fftshift(np.abs(np.fft.fft(PL))**2)
				fff = fff/max(fff)
				fff = fff[len(tau)/2+3:]
				PL_fft [i, :] = fff
			else:
				PL = PL/np.median(PL)
				fff = np.fft.rfft(PL)
				fff = fff[2:]
				if i==0:
					PL_fft = fff 
				else:
					PL_fft = np.vstack ((PL_fft, fff))
			i += 1


		df = 1./max(tau*1e-9)
		freq = np.fft.ifftshift(np.fft.fftfreq(PL.size))

		r_frq = freq*df*len(tau)
		r_frq = r_frq[len(tau)/2+3:]
		X, Y = plt.meshgrid (r_frq*1e-6, self.data_dict['RF_frq'])

		do_renorm = 0
		if do_renorm:
			a,b = np.shape (PL_fft)
			for i in np.arange(b):
				y = PL_fft[:, i]
				PL_fft[:, i] = y/np.sum(y)

		plt.figure(figsize = (4,9))
		if cristian:
			plt.pcolor (X, Y, PL_fft, cmap = 'viridis', vmin = np.min(np.min(PL_fft)), vmax= self.v_factor*np.max(np.max(PL_fft)))
		else:
			plt.imshow(np.abs(PL_fft), aspect='auto', interpolation='none',cmap='viridis', clim=(0, 0.015e1), vmin=0, vmax=0.042e1)

		#plt.imshow(PL_fft, cmap = 'viridis')
		plt.xlabel ('Rabi frq [MHz]', fontsize = 18)
		plt.ylabel ('driving frq [MHz]', fontsize = 18)
		plt.axis([0, 30, 160, 190])
		plt.show()
