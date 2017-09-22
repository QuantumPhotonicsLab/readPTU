import numpy as np
import pylab as plt
import os
import matplotlib
import lmfit

from analysis_simulations.libs.spin import rabi
from analysis_simulations.libs.spin import V1_VSi_SiC_spin_sims as V1
from scipy import signal

reload (rabi)
reload (V1)

matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18)

class fit_Rabi ():

	def __init__(self, power):
		self.power = power

	def load_data (self):
		self.analysis = rabi.rabi_analysis_VSi_SiC()
		self.analysis.load_data(self.power)

	def return_driving_frq (self):
		return self.analysis.data_dict ['RF_frq']

	def get_data_driving_frq (self, driving_frq):
		self.driving_frq = driving_frq
		self.t_exp = 1e-9*self.analysis.data_dict[str(driving_frq)+'MHz']['tau']
		PL = self.analysis.data_dict[str(driving_frq)+'MHz']['PL']
		PL = PL-PL[0]
		self.PL_exp = PL/np.median(PL)
		return self.t_exp, self.PL_exp

	def plot_data (self):
		plt.plot (self.t_exp*1e6, self.PL_exp, '--', color='crimson')
		plt.plot (self.t_exp*1e6, self.PL_exp, 'o', color='crimson')
		plt.xlabel ('time [us]', fontsize=18)
		plt.ylabel ('normalized PL', fontsize=18)
		plt.show()


class fit_Rabi_time_domain (fit_Rabi):

	def __init__(self, power, B_Gauss):
		self.power = power
		self.B_Gauss = B_Gauss

	def residual(self, params, x, data):
	    O = params['Omega']
	    t1 = params['t1']
	    t2 = params['t2']
	    I12 = params['I12']
	    PL_incr = params['PL_incr']
	    init_pol = params['init_pol']

	    sim = V1.V1_SiC (B_Gauss = self.B_Gauss, T2_star = 10000, t=x, verbose = False)
	    sim.set_intensity_transitions(I_m32_m12=1, I_m12_p12=I12, I_p12_p32=1)
	    sim.set_ODMR_pars (polariz_array = np.array([init_pol,1-init_pol,1-init_pol,init_pol]), 
	    		ODMR_contrast_array = np.array([1.,PL_incr,PL_incr,1.]))
	    sim.set_decay (t1 = t1, t2 = t2)
	    model = sim.rabi (f = self.driving_frq*1e6, Omega=O, do_plot = False)
	    return np.abs(data-model)

	def do_fit(self, params):

		data = self.PL_exp
		x = self.t_exp
		self.out = lmfit.minimize(self.residual, params, args=(x, data))
		lmfit.printfuncs.report_fit(self.out.params, min_correl=0.5)

		self.t1_fit = self.out.params['t1'].value
		self.t2_fit = self.out.params['t2'].value
		self.Omega_fit = self.out.params['Omega'].value
		self.I12_fit = self.out.params['I12'].value
		self.PL_incr_fit = self.out.params['PL_incr'].value
		self.init_pol_fit = self.out.params['init_pol'].value

	def plot_fit(self):

		t = np.linspace (0,10000,5000)*1e-9
		sim = V1.V1_SiC (B_Gauss = self.B_Gauss, T2_star = 1000, t=t)
		sim.set_intensity_transitions(I_m32_m12=1, I_m12_p12=self.I12_fit, I_p12_p32=1)
		sim.set_ODMR_pars (polariz_array = np.array([self.init_pol_fit,1-self.init_pol_fit, 1-self.init_pol_fit, self.init_pol_fit]), 
					ODMR_contrast_array = np.array([1.,self.PL_incr_fit, self.PL_incr_fit,1.]))
		sim.set_decay (t1 = self.t1_fit, t2 = self.t2_fit)
		PL_th = sim.rabi (f = self.driving_frq*1e6, Omega=self.Omega_fit, do_plot = False)
		self.t = t
		self.PL_th = PL_th

		plt.figure()
		plt.plot (t*1e6, PL_th, color='RoyalBlue', linewidth = 3)
		plt.plot (self.t_exp*1e6, self.PL_exp, 'o', color='crimson')
		plt.axis([0,1,0,2])
		plt.xlabel ('time [us]', fontsize=18)
		plt.ylabel ('normalized PL', fontsize=18)
		plt.show()

class fit_Rabi_all (fit_Rabi):

	def __init__(self, power):
		self.power = power

	def load_data (self):
		self.analysis = rabi.rabi_analysis_VSi_SiC()
		self.analysis.load_data(self.power)

	def residual(self, params):
	    O = params['Omega']
	    t1 = params['t1']
	    t2 = params['t2']
	    I12 = params['I12']
	    PL_incr = params['PL_incr']

	    ind = 0
	    for driving_frq in self.analysis.data_dict ['RF_frq']:
	    	t = self.analysis.data_dict[str(driving_frq)+'MHz']['tau']
	    	data = self.analysis.data_dict[str(driving_frq)+'MHz']['PL']
	    	data = data-data[0]
	    	data = data/np.median(data)
	    	t = t*1e-9

	    	sim = V1.V1_SiC (B_Gauss = 20, T2_star = 1e-3, t=t, verbose = False)
	    	sim.set_intensity_transitions(I_m32_m12=1, I_m12_p12=I12, I_p12_p32=1)
	    	sim.set_ODMR_pars (polariz_array = np.array([.4,0.1,0.1,0.4]), ODMR_contrast_array = np.array([1.,PL_incr,PL_incr,1.]))
	    	sim.set_decay (t1 = t1, t2 = t2)
    		model = sim.rabi (f = driving_frq*1e6, Omega=O, do_plot = False)
    		if (ind==0):
    			error = np.abs(model-data)
    		else:
    			error = error + np.abs(model-data)
    		ind += 1

	    return error

	def do_fit(self):

		params = lmfit.Parameters()
		params.add('Omega', value=2.48e6, min = 0.5e6, max=10e6, vary = False)
		params.add('t1', value=100e-9, min = 10e-9, max = 500e-9, vary = False)
		params.add('t2', value=300e-9, min = 200e-9, max = 2000e-9, vary = False)
		params.add('I12', value=1.4, min = 0.8, max= 2., vary = True)
		params.add('PL_incr', value = 1.2, min = 0.5, max= 3., vary = False)

		self.out = lmfit.minimize(self.residual, params)
		lmfit.printfuncs.report_fit(self.out.params, min_correl=0.5)

		self.t1_fit = self.out.params['t1'].value
		self.t2_fit = self.out.params['t2'].value
		self.Omega_fit = self.out.params['Omega'].value
		self.I12_fit = self.out.params['I12'].value
		self.PL_incr_fit = self.out.params['PL_incr'].value



class fit_Rabi_frq_domain (fit_Rabi):

	def __init__(self):
		pass

	def find_peaks (self, data):
		peakind = signal.find_peaks_cwt(data, np.arange(1,4))
		print peakind
		return peakind

	def residual(self, params, x, data):
	    O = params['Omega']
	    T2 = params['T2']
	    I12 = params['I12']
	    PL_incr = params['PL_incr']

	    sim = V1.V1_SiC (B_Gauss = 20, T2_star = T2, t=x, verbose = False)
	    sim.set_intensity_transitions(I_m32_m12=1, I_m12_p12=I12, I_p12_p32=1)
	    sim.set_ODMR_pars (polariz_array = np.array([.4,0.1,0.1,0.4]), ODMR_contrast_array = np.array([1.,PL_incr,PL_incr,1.]))
	    model = sim.rabi (f = self.driving_frq*1e6, Omega=O, do_plot = False)
	    f_mod, model_fft = self.calc_fft (x = x, y = model, in_fft = 2)
	    f_exp, data_fft = self.calc_fft (x = x, y = data, in_fft = 2)

	    return np.abs(data_fft-model_fft)

	def calc_fft (self, x, y, in_fft = 0):
		N = len(x)
		fff = np.fft.ifftshift(np.abs(np.fft.fft(y, N))**2)
		fff = fff[N/2+in_fft:]
		fff = fff/np.sum(fff)

		self.df = 1./max(x)
		freq = np.fft.ifftshift(np.fft.fftfreq(N))
		r_frq = freq*self.df*N
		r_frq = r_frq[N/2+in_fft:]
		return r_frq, fff

	def find_peaks_exp (self, power, driving_frq):
		a = rabi.rabi_analysis_VSi_SiC()
		a.load_data(power)
		t_us, PL_exp = a.plot_rabi(driving_frq)
		t_exp = t_us*1e-6
		f, PL_fft = self.calc_fft (x=t_exp, y = PL_exp, in_fft = 2)
		ind = self.find_peaks (PL_fft)

		plt.figure (figsize = (20,5))
		plt.plot (f*1e-6, PL_fft, 'royalblue')
		plt.plot (f[ind]*1e-6, PL_fft[ind], 'o', color = 'crimson')
		plt.axis ([0, 30, 0, 0.5])
		plt.xlabel ('Rabi frequency [MHz]', fontsize = 18)
		plt.show()

	def do_fit(self, power, driving_frq):
		a = rabi.rabi_analysis_VSi_SiC()
		a.load_data(power)
		t_us, PL_exp = a.plot_rabi(driving_frq)
		t_exp = t_us*1e-6
		self.t_exp = t_exp
		self.PL_exp= PL_exp
		self.driving_frq = driving_frq

		params = lmfit.Parameters()
		params.add('Omega', value=2.81e6, min = 0.5e6, max=10e6, vary = True)
		params.add('T2', value=400e-9, min = 150e-9, max = 10e-6, vary = False)
		params.add('I12', value=1.5, min = 0., max= 5., vary = True)
		params.add('PL_incr', value = 1.1, min = 0.5, max= 3., vary = True)

		self.out = lmfit.minimize(self.residual, params, args=(t_exp, PL_exp))
		lmfit.printfuncs.report_fit(self.out.params, min_correl=0.5)

		self.T2_fit = self.out.params['T2'].value
		self.Omega_fit = self.out.params['Omega'].value
		self.I12_fit = self.out.params['I12'].value
		self.PL_incr_fit = self.out.params['PL_incr'].value


	def plot_fit(self):
		t = np.linspace (0,10000,2000)*1e-9
		sim = V1.V1_SiC (B_Gauss = 63, T2_star = 400e-9, t=t)
		sim.set_intensity_transitions(I_m32_m12=1, I_m12_p12=self.I12_fit, I_p12_p32=1)
		sim.set_ODMR_pars (polariz_array = np.array([1.,0.0,0.0,1.]), ODMR_contrast_array = np.array([1.,self.PL_incr_fit, self.PL_incr_fit,1.]))
		PL_th = sim.rabi (f = self.driving_frq*1e6, Omega=self.Omega_fit, do_plot = True)

		self.f_fit, self.y_fit = self.calc_fft (x=t, y=PL_th, in_fft=2)
		self.f_exp, self.y_exp = self.calc_fft (x = self.t_exp, y = self.PL_exp, in_fft=2)

		plt.figure(figsize = (20,5))
		plt.plot (self.f_fit*1e-6, self.y_fit, 'royalblue')
		plt.plot (self.f_exp*1e-6, self.y_exp, 'o', color='crimson')
		plt.show()
