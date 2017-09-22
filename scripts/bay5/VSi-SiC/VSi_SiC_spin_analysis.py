import numpy as np
import pylab as plt
import os
import matplotlib
import lmfit

from analysis_simulations.libs.spin import rabi
print rabi
from analysis_simulations.libs.spin import V1_VSi_SiC_spin_sims as V1
from analysis_simulations.libs.spin import fit_rabi 
from scipy import signal

reload (rabi)
reload (V1)
reload (fit_rabi)

matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18)

def plot_2D_rabi():
	a = rabi.rabi_analysis_VSi_SiC()
	a.load_data(26)
	a.v_factor =1
	a.plot_rabi_fft()

def simulate_2D_Rabi():

	#I12_fit = 0.95
	t1 = 70e-9
	t2 = 200e-9
	Omega = 2.e6

	t = np.linspace (0,2000,1000)*1e-9
	sim = V1.V1_SiC (B_Gauss = 63, T2_star = 10000, t=t)
	sim.set_intensity_transitions(I_m32_m12=1, I_m12_p12=1.1, I_p12_p32=1)
	sim.set_ODMR_pars (polariz_array = np.array([1.,0.0,0.0,1.]), ODMR_contrast_array = np.array([1.,1.5, 1.5,1.]))
	sim.set_decay (t1 = t1, t2 = t2)
	sim.rabi_sweep_drive (init_frq=160e6, end_frq=190e6, nr_steps=100, Omega=Omega, do_plot = True)
	sim.calc_fft()
	sim.v_factor = 0.5e-2
	sim.plot_fft(do_renorm=False)

def fit_time_domain(driving_frq):
	fit = fit_rabi.fit_Rabi_time_domain (power=23, B_Gauss=63)
	fit.load_data()
	fit.get_data_driving_frq(driving_frq)

	params = lmfit.Parameters()
	params.add('Omega', value=1.e6, min = 0.e6, max=5e6, vary = True)
	params.add('t1', value=200e-9, min = 50e-9, vary = True)
	params.add('t2', value=100e-9, min = 50e-9, vary = True)
	params.add('I12', value=0, min = 0., max= 2., vary = False)
	params.add('PL_incr', value = 2., min = 1., max= 10., vary = True)
	params.add('init_pol', value = 0.9, min = 0.5, max = 1., vary = False)
	fit.do_fit(params)
	fit.plot_fit()

def fit_Rabi_all ():
	fit = fit_rabi.fit_Rabi_all (power=23)
	fit.load_data()
	fit.do_fit()

def fit_all_separately():
	fit = fit_rabi.fit_Rabi_time_domain(power=23, B_Gauss=63)
	fit.load_data()
	frq_array = fit.return_driving_frq()
	for f in frq_array:
		print " - - - - - Fitting: ", f, " MHz"
		fit.get_data_driving_frq(f)

		params = lmfit.Parameters()
		params.add('Omega', value=3.e6, min = 2.e6, max=5.e6, vary = True)
		params.add('t1', value=200e-9, min = 50e-9, vary = False)
		params.add('t2', value=100e-9, min = 50e-9, vary = True)
		params.add('I12', value=1., min = 0., max= 10., vary = False)
		params.add('PL_incr', value = .8, min = 0., max= 10., vary = True)
		params.add('init_pol', value = 0., min = 0., max = 1., vary = True)
		fit.do_fit(params)
		fit.plot_fit()

def compare_theory_exp_time_domain (power):
	Omega = 2.82e6
	t1 = 40e-9
	t2 = 150e-9
	I12 = 1.15
	PL_incr = 1.3
	init_pol = 0.9

	data = fit_rabi.fit_Rabi(power=power)
	data.load_data()
	for driving_frq in data.return_driving_frq():

		t_exp, PL_exp = data.get_data_driving_frq(driving_frq)

		t = np.linspace (0,1000,2000)*1e-9
		sim = V1.V1_SiC (B_Gauss = 63, T2_star = 1., t=t, verbose = False)
		sim.set_intensity_transitions(I_m32_m12=1, I_m12_p12=I12, I_p12_p32=1)
		sim.set_ODMR_pars (polariz_array = np.array([init_pol,1-init_pol,1-init_pol,init_pol]), ODMR_contrast_array = np.array([1.,PL_incr,PL_incr,1.]))
		sim.set_decay (t1 = t1, t2 = t2)
		PL_th = sim.rabi (f = driving_frq*1e6, Omega=Omega, do_plot = False)

		plt.figure(figsize = (15,5))
		plt.plot (t*1e6, PL_th, color='RoyalBlue', linewidth = 3)
		plt.plot (t_exp*1e6, PL_exp, 'o', color='crimson')
		plt.axis([0,1,0,2])
		plt.xlabel ('time [us]', fontsize=18)
		plt.ylabel ('normalized PL', fontsize=18)
		plt.title ('Power: '+str(power)+' dBm --- driving frq: '+str(driving_frq)+' MHz', fontsize=18)
		plt.show()

fit_all_separately()
#compare_theory_exp_time_domain(power=26)
#simulate_2D_Rabi()
#fit_time_domain(162.5)