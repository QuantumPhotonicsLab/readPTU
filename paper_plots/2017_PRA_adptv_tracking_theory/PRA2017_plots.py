
import numpy as np
import random
from matplotlib import rc, cm
import os, sys
import h5py
import logging, time

from matplotlib import pyplot as plt
from scipy import signal
import lmfit
matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20)

from analysis_simulations.libs.adaptive_sensing import adaptive_tracking as track_lib
from analysis_simulations.libs.adaptive_sensing import adaptive_tracking_overhead as track_libOH
from analysis_simulations.libs.adaptive_sensing import tracking_sequence as trackSeq

reload (track_lib)
reload (track_libOH)
reload (trackSeq)

t0 = 20e-9
G = 5
F = 3

def sens_units (N, G, F):
	T = G*(2**N-1) + F*(2**N-N-1)
	return T

def oh_units(N, G, F):
	R = G*N+0.5*F*N*(N-1)
	return R


def theory_epsilon (df, tOH, fid=1, sweep_par='df', track = False, prop_factor = 1.):
	c = 3**0.5/(2*np.pi)
	t0 =20e-9
	G = 5
	F = 3
	if track:
		c1 = G
	else:
		c1 = (2*(G+F))

	if (sweep_par == 'df'):
		sweep_array = df
	else:
		sweep_array = tOH

	opt_N = np.zeros(len(sweep_array))
	opt_epsilon = np.zeros(len(sweep_array))
	opt_e1 = np.zeros(len(sweep_array))

	eta = 2
	ind = 0
	for i in sweep_array:
		#print "df = ", i*1e-6
		N = 13
		e1 = 0
		e2 = 1
		while prop_factor*e1<=eta*e2:

			if track:
				R = G
			else:
				R = (G*N + 0.5*F*N*(N-1))

			dt = 2**(N-1)*t0
			csi = (1+(1-fid)*N)**0.5
			e1 = c*csi*(1./(c1**0.5*dt))

			if sweep_par=='df':
				e2 = i*(c1*dt+R*tOH)**0.5
			else:
				e2 = df*(c1*dt+R*i)**0.5

			N=N-0.01
		N = N+0.01
		opt_N[ind] = N
		if not(track):
			R = G*N + 0.5*F*N*(N-1)

		if sweep_par == 'df':
			opt_epsilon[ind] = e2
		else:
			opt_epsilon[ind] = e2
		opt_e1[ind] = e1
		ind += 1

	return opt_N, opt_epsilon




def plot_eps_vs_overhead():
	#plot epsilon vs overhead for track and no-track (dB = 2MHz)

	fig = plt.figure(figsize = (12, 8))
	folder = 'D:/Research/WorkData/adptv_tracking_sim/vsOH'

	ax1 = fig.add_subplot (2,1,1)
	#Tracking
	stamp_tr = 'analysis_track'
	a = trackSeq.Analyze(folder = folder)
	a.load_analysis_results (folder = folder, stamp = stamp_tr)
	OH = a.sweep_par
	y_tr = a.waveform_estimation_rmse
	s_tr = a.waveform_estimation_rmse_std
	OH_th = np.linspace (OH[0], OH[-1], 1000)
	n, y_tr_th = theory_epsilon (df=2e6, tOH=OH_th, sweep_par = 'OH', track = True, fid=1.)
	ax1.plot (OH_th*1e6, y_tr_th*1e-3, color = 'crimson', linewidth=2, label='track')
	ax1.errorbar (OH*1e6, y_tr*1e-3, yerr=s_tr*1e-3, fmt='^', markersize=7, color='crimson')
	
	#No-Tracking
	stamp_notr = 'analysis_noTrack'
	a = trackSeq.Analyze(folder = folder)
	a.load_analysis_results (folder = folder, stamp = stamp_notr)
	OH = a.sweep_par
	y_notr = a.waveform_estimation_rmse
	s_notr = a.waveform_estimation_rmse_std
	OH_th = np.linspace (OH[0], OH[-1], 1000)
	n, y_notr_th = theory_epsilon (df=2e6, tOH=OH_th, sweep_par = 'OH', track = False, fid=1.)
	ax1.plot (OH_th*1e6, y_notr_th*1e-3, color = 'RoyalBlue', linewidth=2, label='no-track')
	ax1.set_xlim([-10, 310])
	ax1.errorbar (OH*1e6, y_notr*1e-3, yerr=s_notr*1e-3, fmt='v', markersize=7, color='RoyalBlue')
	y1_max=250
	y1_min=0
	ax1.set_ylim([y1_min,y1_max])
	ax1.set_xlabel ('overhead [us]', fontsize=18)
	ax1.set_ylabel ('epsilon [kHz]', fontsize=18)
	ax1.legend(loc=2, fontsize=16)
	ax2 = ax1.twinx()
	gamma = 28
	ax2.set_ylim([y1_min/gamma, y1_max/gamma])

	
	x = OH*1e6
	plt.subplot (2,1,2)
	s_f = (y_notr/y_tr)*((s_tr/y_tr)**2+(s_notr/y_notr)**2)**0.5
	plt.plot (OH_th*1e6, y_notr_th/y_tr_th, color='lightgray', linewidth=2)
	plt.errorbar (x[:], y_notr[:]/y_tr[:], yerr=s_f[:], fmt='o', markersize=7, color = 'k')
	plt.axis ([-10.,310, 0.,5.])
	plt.xlabel ('overhead [us]', fontsize=18)
	plt.ylabel ('Ratio error [no-track/track]', fontsize=18)

	fig.savefig ('D:/Research/WorkData/adptv_tracking_sim/eps_vs_OH_dB=2MHz.svg')
	
	plt.show()

def plot_vs_fidelity():
	#still to be finished
	stamp = '20161207_204705_analysis_adptvTr_vs_fid0_fid1=100'
	a = trackSeq.Analyze(folder = folder)
	a.load_analysis_results (folder = folder, stamp = stamp)
	fid0 = a.sweep_par
	y = a.waveform_estimation_rmse
	s = a.waveform_estimation_rmse_std
	#fid0_thr = np.linspace (fid0[0], fid0[-1], 1000)
	#n, y_th = theory_epsilon (df=2e6, tOH=OH_th, sweep_par = 'OH', track = True, fid=1.)
	plt.plot (fid0, y*1e-3, color = 'crimson', linewidth=2)
	plt.errorbar (fid0, y*1e-3, yerr=s*1e-3, fmt='o', markersize=7, color='crimson')
	plt.xlabel ('read-out fidelity', fontsize=18)
	plt.ylabel ('error [kHz]', fontsize=18)
	plt.show()

def plot_eps_vs_dB_fid088 ():
	folder = 'D:/Research/WorkData/adptv_tracking_sim/fid_088'
	plt.figure(figsize=(15,6))

	stamp = '20161208_113857_analysis_noTr_vs_dB_fid=88_OH100us'
	a = trackSeq.Analyze(folder = folder)
	a.load_analysis_results (folder = folder, stamp = stamp)
	dB = a.sweep_par
	y_100us_notr = a.waveform_estimation_rmse
	s_100us_notr = a.waveform_estimation_rmse_std
	dB_th = np.linspace (dB[0], dB[-1], 100)
	n, y_th_100us_notr = theory_epsilon (df=dB_th, tOH=100e-6, sweep_par = 'df', track = False, fid=1.)
	plt.plot (dB_th*1e-6, y_th_100us_notr*1e-3, '--',color = 'RoyalBlue', linewidth=1)
	plt.errorbar (dB*1e-6, y_100us_notr*1e-3, yerr=s_100us_notr*1e-3, fmt='o', markersize=7, color='RoyalBlue')

	stamp = '20161208_114212_analysis_adptvTr_vs_dB_fid=88_OH100us'
	a = trackSeq.Analyze(folder = folder)
	a.load_analysis_results (folder = folder, stamp = stamp)
	dB = a.sweep_par
	y_100us_tr = a.waveform_estimation_rmse
	s_100us_tr = a.waveform_estimation_rmse_std
	dB_th = np.linspace (dB[0], dB[-1], 100)
	n, y_th_100us_tr = theory_epsilon (df=dB_th, tOH=100e-6, sweep_par = 'df', track = True, fid=1.)
	plt.plot (dB_th*1e-6, y_th_100us_tr*1e-3, '--',color = 'crimson', linewidth=1, label='track')
	plt.errorbar (dB*1e-6, y_100us_tr*1e-3, yerr=s_100us_tr*1e-3, fmt='o', markersize=7, color='crimson')
	plt.xlabel ('Delta_f [MHz*sqrt(Hz)]', fontsize=18)
	plt.ylabel ('error [kHz]', fontsize=18)
	plt.axis([0, 20, 0, 800])
	plt.show()



def plot_eps_vs_dB(stamp_track, stamp_notrack, 
			folder = 'D:/Research/WorkData/adptv_tracking_sim/',
			ylim =[], xlim =[], ylim_eta =[],
			OH = 100e-6, fid = 1.,
			do_save = False, filename = 'plot.png'):

	fig = plt.figure(figsize = (12, 8))

	ax1 = fig.add_subplot (2,1,1)
	#Tracking
	a = trackSeq.Analyze(folder = folder)
	a.load_analysis_results (folder = folder, stamp = stamp_track)
	dB = a.sweep_par
	y1 = a.waveform_estimation_rmse
	s1 = a.waveform_estimation_rmse_std
	dB_th = np.linspace (dB[0], dB[-1], 1000)
	n, y_tr = theory_epsilon (df=dB_th, tOH=OH, sweep_par = 'df', track = True, fid=fid)
	ax1.plot (dB_th*1e-6, y_tr*1e-3, color = 'crimson', linewidth=2, label='track')
	ax1.errorbar (dB*1e-6, y1*1e-3, yerr=s1*1e-3, fmt='^', markersize=7, color='crimson')
	
	#No-Tracking
	a = trackSeq.Analyze(folder = folder)
	a.load_analysis_results (folder = folder, stamp = stamp_notrack)
	dB = a.sweep_par
	y2 = a.waveform_estimation_rmse
	s2 = a.waveform_estimation_rmse_std
	G = 5
	dB_th = np.linspace (dB[0], dB[-1], 1000)
	n, y_notr = theory_epsilon (df=dB_th, tOH=OH, sweep_par = 'df', track = False, fid=fid)
	ax1.plot (dB_th*1e-6, y_notr*1e-3, color = 'RoyalBlue', linewidth=2, label='no-track')
	ax1.errorbar (dB*1e-6, y2*1e-3, yerr=s2*1e-3, fmt='v', markersize=7, color='RoyalBlue')
	ax1.set_xlim (xlim)
	ax1.set_ylim(ylim)
	ax1.set_xlabel ('Delta_f [MHz*sqrt(Hz)]', fontsize=18)
	ax1.set_ylabel ('epsilon [kHz]', fontsize=18)
	ax1.legend(loc=2, fontsize=16)

	ax2 = ax1.twinx()
	ax2.set_ylim ([ylim[0]/28., ylim[1]/28.])

	ax3 = fig.add_subplot (2,1,2)
	eta = y2/y1
	s_eta = (y2/y1)*((s1/y1)**2+(s2/y2)**2)**0.5
	ax3.plot (dB_th*1e-6, y_notr/y_tr, color='lightgray', linewidth=4)
	ax3.errorbar (dB*1e-6, eta, yerr=s_eta, fmt='o', markersize=7, color = 'k')
	ax3.set_xlim(xlim)
	ax3.set_ylim (ylim_eta)
	ax3.set_xlabel ('Delta_f [MHz*sqrt(Hz)]', fontsize=18)
	ax3.set_ylabel ('Ratio error [no-track/track]', fontsize=18)

	if do_save:
		fig.savefig (filename, dpi = 400)
	plt.show()


def plot_notracking_fig1(do_save = False):

	folder = 'D:/Research/WorkData/adptv_tracking_sim/OH_0.01us/'
	fig = plt.figure(figsize = (8, 3))

	ax1 = fig.add_subplot (1,1,1)
	a = trackSeq.Analyze(folder = folder)
	a.load_analysis_results (folder = folder, 
				stamp = '20170110_115007_analysis_noTracking')
	dB = a.sweep_par
	y1 = a.waveform_estimation_rmse
	s1 = a.waveform_estimation_rmse_std
	dB_th = np.linspace (dB[0], dB[-1], 1000)
	n, y_tr = theory_epsilon (df=dB_th, tOH=10e-9, sweep_par = 'df', track = False, fid=1.)
	ax1.plot (dB_th*1e-6, y_tr*1e-3, color = 'RoyalBlue', linewidth=2)
	ax1.errorbar (dB*1e-6, y1*1e-3, yerr=s1*1e-3, fmt='^', markersize=8, color='RoyalBlue')
	for xlabel_i in ax1.get_xticklabels():
		xlabel_i.set_fontsize(20) 
	for xlabel_i in ax1.get_yticklabels():
		xlabel_i.set_fontsize(20) 

	ax1.set_xlim ([-1, 51])
	ax1.set_ylim([0, 140])
	ax1.set_xlabel ('Delta_f [MHz*sqrt(Hz)]', fontsize=18)
	ax1.set_ylabel ('epsilon [kHz]', fontsize=18)
	if do_save:
		fig.savefig (folder+'plot.svg', dpi = 400)
	plt.show()


def fidelity_plots (do_save=False):

	OH = 100.e-6
	folder = 'D:/Research/WorkData/adptv_tracking_sim/lower_fidelity/'
	stamp_track75 = '20170118_134158_analysis_adptvTracking_fid75'
	stamp_track88 = 'analysis_adptvTracking_fid88'
	stamp_notrack75 = '20170118_070732_analysis_noTracking_fid75'
	stamp_notrack88 = '20170118_070453_analysis_noTracking_fid88'

	fig = plt.figure(figsize = (8, 10))

	ax1 = fig.add_subplot (2,1,1)
	#Tracking-fid75	
	a = trackSeq.Analyze(folder = folder)
	a.load_analysis_results (folder = folder, stamp = stamp_track75)
	dB = a.sweep_par
	y_tr75 = a.waveform_estimation_rmse
	s_tr75 = a.waveform_estimation_rmse_std
	ind = np.argsort(dB)
	dB = dB[ind]
	y_tr75 = y_tr75[ind]
	s_tr75 = s_tr75[ind]
	dB0 = np.linspace (min(dB), max(dB), 100)
	n, y0_tr75 = theory_epsilon (df=dB0, tOH=OH, sweep_par = 'df', track = True, fid=1, prop_factor = 2000)
	ax1.plot (dB0*1e-6, y0_tr75*1e-3, color = [250/255., 122/255., 109/255., 1.], linewidth=2, label='fid=0.75')
	ax1.errorbar (dB*1e-6, y_tr75*1e-3, yerr=s_tr75*1e-3, fmt='o', markersize=7, color=[250/255.,122/255., 109/255., 1.])
	
	
	#Tracking-fid88
	a = trackSeq.Analyze(folder = folder)
	a.load_analysis_results (folder = folder, stamp = stamp_track88)
	dB = a.sweep_par
	y_tr88 = a.waveform_estimation_rmse
	s_tr88 = a.waveform_estimation_rmse_std
	ind = np.argsort(dB)
	dB = dB[ind]
	y_tr88 = y_tr88[ind]
	s_tr88 = s_tr88[ind]

	n, y0_tr88 = theory_epsilon (df=dB0, tOH=OH, sweep_par = 'df', track = True, fid=1., prop_factor = 300)
	ax1.plot (dB0*1e-6, y0_tr88*1e-3, color = [200/255., 53/255., 35/255., 1.], linewidth=1, label='fid=0.88')
	ax1.errorbar (dB*1e-6, y_tr88*1e-3, yerr=s_tr88*1e-3, fmt='s', markersize=7, color=[200/255., 53/255., 35/255., 1.])
	#Tracking-fid100
	n, y0_tr100 = theory_epsilon (df=dB0, tOH=OH, sweep_par = 'df', track = True, fid=1.)
	ax1.plot (dB0*1e-6, y0_tr100*1e-3, '--', color = 'crimson', linewidth=1, label='fid = 1.00')
	

	#NO-tracking-fid75	
	a = trackSeq.Analyze(folder = folder)
	a.load_analysis_results (folder = folder, stamp = stamp_notrack75)
	dB = a.sweep_par
	y_notr75 = a.waveform_estimation_rmse
	s_notr75 = a.waveform_estimation_rmse_std
	ind = np.argsort(dB)
	dB = dB[ind]
	y_notr75 = y_notr75[ind]
	s_notr75 = s_notr75[ind]

	dB0 = np.linspace (min(dB), max(dB), 100)
	n, y0_notr75 = theory_epsilon (df=dB0, tOH=OH, sweep_par = 'df', track = False, fid=1., prop_factor=3)
	ax1.plot (dB0*1e-6, y0_notr75*1e-3, color = [129/255., 175/255., 254/255., 1.], linewidth=1, label='fid=0.75')
	ax1.errorbar (dB*1e-6, y_notr75*1e-3, yerr=s_notr75*1e-3, fmt='o', markersize=7, color=[129/255., 175/255., 254/255., 1.])
	#NO-racking-fid88
	a = trackSeq.Analyze(folder = folder)
	a.load_analysis_results (folder = folder, stamp = stamp_notrack88)
	dB = a.sweep_par
	y_notr88 = a.waveform_estimation_rmse
	s_notr88 = a.waveform_estimation_rmse_std
	ind = np.argsort(dB)
	dB = dB[ind]
	y_notr88 = y_notr88[ind]
	s_notr88 = s_notr88[ind]

	n, y0_notr88 = theory_epsilon (df=dB0, tOH=OH, sweep_par = 'df', track = False, fid=1., prop_factor=1.5)
	ax1.plot (dB0*1e-6, y0_notr88*1e-3, color = [56/255., 118/255., 224/255., 1.], linewidth=1, label='fid=0.88')
	ax1.errorbar (dB*1e-6, y_notr88*1e-3, yerr=s_notr88*1e-3, fmt='s', markersize=7, color=[56/255., 118/255., 224/255., 1.])
	#Tracking-fid100
	n, y0_notr100 = theory_epsilon (df=dB0, tOH=OH, sweep_par = 'df', track = False, fid=1.)
	ax1.plot (dB0*1e-6, y0_notr100*1e-3, '--', color = 'RoyalBlue', linewidth=2, label='fid=1.00')

	xlim = [-1, 21]
	ax1.set_xlim (xlim)
	ax1.set_ylim ([0, 1200])
	ax1.set_xlabel ('Delta_f [MHz*sqrt(Hz)]', fontsize=18)
	ax1.set_ylabel ('epsilon [kHz]', fontsize=18)
	ax1.legend(loc=2, fontsize=16)

	#ax2 = ax1.twinx()
	#ax2.set_ylim ([ylim[0]/28., ylim[1]/28.])

	
	ax3 = fig.add_subplot (2,1,2)
	eta75 = y_notr75/y_tr75
	eta0_75 = y0_notr75/y0_tr75
	s_75 = eta75*((s_tr75/y_tr75)**2+(s_notr75/y_notr75)**2)**0.5
	ax3.plot (dB0*1e-6, eta0_75, color='lightgray', linewidth=2, label = 'fid=0.75')
	ax3.errorbar (dB*1e-6, eta75, yerr=s_75, fmt='o', markersize=7, color='lightgray')
	eta88 = y_notr88/y_tr88
	eta0_88 = y0_notr88/y0_tr88
	s_88 = eta88*((s_tr88/y_tr88)**2+(s_notr88/y_notr88)**2)**0.5
	ax3.plot (dB0*1e-6, eta0_88, color='darkgray', linewidth=2, label = 'fid=0.88')
	ax3.errorbar (dB*1e-6, eta88, yerr=s_88, fmt='s', markersize=7, color='darkgray')
	eta100 = y0_notr100/y0_tr100
	ax3.plot (dB0*1e-6, eta100, '--', color='k', linewidth=2)
	xlim = [-1, 21]
	ax3.set_xlim (xlim)
	ylim = [1, 4]
	ax3.set_ylim (ylim)

	#ax3.errorbar (dB*1e-6, eta, yerr=s_eta, fmt='o', markersize=7, color = 'k')
	#ax3.set_xlim(xlim)
	#ax3.set_ylim (ylim_eta)
	ax3.set_xlabel ('Delta_f [MHz*sqrt(Hz)]', fontsize=18)
	ax3.set_ylabel ('Ratio error [no-track/track]', fontsize=18)
	ax3.legend(loc=1, fontsize=16)


	if do_save:
		fig.savefig ('D:/low_fid_sims.svg', dpi = 400)
	
	plt.show()


#figures for fid=1
#plot_notracking_fig1(do_save = True)
'''
#OH = 100us
plot_eps_vs_dB (stamp_track = '20170106_125208', stamp_notrack = '20170106_140251',
					xlim = [0, 26.], ylim = [0, 1100], ylim_eta = [0, 5],
					OH = 100.e-6, fid = 1.,
					do_save = True, filename = 'eps_vs_dB_OH=100us.svg')

#OH = 10ns
plot_eps_vs_dB (stamp_track = '20161205_091259_analysis_adptvTrack_OH=10ns_G5F3flex',
					stamp_notrack = '20161205_094610_analysis_noTrack_OH=10ns_G5F3',
					xlim = [0, 51.], ylim = [0, 150], ylim_eta = [0, 2.],
					folder = 'D:/Research/WorkData/adptv_tracking_sim/OH_0.01us/_old_sims/',
					OH = 10e-9, fid =1.,
					do_save = True, filename = 'eps_vs_dB_OH=10ns.svg')
'''
#dB = 2MHz
#plot_eps_vs_overhead()
fidelity_plots(do_save=True)

#figures and tests for fid=0.88
#plot_eps_vs_dB_fid088()
#plot_vs_fidelity()


