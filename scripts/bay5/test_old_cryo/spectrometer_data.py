import numpy as np
import pylab as plt
import os


def load_data (data_folder, angle):
	f = 'D:/Research/__current_work/__ongoing_research/SiC_spectroscopy/resonant_absorption_polarResolved/'
	folder = f+data_folder
	fName = '2016-10-18_CREE_2016-10-18_CREE5E17_858nm_600uW_no_filter_V1 Prime Calibration__'+str(angle)+'_degree_'
	#fName = '2016-10-18_CREE_2016-10-18_CREE5E17_858nm_600uW_870LP_V1 Calibration__'+str(angle+0.)+'_degree_'
	#fName = '2016-10-18_CREE_2016-10-18_CREE5E17_858nm_600uW_870LP_V1_Prime_Calibration_2__'+str(angle+0.0)+'_degree_'
	filename = os.path.join (folder, fName+'.txt')
	lines = np.loadtxt(filename, comments="#", delimiter=" ", unpack=False)
	data_array = np.asarray(lines)
	wavelength = data_array [:, 0]
	intensity = data_array [:,1]
	laser_int = max(intensity[400:500])
	#laser_int = np.mean(intensity[0:200])
	return wavelength, intensity/laser_int

angle_res = 5
n_datasets = 180/angle_res
n_bins = 1340
A = np.zeros ((n_datasets, n_bins))
for i in np.arange (n_datasets):
	x,y = load_data (data_folder = 'V1 Prime Calibration/First_Measurement/', angle=i*angle_res)
	A[i, :] = y

plt.figure(figsize = (20,7))
plt.pcolor (A)
plt.show()