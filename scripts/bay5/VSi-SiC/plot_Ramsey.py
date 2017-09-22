
import numpy as np
import pylab as plt
import os
import cPickle

def get_data(folder, fname):
	filename = os.path.join (folder, fname)
	with open(filename,'rb') as f:
		d=cPickle.load(f)
		n = len(d['spin_state'])
		PL1 = d['spin_state'][:n/2]
		PL2 = d['spin_state'][n/2:]
		PL = PL1-PL2
		PL = PL - np.mean(PL)
		PL = PL/max(PL)
		return d['tau'].transpose(), PL

folder = 'D:/Research/__current_work/__ongoing_research/SiC_spectroscopy/ramsey/'
f1 = '-10dBm_861.59nm_730nm_875LP_no900SP_3600ns_2.820e9_72nspihalf.pys'
f2 = '2016-10-19_CREE5E175.5K_0.5mW_875LP_900SP_861nm_730nm-19.0dBm1000.0sLaser500.0nsAOM996.0nsInit96.0nsRead_165.00MHz.pys'

t1, PL1 = get_data (folder=folder, fname=f1)
t2, PL2 = get_data (folder=folder, fname=f2)

fig = plt.figure(figsize = (6,3))
plt.plot (t1*1e-3, PL1, 'o', color = 'crimson', markersize=4)
plt.plot (t1*1e-3, PL1, color = 'crimson', linewidth = 1)
plt.plot (t2*1e-3, PL2, '^', color = 'RoyalBlue', markersize=6)
plt.plot (t2*1e-3, PL2, color = 'royalblue', linewidth = 1)
plt.axis('tight')
plt.xlabel ('free evolution time [us]', fontsize=16)
fig.savefig ('D:/Research/ramsey.svg')
plt.show()
