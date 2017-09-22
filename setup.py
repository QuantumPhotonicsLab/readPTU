
import sys
import os
import numpy as np
import h5py
import socket


#hostpc = os.getenv ('computername')
hostpc = socket.gethostname()

settings = {}

settings ['DaleMac.local'] = {
    'description': 'Dale Mac',
    'folder': "/Users/dalescerri/Documents/GitHub/Spins/",
    }
settings['BAY3-HP'] = {
	'description': 'Bay3 pc in DB2.17',
	'folder': "C:/Research/",
	}
settings['cristian-mint'] = {
	'description': 'Cristian - old TU Delft laptop',
	'folder': "/home/cristian/Work/bonato-lab/",
	}
settings['cristian-PC'] = {
	'description': 'Cristian - Thinkpad',
	'folder': "C:/Users/cristian/Research/bonato-lab/",	
}	
settings['HWPC0526-EPS'] = {
	'description': 'Cristian - office pc',
	'folder': "H:/Research/bonato-lab/",
	}

print 'Loaded settings for: ', settings[hostpc]['description']

folder = settings[hostpc]['folder']
sys.path.append (folder)

sys.path.append (folder+ '/analysis_simulations/')
os.chdir (folder+ '/analysis_simulations/')

global root_folder 
root_folder = folder

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
mpl.rc('xtick', labelsize=18) 
mpl.rc('ytick', labelsize=18)

try:
	_viridis_data = np.load (folder+'/viridis_cmp.npy')
	viridis = ListedColormap(_viridis_data, name='viridis')
	plt.register_cmap(name='viridis', cmap=viridis)
except:
	pass
