import numpy as np
import pylab as plt
import os


class txtData ():

	def __init__ (self, data_folder):
		self.folder = data_folder

	def skip_lines (self, value):
		self.skip = value

	def load_data (self, file_name, nr_chans):
		fName = os.path.join (self.folder, file_name)
		lines = np.loadtxt(fName, comments="#", dtype='str', delimiter="\n", unpack=False)
		for i in np.arange(self.skip):
			a =  lines[i].split()
			a[0] = a[0][:-1]
			if (isinstance (a[1], basestring) and len(a)>2):
				a [0] = a[0]+'_'+a[1]
				a[1] = a[2]
			try:
				setattr (self, a[0], float(a[1]))
			except:
				pass	

		for i in np.arange(nr_chans):
			setattr (self, 'x'+str(i), np.zeros(len(lines)-self.skip))
			setattr (self, 'y'+str(i), np.zeros(len(lines)-self.skip))

		for i in np.arange (nr_chans):
			for j in np.arange(len(lines)-self.skip):
				a = lines[j+self.skip].split()
				x = getattr (self, 'x'+str(i))
				self.x[j] = float(a[0])
				self.y[j] = float(a[1])

	def plot_data (self):

		plt.figure(figsize=(15, 8))
		plt.plot (self.x*1e9, self.y)
		plt.show()



txt = txtData ('C:/Research/Data/20170531_AWG_M3202A_tests/')
txt.skip_lines (18)
txt.load_data ('wfA.txt')
txt.plot_data()

