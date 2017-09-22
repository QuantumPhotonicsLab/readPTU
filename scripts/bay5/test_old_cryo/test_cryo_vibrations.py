import pylab as py
import matplotlib.colors
import json

from analysis_simulations.libs.drivers import SPE3read as spe3


#
#

fName = '/home/carlo/Scrivania/Heriot-Watt/Cryostat_spectroscopy_measurements/2017 May 17 17_18_51.spe'
data = spe3.SPE3map(fName)

wavelength = data.wavelength
wavelengthRangeMeas = (wavelength[0], wavelength[-1])
nbOfFrames = data.nbOfFrames
counts = py.array([dataRow[0] for dataRow in data.data])
exposureTime = data.exposureTime  # exposure time in ms

print wavelength
print counts
print len(wavelength)
print exposureTime

for j in range(0, len(counts)):
	for i in range(0,len(wavelength)):
		if ((wavelength[i] > 927.3)and(wavelength[i] < 928.3)):
			print i,  ' ', wavelength[i]

for j in range(0,10):
	py.figure(j)
	py.plot(wavelength[328:354], counts[j][328:354])





