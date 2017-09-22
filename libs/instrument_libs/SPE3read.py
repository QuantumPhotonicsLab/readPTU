# -*- coding: utf-8 -*-
"""
Created on Tue May 24 16:18:39 2016

@author: raphaelproux

Partly based on pyspec package by (c) Stuart B. Wilkins 2010
Adapted to use with SPE3 and specific use for a map
(A map is a series of frames with one region of interest only)

IMPORTANT NOTE: needs the xmltodict package to work. 
If not installed install it using PIP typing "pip install xmltodict" in a terminal

NOTE 2: in Spyder, use "%matplotlib TkAgg" in the iPython console to have independent window
"""
import pylab as py
import xmltodict
from tkFileDialog import askopenfilename
#from IPython import get_ipython
#ipython = get_ipython()
#
#ipython.magic("matplotlib qt")

class SPE3map:
    
    SPEVERSIONOFFSET = 1992  # offset for obtaining SPE version (float32)    
    XMLFOOTEROFFSETPOS = 678  # position in bytes giving the offset to the XML footer (UnsignedInteger64)
    DATAOFFSET = 4100  # offset to the binary data is fixed (4100 bytes) in SPE 2.X/3 file format
    
    def __init__(self, fname = None, fid = None):
        """Initialize class.
        Parameters:
           fname = Filename of SPE file
           fid = File ID of open stream

        This function initializes the class and, if either a filename or fid is
        provided opens the datafile and reads the contents"""
        
        self._fid = None
        self.fname = fname
        if fname is not None:
            self.openFile(fname)
        elif fid is not None:
            self._fid = fid

        if self._fid:
            self.readData()
            self._fid.close()
            
    def readData(self):
        """Read all the data into the class"""
        self._readSPEversion()
        try:
            assert(self.SPEversion >= 3)# or print 'This file is not a SPE 3.x file.'
        except:
            raise
#            print 'This program should be used only with SPE3.x files and higher. The file given has a {} version number.'.format(self.SPEversion)
        self._readXMLfooter()
        self._readWavelengths()
        self._readFramesInfo()
        self._readRegionSize()
        self._readExposureTime()
        self._readArray()
        
    def openFile(self, fname):
        """Open a SPE file"""
        self._fname = fname
        self._fid = open(fname, "rb")

    def _readAtNumpy(self, pos, size, ntype):
        """Return numpy array of 'size' elements of data type 'ntype', starting from position 'pos' bytes
        ATTENTION: 'pos' is in bytes while 'size' is in number of items of type ntype"""
        self._fid.seek(pos)
        return py.fromfile(self._fid, ntype, size)
        
    def _readSPEversion(self):
        """Determines SPE file version (always there in SPE 2.x or 3.0 files)"""
        self.SPEversion = self._readAtNumpy(self.SPEVERSIONOFFSET, 1, py.float32)[0]
        
    def _readXMLfooter(self):
        """Extracts the XML footer and puts it in _footerInfo as an ordered dictionnary (cf. xmltodict package)"""
        XMLfooterPos = self._readAtNumpy(self.XMLFOOTEROFFSETPOS, 1, py.uint64)[0]
        self._fid.seek(XMLfooterPos)
        self._footerInfo = xmltodict.parse(self._fid.read())

    def _readWavelengths(self):
        """Extracts the wavelength vector determined by spectrometer calibration"""
        wavelengthStr = self._footerInfo['SpeFormat']['Calibrations']['WavelengthMapping']['Wavelength']['#text']
        self.wavelength = py.array([float(w) for w in wavelengthStr.split(',')])
        
    def _readFramesInfo(self):
        """Extracts frames info from XML footer (number of frames, data type, frame size, frame stride)
        MUST BE CALLED AFTER _readXMLfooter()"""
        assert(self._footerInfo['SpeFormat']['DataFormat']['DataBlock']['@type'] == 'Frame')
        self.nbOfFrames = int(self._footerInfo['SpeFormat']['DataFormat']['DataBlock']['@count'])
        dataTypeName = self._footerInfo['SpeFormat']['DataFormat']['DataBlock']['@pixelFormat']
        possibleDataTypes = {'MonochromeUnsigned16': py.uint16,
                             'MonochromeUnsigned32': py.uint32,
                             'MonochromeFloat32': py.float32,
                             'MonochromeFloating32': py.float32}
        self.dataType = possibleDataTypes[dataTypeName]
        self.frameSize = int(self._footerInfo['SpeFormat']['DataFormat']['DataBlock']['@size'])
        self.frameStride = int(self._footerInfo['SpeFormat']['DataFormat']['DataBlock']['@stride'])
        
    def _readRegionSize(self):
        """Extracts width and height of the region of interest
        MUST BE CALLED AFTER _readXMLfooter()"""
        assert(self._footerInfo['SpeFormat']['DataFormat']['DataBlock']['DataBlock']['@type'] == 'Region')
        height = int(self._footerInfo['SpeFormat']['DataFormat']['DataBlock']['DataBlock']['@height'])
        width = int(self._footerInfo['SpeFormat']['DataFormat']['DataBlock']['DataBlock']['@width'])
        self.regionSize = (height,width)
        
    def _readExposureTime(self):
        """Extracts the camera exposure time
        MUST BE CALLED AFTER _readXMLfooter()"""
        self.exposureTime = float(self._footerInfo['SpeFormat']['DataHistories']['DataHistory']['Origin']['Experiment']['Devices']['Cameras']['Camera']['ShutterTiming']['ExposureTime']['#text'])
        
    def _readArray(self):
        """Reads the binary data contained in the file"""
        self.data = []
        for frameNb in xrange(self.nbOfFrames):
            frameData = self._readAtNumpy(self.DATAOFFSET + frameNb * self.frameStride, self.frameSize / self.dataType().nbytes, self.dataType)
            self.data.append(frameData.reshape(self.regionSize))
        self.data = py.array(self.data)
        
    def saveXMLinfo(self, filePath):
        """allows the user to save the XML footer to a file of his choice"""
        text_file = open(filePath, "w")
        text_file.write(xmltodict.unparse(self._footerInfo))
        text_file.close()


def range_to_edge(list):
    """Converts from a list of values to list of boundaries
    i.e. (1, 2, 3) becomes (0.5, 1.5, 2.5, 3.5)
    Used for pcolorfast plotting"""
    edges = py.array([])
    edges = py.append(edges,list[0] - (list[1]-list[0])/2)
    i=0
    while i < len(list):
        edges = py.append(edges,list[i] + (list[i]-edges[i]))
        i+=1
    return edges

if __name__ == "__main__":
    
#    data = SPE3map("/Users/raphaelproux/Desktop/PL-map/2016-05-18_1200GR_int_1s_center_wvl_935nm_Bias_-0,7V_to_0,3V_101steps_Exc_550mV_1E4_off_SIL_Sheffield_PL_Map.spe")
#    data = SPE3map("/Users/raphaelproux/Desktop/PL-map/pol290_P=300mV_wlen=969_38 2016 May 13 19_39_35.spe")
    
    filename = askopenfilename()
    data = SPE3map(filename)
    dataArray = py.array([dataRow[0] for dataRow in data.data])

    voltageRange = (-0.3,0.7)  # real voltage range (parameter of measurement)
    voltage = py.linspace(voltageRange[0], voltageRange[1], data.nbOfFrames)
    
    voltagePlot = range_to_edge(voltage) 
    voltagePlotRange = (voltagePlot[0],voltagePlot[-1])
    
    wavelengthPlot = range_to_edge(data.wavelength) 
    wavelengthPlotRange = (wavelengthPlot[0],wavelengthPlot[-1])
    
    colorPlotRange = (00,500)
    
    fig = py.figure(figsize=(8.*4./3.*0.8,8.*0.8), dpi=100)
    ax = fig.add_subplot(111)
    pcf = ax.pcolormesh(voltagePlot,
                        wavelengthPlot,
                        dataArray.transpose(),
                        vmin=colorPlotRange[0],
                        vmax=colorPlotRange[1])
    ax.set_xlim(voltagePlotRange)
    ax.set_ylim(wavelengthPlotRange)
    ax.set_xlabel('Voltage (V)')
    ax.set_ylabel('Wavelength (nm)')
#
#    self.ax.grid(True)
#    self.fig.subplots_adjust(left=0.10,right=0.92) # before colorbar
#    cb=self.fig.colorbar(self.pcf,fraction=0.05,pad=0.015)
#    cb.set_label("Counts in %0.1f s" % self.spe.time)
#    self.ax.set_title(os.path.basename(self.filename))
#    self.ax.get_xaxis().set_label_text('Gate  Voltage (V)')
#    self.ax.get_yaxis().set_label_text('Wavelength (nm)')
#    self.ax.get_xaxis().get_major_formatter().set_useOffset(False)
#    self.ax.get_yaxis().get_major_formatter().set_useOffset(False)
#    self.build()