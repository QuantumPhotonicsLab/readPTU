'''
Author: Z.X. Koong 
Date : 31/12/2020

Ver 1 : Open a ptu file, analyze the data, do post-selection, save the figure (both g2 and timetrace) and data (in npz).
Ver 2 : Add options to analyze pulsed g2. Fix some minor issues with figure popping out when it is not needed. Save the previous parameters in a json cache file so one do not have to re-adjust the zero-delay for each measurements.
'''

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QWidget, QFileDialog
from PyQt5.QtCore import QTimer
import matplotlib.pyplot as plt
try:
    from readPTU import PTUmeasurement,PTUfile
except Exception as ep:
    print (ep)
    print ('Fix this error before proceeding!')
import numpy as np
from ui_analyse_ptu import Ui_Dialog
import sys, os

import matplotlib
matplotlib.use('Qt5Agg')

matplotlib.rcParams.update({'font.size': 30})
matplotlib.rcParams['lines.linewidth'] = 4
matplotlib.rcParams['font.family']='arial'
matplotlib.rcParams['mathtext.fontset'] = 'stix'


def construct_postselect_vector(timetrace_y, timetrace_recnum, threshold, constraint = "above",above = True):
    """Constructs a postselection vector based on a threshold condition, selecting items above or below as specified by user.
    
    Args:
        timetrace_y (numpy array): number of photons in the time bin.
        timetrace_recnum (numpy array): array of record numbers corresponding to timebin start (for each time bin of the time trace)
        threshold (number): number of photons threshold used for selection.
        above (string, optional): Default : "above". Options ="above", "below" and "between".
        "above"  : return time ranges where timetrace_y > threshold
        "below"  : return time ranges where timetrace_y < threshold
        "between": return time ranges where timetrace_y > threshold[0] and timetrace_y <= threshold[1]
        Anything apart from this would result in no post-selection.
        # Descriptions below are not valid anymore: 
        # above (bool, optional): if True, will return time ranges where timetrace_y > threshold. If False, will return time ranges where timetrace_y < threshold.
    
    Returns:
        (list, list): first: List of 2-item lists being [start_index, stop_index] array indices in timetrace_y for each post-selected range of points.
                      second: corresponding record numbers allowing for direct use as post-selection array in the calculate_g2() function.
    """
    if constraint == 'above':
        
        if not np.isscalar(threshold):
            print ('Threshold should be a scalar. Threshold is set to %s instead.'%threshold[0])
            threshold = threshold[0]
        select = timetrace_y > threshold
        if timetrace_y[0] > threshold:
            add_first_time = True
        else:
            add_first_time = False
    elif constraint == 'below':
        if not np.isscalar(threshold):
            print ('Threshold should be a scalar. Threshold is set to %s instead.'%threshold[0])
            threshold = threshold[0]
        select = timetrace_y < threshold
        if timetrace_y[0] < threshold:
            add_first_time = True
        else:
            add_first_time = False
    elif constraint == 'between':
        if np.isscalar(threshold):
            print ('Threshold should be an array. Threshold is set to %s instead.'%0)
            threshold = [0,threshold[0]]
        select = (timetrace_y > threshold[0]) == (timetrace_y <= threshold[-1])
        
        if timetrace_y[0] > threshold[0] or timetrace_y[0] <= threshold[-1]:
            add_first_time = True
        else:
            add_first_time = False
    else:
        print ('Unknown definition for constraint : %s. No post-selection'%constraint)
        threshold = 0
        select = timetrace_y > threshold
        if timetrace_y[0] > threshold:
            add_first_time = True
        else:
            add_first_time = False
    nselect = ~select
    post_selec_ranges = ((select[:-1] & nselect[1:]) | (nselect[:-1] & select[1:])).nonzero()[0]
    if add_first_time:
        post_selec_ranges = np.insert(post_selec_ranges, 0, 0).astype(int)
    if len(post_selec_ranges) % 2 != 0:  # odd length means we need to add the end time
        post_selec_ranges = np.append(post_selec_ranges, len(timetrace_y) - 2).astype(int)  # - 2 because 1 will be added below

    # # take into account the end of the interval (otherwise stops 1 time bin before)
    # for i, value in enumerate(post_selec_ranges[1::2]):
    #     post_selec_ranges[2*i - 1] += 1

    # # the previous step will introduce [0,11], [11,15] kind of situation (where ranges are contiguous). We remove these intermediate duplicates.
    # double_index = (post_selec_ranges[:-1] != post_selec_ranges[1:])
    # double_index = (np.append(double_index, True) & np.insert(double_index, 0, True))
    # post_selec_ranges = post_selec_ranges[double_index]

    post_selec_ranges = post_selec_ranges.reshape((-1, 2))

    recnum_post_selec_ranges = [[timetrace_recnum[post_selec_range[0]], timetrace_recnum[post_selec_range[1]]] for post_selec_range in post_selec_ranges]

    return post_selec_ranges, recnum_post_selec_ranges


class AnalysePTU_window(QWidget):
    def __init__(self,parent=None):
        super().__init__(parent=parent)
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)

        # flags
        self.calculatetimetrace = self.ui.istimetrace.isChecked()
        self.calculateg2 = self.ui.isG2.isChecked()
        self.saveData =  self.ui.isSavedData.isChecked()
        self.saveFigure = self.ui.isSaveFigure.isChecked()
        self.f_mapOpen  = False

        # initialization
        self.timetracex = None
        self.countsch0 = None
        self.countsch1 = None
        self.delay = None
        self.coincidences = None
        self.filename = ''        

        # connect to the buttons
        self.ui.openfile.clicked.connect(self.selectFile)
        self.ui.savebutton.clicked.connect(self.saveFile)
        self.ui.CalculateG2.clicked.connect(self.analysePTU)

        # connect to the toggles
        self.ui.istimetrace.toggled.connect(self.timetraceToggled)
        self.ui.isG2.toggled.connect(self.g2Toggled)
        self.ui.isSavedData.toggled.connect(self.savedataToggled)
        self.ui.isSaveFigure.toggled.connect(self.savefigureToggled)

        # open file automatically
        QTimer.singleShot(100, self.selectFile)  # select file at start

    def saveFile(self):
        self.ui.saveprogress.setProperty("value", 0)
        if self.saveFigure:
            if self.calculateg2:
                if not self.f_mapOpen:
                    self.fig = plt.figure(figsize=(10,8))
                    self.mainAxes = self.fig.add_subplot(111)
                    self.f_mapOpen = True
                else:
                    self.mainAxes.clear()  # fresh start for replotting
                self.mainAxes.plot(self.delay*1e9 ,self.coincidences,color='red')
                self.mainAxes.set_xlabel('Delay (ns)')
                self.mainAxes.set_xlim([-30,30])
                self.mainAxes.set_ylabel('Coincidence')
                self.mainAxes.set_ylim(bottom=0)
                self.fig.savefig(self.filename[:-4]+'_%.0fps_g2.png'%(float(self.ui.g2timeres.text())),dpi=100,bbox_inches='tight',facecolor='white')
                plt.close('all')
                self.f_mapOpen = False
            if self.calculatetimetrace:
                if not self.f_mapOpen:
                    self.fig = plt.figure(figsize=(10,8))
                    self.mainAxes = self.fig.add_subplot(111)
                    self.f_mapOpen = True
                else:
                    self.mainAxes.clear()  # fresh start for replotting
                valid = self.countsch0!=0
                self.mainAxes.plot(self.timetracex[valid],self.countsch0[valid],label='ch0')
                self.mainAxes.plot(self.timetracex[valid],self.countsch1[valid],label='ch1')
                self.mainAxes.set_xlabel('Time (s)')
                self.mainAxes.set_ylabel('Counts/%.3fs'%self.readparameters()[0])
                self.mainAxes.legend(loc='best')
                self.mainAxes.set_ylim(bottom=0)
                self.fig.savefig(self.filename[:-4]+'_timetrace.png',dpi=100,bbox_inches='tight',facecolor='white')
                plt.close('all')
                self.f_mapOpen = False

        if self.saveData:
            saved_file = self.filename[:-4]+'_%.0fps_g2__%.2gms_timetrace_coincidences.npz'%(float(self.ui.g2timeres.text()),1e3*float(self.ui.timetraceres.text()))
            
            np.savez(saved_file,timetracex=self.timetracex,countsch0=self.countsch0,countsch1=self.countsch1,g2_delay=self.delay,g2_coin=self.coincidences)
        if self.saveData or self.saveFigure:
            self.ui.saveprogress.setProperty("value", 100)

    def readparameters(self):
        try:
            return [float(self.ui.timetraceres.text()),float(self.ui.g2timeres.text()),float(self.ui.g2timewindow.text()),float(self.ui.zerodelaypos.text()),self.ui.g2algo.currentText(),self.ui.g2postselectionopts.currentText(),float(self.ui.postselect0.text()),float(self.ui.postselect1.text())]
        except Exception as ep:
            errorMessageWindow(self,'ReadParameters Error','%s'%(ep))
    
    def mapClose(self, event):
        self.f_mapOpen = False

    def plotFigure(self):
        plt.ion()
        both_figure = False
        if not self.f_mapOpen:
            if not (self.calculateg2 and self.calculatetimetrace):
                self.fig = plt.figure(figsize=(20,8))
                self.mainAxes = self.fig.add_subplot(121)
            else:
                self.fig = plt.figure(figsize=(20,8))
                self.mainAxes = self.fig.add_subplot(121)
                self.mainAxes1 = self.fig.add_subplot(122)
                both_figure = True
            # self.mainAxes1 = plt.subplot2grid((4, 1), (1, 0), rowspan=1, colspan=1) # self.fig.add_subplot(111)
            # self.mainAxes2 = plt.subplot2grid((4, 1), (2, 0), rowspan=2, colspan=1) # self.fig.add_subplot(111)
            self.f_mapOpen = True
        else:
            self.mainAxes.clear()  # fresh start for replotting
            try:
                self.mainAxes1.clear()
            except:
                pass

        if self.calculateg2:
            self.mainAxes.plot(self.delay*1e9 ,self.coincidences,color='red')
            self.mainAxes.set_xlabel('Delay (ns)')
            self.mainAxes.set_xlim([-30,30])
            self.mainAxes.set_ylabel('Coincidence')
            self.mainAxes.set_ylim(bottom=0)
        if self.calculatetimetrace:
            valid = self.countsch0!=0
            if both_figure:
                out = self.mainAxes1
            else:
                out = self.mainAxes
            out.plot(self.timetracex[valid],self.countsch0[valid],label='ch0')
            out.plot(self.timetracex[valid],self.countsch1[valid],label='ch1')
            out.set_xlabel('Time (s)')
            out.set_ylabel('Counts/%.3fs'%self.readparameters()[0])
            out.set_ylim(bottom=0)
            out.legend(loc='best')
            

        plt.tight_layout()
        self.fig.canvas.mpl_connect('close_event', self.mapClose)
        self.fig.canvas.draw()
        self.fig.canvas.show()
        
    def analysePTU(self,event):
        _parameters = self.readparameters()
        _timetraceres = _parameters[0]
        _g2res =  _parameters[1]*1e-12
        _g2algo = _parameters[4]
        if _g2algo =='Symmetric':
            factor = 2e-9
        else:
            factor = 1e-9
        _g2window = factor * _parameters[2]
        _zero_delay = _parameters[3]*1e-9
        _postselectionopts = _parameters[5]
        _post0 = _parameters[6]
        _post1 = _parameters[7] 
        self.ui.calculatebar.setProperty("value", 0)
        self.ui.saveprogress.setProperty("value", 0)
        try:
            with PTUfile(self.filename) as ptu_file:
                ptu_file.print_header()
                ptu_meas = PTUmeasurement(ptu_file)
                if self.calculatetimetrace:
                    self.timetracex, self.countsch0, self.timetrace_recnum =\
                        ptu_meas.timetrace(resolution=_timetraceres,channel=0, n_threads=4)
                    _, self.countsch1, self.timetrace_recnum =\
                        ptu_meas.timetrace(resolution=_timetraceres, channel=1,n_threads=4)
                if self.calculateg2:
                    threshold = [_post0,_post1]
                    constraint = _postselectionopts
                    if constraint=='None':
                        recnum_post_selec_ranges  = None
                    else:
                        try:
                            post_selec_ranges, recnum_post_selec_ranges = construct_postselect_vector(self.countsch0,self.timetrace_recnum,threshold,constraint)

                            if (len(recnum_post_selec_ranges))==1:
                                recnum_post_selec_ranges = None
                        except Exception as ep:
                            errorMessageWindow(self,"G2 Post-Selection Error","%s\n You need to run time trace along with G2 to do post-selection."%ep)
                            recnum_post_selec_ranges = None
                    self.delay, self.coincidences = ptu_meas.calculate_g2(_g2window,_g2res,
                                                         post_selec_ranges=recnum_post_selec_ranges ,
                                                         n_threads=4,
                                                         mode=_g2algo)
                    self.delay -= _zero_delay 
            self.ui.calculatebar.setProperty("value", 100)
        except Exception as ep:
            errorMessageWindow(self,'Analyse PTU Error','%s\n Check if the file format is correct (.ptu, T2 mode).'%ep)
        
        if self.saveFigure or self.saveData:
            self.saveFile()
            plt.close('all')
        else:
            self.plotFigure()

    def timetraceToggled(self, event):
        self.calculatetimetrace = self.ui.istimetrace.isChecked()
        
    def g2Toggled(self, event):
        self.calculateg2 = self.ui.isG2.isChecked()

    def savedataToggled(self, event):
        self.saveData =  self.ui.isSavedData.isChecked()
        

    def savefigureToggled(self, event):
        self.saveFigure = self.ui.isSaveFigure.isChecked()

    def selectFile(self, event=0):
        """ 
        Generates the selection dialog window until cancelled or valid selected file.
        This function calls self.openFile() to open the file and validate it.
        """
        while True:
            
            self.filename = QFileDialog.getOpenFileName(self, 'Open T2 PTU file')[0]
            if self.filename == '':  # user cancelled
                break
            else:
                try: 
                    self.ui.filepath.setText(self.filename)
                except Exception as ep:
                    errorMessageWindow(self, "Error Message", 
                                            "%s"%ep)
                finally:
                    break
def errorMessageWindow(parentWindow, winTitle, winText):
    """
    Displays a QT error message box, with title, text and OK button
    
    Args:
        parentWindow (QWidget): Parent widget used to display the error window
        winTitle (str): Text displayed as error window title
        winText (str): Text displayed as error message.
    """
    msg = QtWidgets.QMessageBox(parentWindow)
    msg.setIcon(QtWidgets.QMessageBox.Critical)
    msg.setWindowTitle(winTitle)
    msg.setText(winText)
    msg.exec_()

def main():
    app = QApplication(sys.argv)
    window = AnalysePTU_window()
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()