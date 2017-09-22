import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib.tools import plot
from analysis.lib.m2.ssro import sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr; reload(fit)


def analyze_dark_esr(guess_ctr, guess_splitN,
guess_offset = 1,
guess_width = 0.2e-3,
guess_amplitude = 0.3,
min_dip_depth = 0.9, 
timestamp = None,
add_folder = None,
ret='f0',
ssro_calib_folder='',
do_save=True,
sweep_direction = 'right',
**kw):

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('DarkESR')

    if add_folder !=None:
        folder = add_folder
    print 'analysis script!'
    print folder
    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    a.get_electron_ROC()

    x = a.sweep_pts # convert to MHz
    y = a.p0.reshape(-1)[:]

    # Find the esr resonance
    j=0
    print 'j = '+str(j)
    print folder
    print y[21]
    k = len(y)
    # min_dip_depth = np.amin(y)
    # min_index = np.argmin(y)
    print 'min_dip_depth = ' + str(min_dip_depth)
    
    ### Option to make the dip search sweep towards left or right, usefull in case of N polarization
    if sweep_direction == 'left':
        y1 = y[::-1]
        x1 = x[::-1]
        guess_splitN = -1*guess_splitN
    elif sweep_direction == 'right':
        y1 = y
        x1 = x

    while y1[j]>min_dip_depth and j < len(y)-2:  #y[j]>0.93*y[j+1]: # such that we account for noise
        k = j
        j += 1
    #j = len(y)-2
    if k > len(y)-3:
        print 'Could not find dip'
        return
    else:
        guess_ctr = x1[k]+ guess_splitN #convert to GHz and go to middle dip

        print 'guess_ctr= '+str(guess_ctr)
        print 'k'+str(k)
    ## I added this to be more robust for SSRO calibration.Please monitor if this is better - Machiel may-2014

    fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
            guess_amplitude, guess_width, guess_ctr,
            (3, guess_splitN),
            do_print=True, ret=True, fixed=[])
    print 'fit finished'
    if do_save:
        print 'saving data and fit'
        if fit_result:
           
            fig, ax = plt.subplots(1,1)
            plot.plot_fit1d(fit_result, np.linspace(min(x), max(x), 1000), ax=ax, plot_data=True, **kw)
            plt.savefig(os.path.join(folder, 'darkesr_analysis.png'),
            format='png')
            #plt.show()
            plt.close(fig)
    if ret == 'f0':
        f0 = fit_result['params_dict']['x0']
        u_f0 = fit_result['error_dict']['x0']
        return f0, u_f0

def analyze_dark_esr_2(guess_ctr, guess_splitN,
guess_offset = 1,
guess_width = 0.2e-3,
guess_amplitude = 0.3,
min_dip_depth = 0.86, 
timestamp = None,
add_folder = None,
ret='f0',
ssro_calib_folder='',
do_save=True,
sweep_direction = 'right',
**kw):

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('DarkESR')

    if add_folder !=None:
        folder = add_folder
    print 'analysis script!'
    print folder
    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    a.get_electron_ROC()

    x = a.sweep_pts # convert to MHz
    y = a.p0.reshape(-1)[:]

    # Find the esr resonance
    j=0
    print 'j = '+str(j)
    print folder
    print y[21]
    k = len(y)
    print 'min_dip_depth = ' + str(min_dip_depth)
    
    ### Option to make the dip search sweep towards left or right, usefull in case of N polarization
    if sweep_direction == 'left':
        y1 = y[::-1]
        x1 = x[::-1]
        guess_splitN = -1*guess_splitN
    elif sweep_direction == 'right':
        y1 = y
        x1 = x

    while y1[j]>min_dip_depth and j < len(y)-2:  #y[j]>0.93*y[j+1]: # such that we account for noise
        k = j
        j += 1
    #j = len(y)-2
    if k > len(y)-3:
        print 'Could not find dip'
        return
    else:
        guess_ctr = x1[k]+ guess_splitN #convert to GHz and go to middle dip

        print 'guess_ctr= '+str(guess_ctr)
        print 'k'+str(k)
    ## I added this to be more robust for SSRO calibration.Please monitor if this is better - Machiel may-2014

    fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
            guess_amplitude, guess_width, guess_ctr,
            (3, guess_splitN),
            do_print=True, ret=True, fixed=[])
    print 'fit finished'
    if do_save:
        print 'saving data and fit'
        if fit_result:
         
            fig, ax = plt.subplots(1,1)
            print min(x)
            print max(x)
            plot.plot_fit1d(fit_result, np.linspace(min(x), max(x), 1000), ax=ax, plot_data=True, print_info=False, figsize=(6,4.7), **kw)
            ax.set_xlim([1.742,1.751])
            ax.set_ylim([0.7,1.05])
            plt.xticks([1.742,1.746,1.750],['1.742','1.746','1.750'])
            plt.yticks([0.7,0.8,0.9,1])
            ax.tick_params(axis='x', which='major', labelsize=20)
            ax.tick_params(axis='y', which='major', labelsize=20)
            ax.set_xlabel('MW frq (GHz)',fontsize = 20)
            ax.set_ylabel(r'Fidelity wrt. $|0\rangle$',fontsize = 20)
            plt.tight_layout()
            plt.savefig(os.path.join(folder, 'darkesr_analysis.pdf'),
            format='pdf')
            #plt.show()
            plt.close(fig)
    if ret == 'f0':
        f0 = fit_result['params_dict']['x0']
        u_f0 = fit_result['error_dict']['x0']
        return f0, u_f0


def analyze_dark_esr_single( 
guess_offset = 1,
guess_width = 0.2e-3,
guess_amplitude = 0.3,
timestamp = None,
add_folder = None,
ret='f0',
ssro_calib_folder='',
**kw):

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('DarkESR')


    if add_folder !=None:
        folder = add_folder

    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    a.get_electron_ROC()

    x = a.sweep_pts # convert to MHz
    y = a.p0.reshape(-1)[:]

    guess_ctr = x[y.argmin()]
    print 'guess_ctr = '+str(guess_ctr)
    ## I added this to be more robust for SSRO calibration.Please monitor if this is better - Machiel may-2014
    guess_offset=np.average(y)
    dip_threshold=guess_offset-1.5*np.std(y)
    print guess_offset
    print dip_threshold

    # if min(y) > dip_threshold:
    #     print 'Could not find dip'
    #     return

    fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
            guess_amplitude, guess_width, guess_ctr,
            do_print=True, ret=True, fixed=[])

    if ret == 'f0':
        f0 = fit_result['params_dict']['x0']
        u_f0 = fit_result['error_dict']['x0']
        return f0, u_f0


def analyze_dark_esr_double( 
guess_offset = 1,
guess_A_min = 0.3,
guess_A_plus = 0.3,
guess_width = 0.2e-3,
guess_amplitude = 0.3,
guess_Csplit = 0.100e-3,
timestamp = None,
add_folder = None,
ret='f0',
ssro_calib_folder='',
do_ROC = True,
ret_folder = False,
do_plot = True,
**kw):

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('DarkESR')

    print folder

    if add_folder !=None:
        folder = add_folder

    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    if do_ROC == True:
        a.get_electron_ROC()
        y = a.p0.reshape(-1)[:]
    else:
        y = a.get_readout_results('ssro')

    x = a.sweep_pts # convert to MHz
    

    ### create a guess for the center frequency
    guess_ctr1 = x[y.argmin()]+guess_Csplit*1.2
    guess_ctr2 = x[y.argmin()]-guess_Csplit*1.2 #

    # guess_ctr1 =x[len(x)/2.]
    print 'guess_ctr1 = '+str(guess_ctr1)
    print 'guess_ctr2 = '+str(guess_ctr2)

    
    ### First fit attempt: guess_ctr1

    ### fitfunction
    A_min = fit.Parameter(guess_A_min, 'A_min')
    A_plus = fit.Parameter(guess_A_plus, 'A_plus')
    o = fit.Parameter(guess_offset, 'o')
    ctr = fit.Parameter(guess_ctr1, 'ctr')
    width = fit.Parameter(guess_width, 'width')
    Csplit = fit.Parameter(guess_Csplit, 'Csplit')

    def fitfunc(x):
        return o() - A_min()*np.exp(-((x-(ctr()-Csplit()))/width())**2) \
                - A_plus()*np.exp(-((x-(ctr()+Csplit()))/width())**2) \


    print 'running fit1'
    fit_result = fit.fit1d(x, y, None, p0 = [A_min, A_plus, o, ctr, width, Csplit],
            fitfunc = fitfunc, do_print=False, ret=True, fixed=[])
    do_another_fit = False

    if fit_result['success'] == False:
            do_another_fit = True
    elif fit_result['error_dict']['ctr'] > 8e-6:
            do_another_fit = True

    if do_another_fit == True:
            print 'first fit failed'
            ### Second fit attempt: guess_ctr2
            A_min = fit.Parameter(guess_A_min, 'A_min')
            A_plus = fit.Parameter(guess_A_plus, 'A_plus')
            o = fit.Parameter(guess_offset, 'o')
            ctr = fit.Parameter(guess_ctr2, 'ctr')
            width = fit.Parameter(guess_width, 'width')
            Csplit = fit.Parameter(guess_Csplit, 'Csplit')

            def fitfunc(x):
                return o() - A_min()*np.exp(-((x-(ctr()-Csplit()))/width())**2) \
                        - A_plus()*np.exp(-((x-(ctr()+Csplit()))/width())**2) \

            print 'running fit2'
            fit_result = fit.fit1d(x, y, None, p0 = [A_min, A_plus, o, ctr, width, Csplit],
                    fitfunc = fitfunc, do_print=False, ret=True, fixed=[])
    
    if do_plot == True:        
        fig, ax = plt.subplots(1,1)
        plot.plot_fit1d(fit_result, np.linspace(min(x), max(x), 1000), ax=ax, plot_data=True, **kw)

        plt.savefig(os.path.join(folder, 'darkesr_analysis.png'), format='png')
        #plt.show()
        plt.close(fig)

    if ret == 'f0' and ret_folder == False:
        if fit_result['success'] == False:
            print 'Fit failed, returned dummy values for frequency'
            f0 = 0
            u_f0 = 1e-3         #large arbitrary uncertainty makes sure the data is filtered out later
            return f0, u_f0
        else: 
            f0 = fit_result['params_dict']['ctr']
            u_f0 = fit_result['error_dict']['ctr']
            return f0, u_f0
      

    elif ret == 'f0' and ret_folder == True:
        return f0, u_f0, folder
    












