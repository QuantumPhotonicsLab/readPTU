import numpy as np
import os 
from numpy import *
from scipy import optimize
import h5py
# import pylab

# taken from the scipy fitting cookbook:
class Parameter:
    def __init__(self, value, name=''):
        self.value = value
        self.name = name

    def set(self, value):
        self.value = value

    def __call__(self):
        return self.value

###############################################################################
# fitting a 1d function
################################################################################

# TODO 
# - fitmethods should be classes
# - fit should actually also be a class, and we want a simple function as
# wrapper for interactive work; then we still need sth better for fixing,
# though; good maybe: generally identify parameters by names

# added capability to fit data using weights (error bars in datapoints): 
# errors can be passed by a list err_y - Cristian 24/11/2014
# THT: I had to remove this because it crashes many measurements, including the optimizOr




def fit1d(x, y, fitmethod, *arg, **kw):
    """
    example: from analysis.lib.fitting import fit,common
             x=np.array([0,1,2,3,4])
             y=np.array([2,12,22,32,42])
             fit_result=fit.fit1d(x,y,common.fit_line,2,8,ret=True,
                    fixed=[0],do_print=True)
    Returns a dictionary with the results if ret=True (regardles of the sit success), None otherwise.
             
    """
    # process known kws
    do_print = kw.pop('do_print', False)
    ret = kw.pop('ret', True)
    fixed = kw.pop('fixed', [])
    VERBOSE= kw.pop('VERBOSE',False)


    # err_y = kw.pop ('err_y', None)
	#if False :
	#    if (len(err_y) != len(y)):
	#    	print 'Data and error arrays have non-matching lengths!'
	#    	err_y = None

    # if (len(err_y) != len(y)):
    # 	print 'Data and error arrays have non-matching lengths!'
    # 	err_y = None

    # use the standardized fitmethod: any arg is treated as initial guess
    if fitmethod != None:
        p0, fitfunc, fitfunc_str = fitmethod(*arg)
    else:
        p0 = kw.pop('p0')
        fitfunc = kw.pop('fitfunc')
        fitfunc_str = kw.pop('fitfunc_str', '')        
    
    # general ability to fix parameters
    fixedp = []
    for i,p in enumerate(p0):
        if i in fixed:
            fixedp.append(p)
    for p in fixedp:
        p0.remove(p)
   
    # convenient fitting method with parameters; see scipy cookbook for details
    def f(params):
        i = 0
        for p in p0:
            p.set(params[i])
            i += 1

        # if (err_y != None):
        # 	return ((y-fitfunc(x))/(err_y))
        # else:

        return y - fitfunc(x)


    if x is None: x = arange(y.shape[0])
    p = [param() for param in p0]
    
    # do the fit and process
    p1, cov, info, mesg, success = optimize.leastsq(f, p, full_output=True, maxfev=len(x)*100)
    if not success or cov == None: # FIXME: find a better solution!!!
        success = False
        if VERBOSE:
            print 'ERROR: Fit did not converge !'
            print 'reason: ',mesg
        # return success    #commented out by THT and MA because it bvreaks all old automatic fitting code. 160802
    
    result = result_dict(p1, cov, info, mesg, success, x, y, p0, 
            fitfunc, fitfunc_str)
    # package the result neatly
    if do_print and success:
        print_fit_result(result)
    if ret: #and success:
        return result


def fit2d((meshx,meshy),z,fitmethod, *arg, **kw):
    """
    This function fits 2d data. It does so by raveling the data down to a 1 d array.
    """
    # process known kws
    do_print = kw.pop('do_print', False)
    ret = kw.pop('ret', False)
    fixed = kw.pop('fixed', [])
    VERBOSE= kw.pop('VERBOSE',False)

    # use the standardized fitmethod: any arg is treated as initial guess
    if fitmethod != None:
        p0, fitfunc, fitfunc_str = fitmethod(*arg)
    else:
        p0 = kw.pop('p0')
        fitfunc = kw.pop('fitfunc')
        fitfunc_str = kw.pop('fitfunc_str', '')        
 
 
    # general ability to fix parameters
    fixedp = []
    for i,p in enumerate(p0):
        if i in fixed:
            fixedp.append(p)
    for p in fixedp:
        p0.remove(p)

    if ((meshx is None) or (meshy is None)): meshx,meshy = mgrid[0:z.shape[0], 0:z.shape[1]]
   
    # convenient fitting method with parameters; see scipy cookbook for details
    def f(params):
        i = 0
        for p in p0:
            p.set(params[i])
            i += 1

        return ravel(z) - ravel(fitfunc(ravel(meshx),ravel(meshy)))



    p = [param() for param in p0]

    # do the fit and process
    p1, cov, info, mesg, success = optimize.leastsq(f, p, full_output=True, maxfev=len(ravel(meshx))*20)
    if not success or cov == None: # FIXME: find a better solution!!!
        if VERBOSE:
            print 'ERROR: Fit did not converge !'
            print 'reason: ',mesg
        return success
    
    result = result_dict(p1, cov, info, mesg, success, meshx, z, p0, 
            fitfunc, fitfunc_str)    #if this 2d function becomes more commonly used, the dictionary should also have meshy as a key.

    # package the result neatly
    if do_print and success:
        print_fit_result(result)
    if ret and success:
        return result


###############################################################################
# tools, for formatting, printing, etc.
###############################################################################

# put all the fit results into a dictionary, calculate some more practical 
# numbers
def result_dict(p1, cov, info, mesg, success, x, y, p0, fitfunc, fitfunc_str):
    chisq = 1
    dof = 1
    error_dict = {}
    error_list = []
    params_dict = {}
    
    # print cov, success, mesg, info
    if success:
        chisq = sum(info['fvec']*info['fvec'])
        dof = len(y)-len(p0)
        for i,pmin in enumerate(p1):
            error_dict[p0[i].name] = sqrt(cov[i,i])*sqrt(chisq/dof)
            #print chisq
            #print dof   
            error_list.append(sqrt(cov[i,i])*sqrt(chisq/dof))
            params_dict[p0[i].name] = pmin

    result = {
        'success' : success,
        'params' : p1,
        'params_dict' : params_dict,
        'chisq': chisq,
        'dof': dof,
        'residuals_rms': sqrt(chisq/dof),
        'reduced_chisq': chisq/dof,
        'error' : error_list,
        'error_dict' : error_dict, 
        'cov' : cov,
        'p0' : p0,
        'fitfunc' : fitfunc,
        'fitfunc_str' : fitfunc_str,
        'x' : x,
        'y' : y,
        }
    
    return result

def fitparam_table(result):
    """
    returns the list of parameters and a table with columns 
    fitted value and the uncertainty (1 row for each parameter).
    """
    params = np.empty((0,2))
    param_names = []

    for n in result['params_dict']:
        params = np.vstack((params, 
            np.array([result['params_dict'][n], result['error_dict'][n]])))
        param_names.append(n)

    return param_names, params


# convenient for pylab usage (or other interactive)
def do_fit_func(fitfunc, p0, y, x):
    result = _fit_return(fit(fitfunc, p0, y, x), y, p0, fitfunc(x))
    print_fit_result(result)
    return result

# make a string that contains the fit params in a neat format
def str_fit_params(result):
    
    # uncertainties are calculated as per gnuplot, "fixing" the result
    # for non unit values of the reduced chisq.
    # values at min match gnuplot
    
    str = "fitted parameters at minimum, with 68% C.I.:\n"
    for i,pmin in enumerate(result['params']):
        str += "%2i %-10s %12f +/- %10f\n" % \
            (i, result['p0'][i].name, pmin, result['error'][i])
    return str

def str_correlation_matrix(result):
    str = "correlation matrix:\n"
    str += "               "
    for i in range(len(result['p0'])): 
        str+= "%-10s" % (result['p0'][i].name,)
    str += "\n"
    
    for i in range(len(result['params'])):
        str += "%10s" % result['p0'][i].name
        for j in range(i+1):
            str+= "%10f" % \
                (result['cov'][i,j] / \
                     sqrt(result['cov'][i,i] * result['cov'][j,j]),)
        str+='\n'

    return str
    
def print_fit_result(result):
    if result == False:
       print "Could not fit data"
       return
    
    print "Converged with chi squared ", result['chisq']
    print "degrees of freedom, dof ", result['dof']
    print "RMS of residuals (i.e. sqrt(chisq/dof)) ", \
        sqrt(result['chisq']/result['dof'])
    print "Reduced chisq (i.e. variance of residuals) ", \
        result['chisq']/result['dof']
    print

    print str_fit_params(result)
    print str_correlation_matrix(result) 

def write_to_file(fitresult,folder, filename='fit_results.txt', fitname = 'Name not specified'):
    print 'Writting to File!'
    text_file = open( os.path.join(folder, filename), 'w')  
    print 'path joined'
    if fitresult == False:
        text_file.write("Could not fit data") 
    else: 
        text_file.write('''
Fit results of: %s 

        Converged with chi squared: %s
        Degrees of freedom, dof %s 
        RMS of residuals (i.e. sqrt(chisq/dof)) %s 
        Reduced chisq (i.e. variance of residuals) %s

                ''' %(fitname,fitresult['chisq'],fitresult['dof'],sqrt(fitresult['chisq']/fitresult['dof']),fitresult['chisq']/fitresult['dof']) )
    
        text_file.write(str_fit_params(fitresult))
        text_file.write(str_correlation_matrix(fitresult) )
    print 'Writing fit results to file (%s) at (%s) succesfull' %(filename,folder) 


    text_file.close() 
    
def write_to_hdf(fitresult,fp):
    f=h5py.File(fp, 'a')
    g = f.create_group('fit_result')
    for attr_key in ['success','dof', 'chisq','fitfunc_str','residuals_rms','reduced_chisq']:
        g.attrs[attr_key] = fitresult[attr_key]
    g['x'] = fitresult['x']
    g['y'] = fitresult['y']
    g_p = g.create_group('params_dict')
    for k in fitresult['params_dict']:
        g_p[k] = fitresult['params_dict'][k]
    g_e = g.create_group('error_dict')
    for k in fitresult['error_dict']:
        g_e[k] = fitresult['params_dict'][k]
    f.close()
