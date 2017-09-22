import numpy
import scipy
import scipy.linalg.basic as slb
from  scipy.optimize import leastsq

def AOMfunction(x,a,xc,k):
    y = a * numpy.exp( -numpy.exp( -k*(x-xc)))
    return y

def SaturationFunction(x,max,p_sat):
    y = max * x / (p_sat + x)
    return y

def lorentzian(x,x0,A,w,h):
    '''
    x0 is center of peak
    A is amplitude of peak
    w is Full Width Half Maximum
    h is baselevel
    '''
    y = h + A / (1.0 + ((x - x0)/(w/2.0))**2)
    return y

def lorentzian2(x,x0,A,w,h):
    '''
    x0 is center of peak
    A is amplitude of peak
    w is Full Width Half Maximum
    h is baselevel
    '''
    y = h + A / (1.0 + ((x - x0)/(w/2.0))**2)
    return y

def double_lorentzian(x,xc1, xc2,A1, A2, w1, w2 ,h):
    '''
    x0 is center of peak
    A is amplitude of peak
    w is Full Width Half Maximum
    h is height of baselevel
    '''
    y =h + A1 / (1.0 + ((x - xc1)/(w1/2.0))**2) + A2 / (1.0 + ((x - xc2)/(w2/2.0))**2)

    return y

def gaussian(x,x0,A,sd):
    '''
    x0 is center of peak
    A is amplitude of peak
    sd is standard dev
    '''
    y = A * numpy.exp(-((x-x0)/sd)**2/2.0)
    return y

def gaussian2(x,x0,A,h,sd):
    '''
    x0 is center of peak
    A is amplitude of peak
    h is baselevel
    sd is standard dev
    '''
    y = h + A * numpy.exp(-((x-x0)/sd)**2/2.0)
    return y

def damped_sin(x,y0,A,f,d):
    '''
    y0 is final level
    A is initial amplitude
    f is frequency
    d is damping
    '''
    y = y0 + numpy.exp(-x/d) * A * numpy.sin(2*numpy.pi*f*x)
    return y

def sin(x,A,f,phi,y0):
    '''
    A is initial amplitude
    f is frequency
    phi is phase
    y0=offset
    '''
    y = y0 + A * numpy.sin(2*numpy.pi*f*x+phi)
    return y

def cos(x,A,f,phi=0,y0=0):
    '''
    A is initial amplitude
    f is frequency
    phi is phase
    y0=offset
    '''
    y = y0 + A * numpy.cos(2*numpy.pi*f*x+phi)
    return y

def cos_fixed_phase(x,omega,y0):
    y = y0 - 0.5 * numpy.cos(omega*x)
    return y

def cos_fixed_phase_offset(x,omega):
    y = - 0.5 * numpy.cos(omega*x)
    return y

def damped_cos_fixed(x,omega,d,phi):
    y = numpy.exp(-d*x) * (-0.5) * (numpy.cos((omega + phi)*x))
    return y

def damped_cos_fixed_0(x,omega,d):
    y = numpy.exp(-d*x) * (-0.5) * numpy.cos(omega*x)
    return y


def abs_linear(x,x0, mr,ml):
    y = mr*(x-x0)*(numpy.sign(x-x0)+1)/2+ml*(x0-x)*(numpy.sign(x0-x)+1)/2
    return y

def abs_linear_0(x,x0, m):
    y = m*abs(x-x0)
    return y

def linear(x,x0,m):
    y = m*(x-x0)
    return y


def damped_cos(x,y0,A,f,d):
    '''
    y0 is final level
    A is initial amplitude
    f is frequency
    d is damping
    '''
    y = y0 + numpy.exp(-x/d) * A * numpy.cos(2*numpy.pi*f*x)
    return y

def exp_decay(x, x0, y0, A, T):
    '''
    y0 is saturation value
    T is decay rate
    '''
    y = y0 + A*numpy.exp(-(x-x0) / T)
    return y

def integer_phase(n , phase, A, y0):
    '''
    y0 is value at 0
    phase is the constant phase per increase n
    n is x-coord
    A is amplitude
    '''
    y = y0 - A*sin(n*phase)
    return y

def _residuals_AOMfunction(p,x,y):
    err = y - AOMfunction(x,p[0],p[1],p[2])
    return err

def fit_AOMfunction(xdata,ydata,p0):
    plsq = leastsq(_residuals_AOMfunction, p0, args=(xdata,ydata))
    return plsq

def _residuals_SaturationFunction(p,x,y):
    err = y - SaturationFunction(x,p[0],p[1])
    return err

def fit_SaturationFunction(xdata,ydata,p0):
    plsq = leastsq(_residuals_SaturationFunction, p0, args=(xdata,ydata))
    return plsq

def _residuals_gaussian(p, x, y):
    err = y - gaussian(x,p[0],p[1],p[2])
    return err

def fit_gaussian(xdata, ydata, p0):
    '''
    p0 is a vector of starting estimates:
        
        x0
        A
        sd
    '''
    plsq = leastsq(_residuals_gaussian, p0, args=(xdata,ydata))
    return plsq

def _residuals_gaussian2(p, x, y):
    err = y - gaussian2(x,p[0],p[1],p[2],p[3])
    return err

def fit_gaussian2(xdata, ydata, p0):
    '''
    p0 is a vector of starting estimates:
        
        x0
        A
        h
        sd
        
    '''
    plsq = leastsq(_residuals_gaussian2, p0, args=(xdata,ydata))
    return plsq


def fit_plot_gaussian2(xdata, ydata, p0):
    qt.plots.clear()
    plot(xdata,ydata, name='fit_gaussian')
    fits = fit_gaussian2(xdata, ydata, p0)
    fits = fits[0]
    fitdata = gaussian2(xdata, fits[0], fits[1], fits[2], fits[3])
    plot(xdata, fitdata, name='fit_gaussian')
    return fits

def _residuals_double_gaussian(p, x, y):
    err = y - (gaussian(x,p[0],p[1],p[2]) + gaussian(x,p[3],p[4],p[5]))
    return err

def fit_double_gaussian(xdata, ydata, p0):
    '''
    p0 is a vector of starting estimates:
        
        x0_1
        A_1
        sd_1
        x0_2
        A_2
        sd_2
    '''
    plsq = leastsq(_residuals_double_gaussian, p0, args=(xdata,ydata), full_output=1, col_deriv=1)
    return plsq


def _residuals_lorentzian(p, x, y):
    err = y - lorentzian(x,p[0],p[1],p[2], p[3])
    return err

def fit_lorentzian(xdata, ydata, p0):
    '''
    p0 is a vector of the starting estimates:

        x0
        A
        w
        h
    '''
    plsq = leastsq(_residuals_lorentzian, p0, args=(xdata,ydata))
    return plsq[0]

def _residuals_lorentzian2(p, x, y):
    err = y - lorentzian2(x,p[0],p[1],p[2],p[3])
    return err

def fit_lorentzian2(xdata, ydata, p0):
    '''
    p0 is a vector of the starting estimates:

        x0 - peak position
        A - peak amplitude
        w - peak width
        b - baselevel

    '''
    plsq,cov,info,mesg,success = leastsq(_residuals_lorentzian2, p0, args=(xdata,ydata), full_output = 1)
    dof = len(xdata) - len(plsq)
    error = len(plsq)*[0]
    chisq=sum(info["fvec"]*info["fvec"])
    for i in arange(plsq):
        error[i] = sqrt(cov[i,i])*sqrt(chisq/dof)
    return plsq, error

def _residuals_double_lorentzian(p, x, y):
    err = y - double_lorentzian(x,p[0],p[1],p[2],p[3],p[4],p[5], p[6])
    return err


def fit_double_lorentzian(xdata, ydata, p0):
    '''
    p0 is a vector of the starting estimates:

        xc1 - peak position
        xc2 - peak position
        A1 - peak amplitude
        A2 - peak amplitude
        w1 - peak width
        w2 - peak width
        b - baselevel

    '''
    plsq,cov,info,mesg,success = leastsq(_residuals_double_lorentzian, p0, args=(xdata,ydata), full_output = 1)
    print plsq
    print cov
    dof = len(xdata) - len(plsq)
    error = len(plsq)*[0]
    print 'error %s'%error
    chisq=sum(info["fvec"]*info["fvec"])
    #for i in arange(len(plsq)):
        #error[i] = sqrt(cov[i,i])*sqrt(chisq/dof)
    return plsq, error


def fit_plot_lorentzian2(xdata, ydata, p0):
    qt.plots.clear()
    plot(xdata,ydata, name='fit_lorentzian')
    fits = fit_lorentzian2(xdata, ydata, p0)
    fitdata = lorentzian2(xdata, fits[0], fits[1], fits[2], fits[3])
    plot(xdata, fitdata, name='fit_lorentzian')
    return fits

def _residuals_exponential(xdata, ydata, p):
    err = y - exp_decay(xdata, p[0],p[1],p[2],p[3])
    return err

def fit_exp_decay(xdata, ydata, p0):
    '''
    p0 is a vector of the starting estimates:

        x0 - starting position
        y0 - level at saturation
         A - amplitude
         T - decay rate
        

    '''
    plsq = leastsq(_residuals_exponential, p0, args=(xdata,ydata,x0))
    return plsq[0]

def _residuals_damped_sin(p, x, y):
    err = y - damped_sin(x,p[0],p[1],p[2],p[3])
    return err

def fit_damped_sin(xdata, ydata, p0):
    '''
    p0 is vector of starting estimates:
        
        y0
        A
        f
        d
    '''
    plsq = leastsq(_residuals_damped_sin, p0, args=(xdata,ydata))
    return plsq

def _residuals_damped_cos(p, x, y):
    err = y - damped_cos(x,p[0],p[1],p[2],p[3])
    return err

def fit_damped_cos(xdata, ydata, p0):
    '''
    p0 is vector of starting estimates:
        
        y0
        A
        f
        d
    '''
    plsq = leastsq(_residuals_damped_cos, p0, args=(xdata,ydata))
    return plsq

def _residuals_sin(p, x, y):
    err = y - sin(x,p[0],p[1],p[2],p[3])
    return err

def _residuals_cos(p, x, y):
    err = y - cos(x,p[0],p[1],p[2],p[3])
    return err

def _residuals_cos_fixed_phase(p, x, y):
    err = y - cos_fixed_phase(x,p[0],p[1])
    return err

def _residuals_cos_fixed_phase_offset(p, x, y):
    err = y - cos_fixed_phase_offset(x,p[0])
    return err

def _residuals_damped_cos_fixed(p, phi, x, y):
    err = y - damped_cos_fixed(x,p[0],p[1],phi)
    return err

def _residuals_damped_cos_fixed_0(p, x, y):
    err = y - damped_cos_fixed_0(x,p[0],p[1])
    #if abs(p[0])<0.01:
        #err=err*2

    return err


def _residuals_damped_cos_fixed_tau(p, tau, x, y):
     err = y - damped_cos_fixed(x,p,tau)
     return err



def _residuals_abs_linear(p, x, y):
    err = y - abs_linear(x,p[0],p[1],p[2])
    return err

def _residuals_abs_linear_0(p, x, y):
    err = y - abs_linear_0(x,p[0],p[1])
    return err

def _residuals_linear(p, x, y):
    err = y - linear(x,p[0],p[1])
    return err


def fit_cos_fixed_phase(xdata,ydata,p0):
    plsq,cov,info,mesg,success = leastsq(_residuals_cos_fixed_phase, p0, args=(xdata,ydata), full_output = 1)
    dof = len(xdata) - len(plsq)
    errors = estimate_errors(plsq, info, cov, dof)
    return plsq, errors


def fit_cos_fixed_phase_offset(xdata,ydata,p0):
    plsq,cov,info,mesg,success = leastsq(_residuals_cos_fixed_phase_offset, p0, args=(xdata,ydata), full_output = 1)
    plsq=[plsq]
    dof = len(xdata) - 1
    errors = estimate_errors(plsq, info, cov, dof)
    return plsq, errors


def fit_abs_linear(xdata,ydata,p0):
    plsq,cov,info,mesg,success = leastsq(_residuals_abs_linear, p0, args=(xdata,ydata), full_output = 1)
    dof = len(xdata) - len(plsq)
    if cov.any() :
        errors = estimate_errors(plsq, info, cov, dof)
    else:
        errors = [100]
    return plsq, errors

def fit_abs_linear_0(xdata,ydata,p0):
    plsq,cov,info,mesg,success = leastsq(_residuals_abs_linear_0, p0, args=(xdata,ydata), full_output = 1)
    dof = len(xdata) - len(plsq)
    if cov is None:
        errors = [100,100]
    else:
        errors = estimate_errors(plsq, info, cov, dof)
    return plsq, errors

def fit_linear(xdata,ydata,p0):
    plsq,cov,info,mesg,success = leastsq(_residuals_linear, p0, args=(xdata,ydata), full_output = 1)
    dof = len(xdata) - len(plsq)
    if cov is None:
        errors = [100,100]
    else:
        errors = estimate_errors(plsq, info, cov, dof)
    return plsq, errors


def fit_damped_cos_fixed(xdata,ydata,p0,phi):
    plsq,cov,info,mesg,success = leastsq(_residuals_damped_cos_fixed, p0, args=(phi,xdata,ydata), full_output = 1)
    dof = len(xdata) - 1
    if cov is None:
        errors = [100,100]
    else:
        errors = estimate_errors(plsq, info, cov, dof)
    return abs(plsq), errors

def fit_damped_cos_fixed_0(xdata,ydata,p0):
    plsq,cov,info,mesg,success = leastsq(_residuals_damped_cos_fixed_0, p0, args=(xdata,ydata), full_output = 1)
    dof = len(xdata) - 1
    if cov is None:
        errors = [100,100]
    else:
        errors = estimate_errors(plsq, info, cov, dof)
    return plsq, errors


def fit_damped_cos_fixed_tau(xdata,ydata,p0,tau):
    plsq,cov,info,mesg,success = leastsq(_residuals_damped_cos_fixed_tau, p0, args=(tau,xdata,ydata), full_output = 1)
    plsq = [plsq]
    dof = len(xdata) - 1
    if cov:
        errors = estimate_errors(plsq, info, cov, dof)
    else: 
        errors = [100]
    return plsq, errors


def fit_sin(xdata, ydata, p0):
    '''
    p0 is vector of starting estimates:
        
        A
        f
        phi
        y0
    '''
    plsq = leastsq(_residuals_sin, p0, args=(xdata,ydata))
    return plsq

def fit_cos(xdata, ydata, p0):
    '''
    p0 is vector of starting estimates:
        
        A
        f
        phi
        y0
    '''
    print p0
    plsq= leastsq(_residuals_cos, p0, args=(xdata,ydata))

    return plsq

def fit_cos_we(xdata, ydata, p0):
    '''
    p0 is vector of starting estimates:
        
        A
        f
        phi
        y0
    '''
    plsq,cov,info,mesg,success = leastsq(_residuals_cos, p0, args=(xdata,ydata), full_output = 1)
    dof = len(xdata) - len(plsq)
    errors = estimate_errors(plsq, info, cov, dof)
    return plsq, errors




def _residuals_sin2(p, x, y):
    err = y - sin(x,p[0],1,p[1],p[2])
    return err

def fit_integer_phase(ndata, ydata, p0):
    '''
    p = [phase, A, y0]
    '''
    plsq,cov,info,mesg,success = leastsq(_residuals_integer_pahse, p0, args=(ndata,ydata), full_output = 1)

def _residuals_integer_phase(n, y, p):
    err = y - integer_phase(n,p[0],p[1],p[2])
    return err


def fit_sin2(xdata, ydata, p0):
    '''
    p0 is vector of starting estimates:
        
        A
        phi
        y0
    '''
    plsq = leastsq(_residuals_sin2, p0, args=(xdata,ydata))
    return plsq

def fit_plot_sin(xdata, ydata, p0,clear=True):
    plot(xdata,ydata, name='fit_sin',clear=clear)
    fits = fit_sin(xdata, ydata, p0)[0]
    fitdata = sin(xdata, fits[0], fits[1], fits[2], fits[3])
    plot(xdata, fitdata, name='fit_sin')
    return fits

def estimate_errors(plsq, fit_info, cov, dof):
    '''
    plsq = fitparams
    fit_info = full_output of leastsq
    cov = covariance matrix
    dof = degrees of freedom (or len(x_data) - len(plsq)) 
    '''
    error =  len(plsq)*[0]
    chisq=sum(fit_info["fvec"]*fit_info["fvec"])
    for i in arange(len(plsq)):
        error[i] = sqrt(cov[i,i])*sqrt(chisq/dof)
    return error
    

