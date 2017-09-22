from numpy import *

# own stuff
from analysis.lib.fitting import fit

### gaussian decay (FID, yielding T2*)
def fit_FID_gauss(g_tau, g_A, g_a, *arg):
    """
    fitfunction for a gaussian decay,
        y(x) = a + A*exp(-(x/tau)**2)

    Initial guesses (in this order):
        g_tau : decay constant
        g_A : amplitude
        g_a : offset
    """
    
    tau = fit.Parameter(g_tau, 'tau')
    A = fit.Parameter(g_A, 'A')
    a = fit.Parameter(g_a, 'a')
    p0 = [tau, A, a]
    def fitfunc(x): return a() + A()*exp(-(x/tau())**2)
    return p0, fitfunc, fitfunc_str

def fit_ramsey_gaussian_decay(g_tau, g_a, *arg):
    """
    fitfunction for a ramsey modulation, with gaussian decay,
        y(x) = a + A*exp(-(x/tau)**2) * mod,

        where:
        mod = sum_i(cos(2pi*f_i*x +phi) - 1)

    Initial guesses (in this order):
        g_tau : decay const
        g_A : Amplitude
        g_a : offset

        For the modulation:
        an arbitrary no of tuples, in the form
        (g_f, g_A)[i] = (frequency, Amplitude, phase)[i]
    """
    fitfunc_str = 'y(x) = a + exp(-(x/tau)**2)*('
    no_frqs = len(arg)
    if no_frqs == 0:
        print 'no modulation frqs supplied'
        return False
    
    tau = fit.Parameter(g_tau, 'tau')
    # A = fit.Parameter(g_A, 'A')
    a = fit.Parameter(g_a, 'a')
    p0 = [tau, a]

    print 'fitting with %d modulation frequencies' % no_frqs

    frqs = []
    amplitudes = []
    phases = []
    print 'arg2 =' + str(arg)
    for i, m in enumerate(arg):
        fitfunc_str += 'A%d*cos(2pi*f%d*x+phi%d)' % (i, i, i)
        frqs.append(fit.Parameter(m[0], 'f%d'%i))
        phases.append(fit.Parameter(m[2], 'phi%d'%i))
        amplitudes.append(fit.Parameter(m[1], 'A%d'%i))
        p0.append(frqs[i])
        p0.append(amplitudes[i])
        p0.append(phases[i])
    fitfunc_str += ')'

    def fitfunc(x):
        prd = exp(-(x/tau())**2)
        mod = 0
        for i in range(no_frqs):
            mod += amplitudes[i]() * (cos(2*pi*frqs[i]()*x+phases[i]()))
        return a() + prd*mod

    return p0, fitfunc, fitfunc_str

def fit_ramsey_hyperfinelines_fixed(g_tau, g_A, g_a,g_det,g_hf,g_phi1,g_phi2,g_phi3, *arg):
    """
    fitfunction for a gaussian decay,
        y(x) = a + A*exp(-(x/tau)**2)

    Initial guesses (in this order):
        g_tau : decay constant
        g_A : amplitude
        g_a : offset
    """
    
    tau = fit.Parameter(g_tau, 'tau')
    A = fit.Parameter(g_A, 'A')
    a = fit.Parameter(g_a, 'a')
    det = fit.Parameter(g_det, 'det')
    hf = fit.Parameter(g_hf, 'hf')
    phi1 = fit.Parameter(g_phi1, 'phi1')
    phi2 = fit.Parameter(g_phi2, 'phi2')
    phi3 = fit.Parameter(g_phi3, 'phi3')
    p0 = [tau, A, a,det,hf,phi1,phi2,phi3]
    fitfunc_str = 'sumf of three cos and decay'

    def fitfunc(x): return a() + A()*exp(-(x/tau())**2) * (cos(2*pi*det()*x+phi1())+cos(2*pi*(det()-hf())*x+phi2())+cos(2*pi*(det()+hf())*x+phi3()))/3.
    return p0, fitfunc, fitfunc_str






def fit_ramsey_hyperfinelines_fixed_6cos(g_tau, g_A, g_a,g_det,g_hf_N,g_hf_C):
    """
    fitfunction for a gaussian decay,
        y(x) = a + A*exp(-(x/tau)**2)

    Initial guesses (in this order):
        g_tau : decay constant
        g_A : amplitude
        g_a : offset
    """
    
    tau = fit.Parameter(g_tau, 'tau')
    A = fit.Parameter(g_A, 'A')
    a = fit.Parameter(g_a, 'a')
    det = fit.Parameter(g_det, 'det')
    hf_N = fit.Parameter(g_hf_N, 'hf_N')
    hf_C = fit.Parameter(g_hf_C, 'hf_C')
    # phi1 = fit.Parameter(g_phi1, 'phi1')
    # phi2 = fit.Parameter(g_phi2, 'phi2')
    # phi3 = fit.Parameter(g_phi3, 'phi3')
    p0 = [tau,A,a,det,hf_N,hf_C]
    fitfunc_str = 'sumf of six cos and decay'

    def fitfunc(x): 
        return a() + A()*exp(-(x/tau())**2) * (cos(2*pi*(det()+hf_C())*x)+(cos(2*pi*(det()-hf_C())*x)+cos(2*pi*(det()+hf_N()+hf_C())*x)+cos(2*pi*(det()+hf_N()-hf_C())*x)+cos(2*pi*(-det()+hf_N()-hf_C())*x)+cos(2*pi*(-det()+hf_N()+hf_C())*x)))/6.


    return p0, fitfunc, fitfunc_str

def fit_ramsey_hyperfinelines_fixed_6sin(g_tau, g_A, g_a,g_det,g_hf_N,g_hf_C):
    """
    fitfunction for a gaussian decay,
        y(x) = a + A*exp(-(x/tau)**2)

    Initial guesses (in this order):
        g_tau : decay constant
        g_A : amplitude
        g_a : offset
    """
    
    tau = fit.Parameter(g_tau, 'tau')
    A = fit.Parameter(g_A, 'A')
    a = fit.Parameter(g_a, 'a')
    det = fit.Parameter(g_det, 'det')
    hf_N = fit.Parameter(g_hf_N, 'hf_N')
    hf_C = fit.Parameter(g_hf_C, 'hf_C')
    # phi1 = fit.Parameter(g_phi1, 'phi1')
    # phi2 = fit.Parameter(g_phi2, 'phi2')
    # phi3 = fit.Parameter(g_phi3, 'phi3')
    p0 = [tau,A,a,det,hf_N,hf_C]
    fitfunc_str = 'sumf of six sin and decay'

    def fitfunc(x): 
        return a() + A()*exp(-(x/tau())**2) * (sin(2*pi*(det()+hf_C())*x)+(sin(2*pi*(det()-hf_C())*x)+sin(2*pi*(det()+hf_N()+hf_C())*x)+sin(2*pi*(det()+hf_N()-hf_C())*x)+sin(2*pi*(det()-hf_N()+hf_C())*x)+sin(2*pi*(det()-hf_N()-hf_C())*x)))/6.


    return p0, fitfunc, fitfunc_str

def fit_ramsey_14N_fixed_13C_opt(g_tau, g_A, g_a,g_det,g_hf_N,g_phi1,g_phi2,g_phi3, *arg):
    """
    fitfunction for a gaussian decay,
        y(x) = a + A*exp(-(x/tau)**2)* \sum_n cos(2*pi*f_n x)

    Initial guesses (in this order):
        g_tau : decay constant
        g_A : amplitude
        g_a : offset

    The splitting between 14N hyperfine lines is constrained to be equal.
    The amplitude of all hyperfine lines are the same and the phases are set to zero.
    13C hyperfine splitting can be passed in arg. 
    """
    
    tau = fit.Parameter(g_tau, 'tau')
    A = fit.Parameter(g_A, 'A')
    a = fit.Parameter(g_a, 'a')
    det = fit.Parameter(g_det, 'det')
    hf_N = fit.Parameter(g_hf_N, 'hf_N')
    phi1 = fit.Parameter(g_phi1, 'phi1')
    phi2 = fit.Parameter(g_phi2, 'phi2')
    phi3 = fit.Parameter(g_phi3, 'phi3')
    p0 = [tau, A, a,det,hf_N,phi1,phi2,phi3]
    fitfunc_str = 'a+A*exp(-(x/tau)**2)*(cos(2*pi*det)+cos(2*pi*(det+hf_N))+cos(2*pi*(det-hf_N))'


    hf_C = []
    print arg
    for i, m in enumerate(arg):
        fitfunc_str += 'cos(2pi*f%d*x)' % (i)
        hf_C.append(fit.Parameter(m, 'hf_C%d'%i))
        p0.append(hf_C[i])
    fitfunc_str += ')'

    no_hf_C = len(arg)
    print "number 13C split . ", no_hf_C
    hf_str = ''
    #print 'Fit with {:.i} 13C nuclear spin splitting'.format(no_hf_C)
    if no_hf_C !=0:
        fitfunc_str = 'a+A*exp(-(x/tau)**2)*(\sum_i cos(2*pi*(hf_i+det))\n, where hf_i =m_N*hf_N'
        for i in range(no_hf_C):
            fitfunc_str += ' + m_C{:g}*hf_C{:g}'.format(i,i)
            hf_str += ', m_C \in {+/-1/2}'
        fitfunc_str += '\nand m_N \in {-1,0,+1}'+ hf_str



    def fitfunc(x): 
        #return a() + A()*exp(-(x/tau())**2) * (cos(2*pi*det()*x+cos(2*pi*(det()-hf_N())*x+cos(2*pi*(det()+hf_N())*x)/3.
        prd = exp(-(x/tau())**2)
        mod = 0
        phi=[phi1(),phi2(),phi3()]
        if no_hf_C == 0 :
            for j in range(-1,2,1):
                mod += A() * (cos(2*pi*(det()+j*hf_N())*x))
        else:
            for i in range(no_hf_C):
                for j in range(-1,2,1):
                    mod += A() * (cos(2*pi*(det()+j*hf_N() +1/2.*hf_C[i]())*x+phi[j+1]))
                    mod += A() * (cos(2*pi*(det()+j*hf_N() -1/2.*hf_C[i]())*x+phi[j+1]))
        return a() + prd*mod

    return p0, fitfunc, fitfunc_str

def fit_gaussian_decaying_cos_withoffset(g_avg, g_A, g_t, g_a, g_f, g_phi, g_b):
    '''
    Fits cosine with offset modulated by Gaussian decay
    Useful to fit e.g. electron Ramsey without nitrogen MBI
    '''
    fitfunc_str = 'A *exp(-(x/t)**2) * (a * cos(2pi * (f*x + phi/360) ) + b)'

    avg = fit.Parameter(g_avg, 'avg')
    A = fit.Parameter(g_A, 'A')
    t = fit.Parameter(g_t, 't')
    a = fit.Parameter(g_a, 'a')
    f = fit.Parameter(g_f, 'f')
    phi = fit.Parameter(g_phi, 'phi')
    b = fit.Parameter(g_b, 'b')
    
    print 'guessed frequency is '+str(g_f)
    
    p0 = [avg, A, t, a, f, phi, b]

    def fitfunc(x):
        return avg() + A()*exp(-(x/t())**2) * ( a() * cos(2*pi*( f()*x + phi()/360.)) + b() )

    return p0, fitfunc, fitfunc_str