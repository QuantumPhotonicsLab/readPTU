import numpy as np

# own modules
import fit
#=======
#import analysis.lib.fitting.fit as fit
print 'reloaded'

def fit_lorentz(g_a, g_A, g_x0, g_gamma):
    fitfunc_str = 'a + 2*A/np.pi*gamma/(4*(x-x0)**2+gamma**2)'

    a = fit.Parameter(g_a, 'a')
    A = fit.Parameter(g_A, 'A')
    x0 = fit.Parameter(g_x0, 'x0')
    gamma = fit.Parameter(g_gamma, 'gamma')

    p0 = [a, A, x0, gamma]


    def fitfunc(x):
        return a() + 2*A()/np.pi*gamma()/(4*(x-x0())**2+gamma()**2)

    print p0, fitfunc, fitfunc_str


x = np.array([1,3,5])
print type(x)

fit_lorentz(1,2,4,5)
