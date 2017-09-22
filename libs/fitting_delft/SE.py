from numpy import *

# own stuff
from analysis.lib.fitting import fit

### gaussian decay (FID, yielding T2*)
def fit_echo(g_tau, g_A, g_a, g_k, *arg):
    """
    fitfunction for a gaussian decay,
        y(x) = a + A*exp(-(x/tau)**k)

    Initial guesses (in this order):
        g_tau : decay constant
        g_A : amplitude
        g_a : offset
        g_k : exponent
    """
    fitfunc_str='a() + A()*exp(-(x/tau))**k)'
    tau = fit.Parameter(g_tau, 'tau')
    A = fit.Parameter(g_A, 'A')
    a = fit.Parameter(g_a, 'a')
    k = fit.Parameter(g_k,'k')
    p0 = [tau, A, a,k]
    def fitfunc(x): return a() + A()*exp(-(x/tau())**k())

    return p0, fitfunc, fitfunc_str

