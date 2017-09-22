from numpy import *

# own stuff
from analysis.lib.fitting import fit as fit

### gaussian esr dips in the PL
def fit_ESR_gauss(g_a, g_A, g_sigma, g_x0, *arg):
    """
    fitfunction for gaussian esr dips,
        y(x) = a - |A|*sum_i(exp(-((x-x0_i)/sigma)**2))

    Initial guesses (in this order):
        g_a : offset
        g_A : dip depth (all the same)
        g_sigma : std dev of the dip width (all the same)
        g_x0 : zero-point of the esr pattern
    
    For the splittings:
        an arbitrary amount of splits, in the form
        (m, g_s) : where m=multiplicity, g_s=initial guess for the
        splitting of two neighboring dips

    Example:
        giving (2, 2e3), (3, 2e3), would describe two splittings,
        e.g., the two-fold one for the zeeman splitting of the +/-1
        manifolds, and on top of that the 3-fold splitting of the N14.
        This results in 6 peaks in total, where 4 of them appear,
        because of the equal splitting in the two cases, on top
        of each other, yielding 2 low, and two deep dips.
    """
    fitfunc_str = 'a - |A|*sum_i(exp(-((x-x0_i)/sigma)**2))'
    
    no_splits = len(arg)

    a = fit.Parameter(g_a, 'a')
    A = fit.Parameter(g_A, 'A')
    sigma = fit.Parameter(g_sigma, 'sigma')
    x0 = fit.Parameter(g_x0, 'x0')
    p0 = [a, A, sigma, x0]
    g_p0 = [g_a, g_A, g_sigma, g_x0]

    # print 'fitting with %d splittings' % no_splits

    splits = []
    for i, s in enumerate(arg):
        fitfunc_str += '\nsplitting s%d with multiplicity %d' % (i, s[0])
        splits.append(fit.Parameter(s[1], 's%d'%i))
        p0.append(splits[i])

    # remove the fixed params from the p0 array
    # fixedp = []
    # for i,p in enumerate(p0):
    #     if i in fixedsplits:
    #         fixedp.append(p)
    # for p in fixedp:
    #     p0.remove(p)
    
  
    def fitfunc(x):
        pts = [x0()]
        pts_next = []
        
        # generate the points where the dips sit for the current param
        # values
        for i in range(no_splits):
            m = arg[i][0]
            s = arg[i][1]
            for pt in pts:
                j = 0
                while j < m:
                    split = splits[i]()
                    pts_next.append(pt-split*(m-1.)/2. + j*split)
                    j+=1

            pts = pts_next
            pts_next = []

        depth = 0.
        for p in pts:
            # it is very convenient to force the dip to be negative. I changed this allready quite a few times.
            # please leave it in!!! -Machiel apr-2014
            depth += abs(A()) * exp(-((x-p)/sigma())**2)

        return a() - depth

    return p0, fitfunc, fitfunc_str

