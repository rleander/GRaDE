import numpy as np
from tqdm.notebook import tqdm
import numpy.ma as ma


def all_distributions(distfun, parms, levels):
    distribution = []
    for level in levels:
        distribution.append(distfun(level,parms))
    dist_np=np.array(distribution)
    return dist_np

def linfit(xvec, yvec, wvec=None):
    if wvec is None:
        w = np.ones(len(xvec))
    else:
        w = wvec
    sw = ma.sum(w)
    sx = ma.sum(w*xvec)
    sy = ma.sum(w*yvec)
    sxx = ma.sum(w*xvec*xvec)
    sxy = ma.sum(w*xvec*yvec)   
    a = (sxy-sy*sx/sw)/(sxx-sx*sx/sw)
    b = (sy-a*sx)/sw
    return a,b

def outintegrate_equ(parms, levels, distfun):
    # based on a sample average
    # this works if we have an evenly distributed sample (in terms of reference q) in parms
    # each record in parms consists of a tuple of parameters
    # does not have to be ordered
    nx = len(parms['distpar'])
    nl = len(levels)
    cdf = np.array([0.0]*nl)
#   for i in range(nx):
    for i in tqdm(list(range(nx))):
        distribution = all_distributions(distfun, parms['distpar'][i], levels)
        cdf += distribution
    cdf /= nx
    return cdf

def outintegrate(parms, levels, distfun):
    # based on a numerical integral over q probabilities
    # this works if we have an ordered set (in terms of reference q) in parms, each with its non-exceedance probability
    # each record in parms consists of a tuple of parameters AND an associated non-exceedance probability
    # has to be ordered on reference Q (that also means ordered by non-exceedance probability) 
    nx = len(parms['F'])
    nl = len(levels)
    cdf = np.array([0.0]*nl)
    H = parms['F']
    wtsum = 0.0
#   for i in range(nx):
    for i in tqdm(list(range(nx))):
        if i==0:
            wt = 0.5*(H[0]+H[1])
        elif i==nx-1:
            wt = 1.-0.5*(H[nx-1]+H[nx-2])
        else:
            wt = 0.5*(H[i+1]-H[i-1])
        wtsum = wtsum + wt
        distribution = all_distributions(distfun, parms['distpar'][i], levels)
        cdf += distribution * wt
    return cdf


# outintegration with log-linear extrapolation
def outintegrateExtrapLoglin(parms, levels, distfun, nlevels, npoints):
    # extrapolates the discharge curve beyond the highest point with an extra nlevels
    # the additional points are equally spaced on the F-axis between the highest F and 1
    # the additional parametersets are obtained by loglinear extrapolation from the 
    # highest npoints points in the data on the F-axis
    nparam = len(parms['distpar'][0])
    list_a = []
    list_b = []

    # fit the linear coefficients for each of the parameters
    for ipar in range(nparam):
        y = np.array([partuple[ipar] for partuple in parms['distpar']][-npoints:])
        F = np.array(parms['F'][-npoints:])
        logt = np.log(1./(1.-F))
        a,b = linfit(logt,y)
        list_a.append(a)
        list_b.append(b)

    # add a tuple of parameters and the F-value for points equally spaced between F_last and 1.0
    fhighest = np.max(F)
    for ipt in range(nlevels):
        ff_ipt = fhighest + (1.-fhighest)/(nlevels+1)*(ipt+1)
        logt_ipt = np.log(1./(1.-ff_ipt))
        parms['distpar'].append(tuple([list_a[ipar]*logt_ipt+list_b[ipar] for ipar in range(nparam)]))
        parms['F'].append(ff_ipt)
    cdf = outintegrate(parms, levels, distfun)
    return cdf
 
