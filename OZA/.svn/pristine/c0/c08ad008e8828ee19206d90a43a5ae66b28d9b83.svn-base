import numpy as np
from tqdm.notebook import tqdm

def all_distributions(distfun, parms, levels):
    distribution = []
    for level in levels:
        distribution.append(distfun(level,parms))
    dist_np=np.array(distribution)
    return dist_np

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
