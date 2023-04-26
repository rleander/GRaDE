import numpy as np
import scipy.optimize as op
import scipy.special as scs

def qfunweibull(parms,P):
    a, k, c = parms[:3]
    q = 1. - P
    # De weibull inverse implementatie uit ProbLib
    weibullInverse = ( -np.log( q ) )**( 1. / k ) * a + c
    return weibullInverse

def quantilesetweibull(n,parms):
    pvalues = [(i+1.0-0.3)/(n+0.4) for i in range(n)]
    qvalues = [qfunweibull(parms,p) for p in pvalues]
    return pvalues, qvalues

def pdfweibull(xx,parms):
    a, k, c = parms[:3]
    if (xx-c)/a>0:
       pdf = k/a*((xx-c)/a)**(k-1)*np.exp(-((xx-c)/a)**k)
    else:
       pdf = 0.0
    return pdf

def cdfweibull(xx,parms):
    a, k, c = parms[:3]
    if (xx-c)/a>0:
       cdf = 1.0 - np.exp(-((xx-c)/a)**k)
    else:
       cdf = 0.0
    return cdf


def LskewWeibull(k):
    return (1.-3./(2.**(1./k))+2./(3.**(1./k)))/(1.-1./(2.**(1./k)))

def fitWeibullLMoments(sample):
    MAXITER = 10000
    x = np.sort(sample)
    beta0 = 0.0
    beta1 = 0.0
    beta2 = 0.0
    n = 0
    for j in range(len(x)):
        n = n + 1
        beta0 = beta0 + x[j]
        beta1 = beta1 + x[j] * j
        beta2 = beta2 + x[j] * j * (j-1)
    beta0 = beta0/n
    beta1 = beta1/n/(n-1)
    beta2 = beta2/n/(n-1)/(n-2)
        
    L1 = beta0
    L2 = 2*beta1 - beta0
    L3 = 6*beta2 - 6*beta1 + beta0
    Lskew = L3/L2
    
    # Solve equation for Lskew to obtain k
    tol = 0.001
    k0 = 0.1
    k1 = 20
    Lskew0 = LskewWeibull(k0)
    Lskew1 = LskewWeibull(k1)
    for iter in range(MAXITER):
        k = (Lskew-Lskew0)/(Lskew1-Lskew0)*(k1-k0)+k0
        k = (k0+k1)/2
        Lskew2 = LskewWeibull(k)
        if (Lskew2-Lskew)*(Lskew0-Lskew)>0:
            k0 = k    # if Lskew2 and Lskew0 are on the same side of Lskew, replace
        else:
            k1 = k
        if (k1-k0) < tol:
            break
    a = L2/scs.gamma(1./k+1.)/(1.-1./(2.**(1./k)))         
    c = L1 - a*scs.gamma(1./k+1.)
    return a,k,c



    




    

 

