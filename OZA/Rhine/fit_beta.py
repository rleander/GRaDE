import sys
import numpy as np
import scipy.optimize as op
import scipy.special as scs

def qfunbeta(parms,P):
    a, b, alpha, beta = parms[:4]
    return scs.betaincinv(alpha,beta,P)*(b-a)+a

def quantilesetbeta(n,parms):
    a, b, alpha, beta = parms[:4]
    pvalues = [(i+1.0-0.3)/(n+0.4) for i in range(n)]
    qvalues = [qfunbeta(alpha,beta,a,b,p) for p in pvalues]
    return pvalues, qvalues

def pdfbeta(y,parms):
    a, b, alpha, beta = parms[:4]
    pdf = ((y-a)/(b-a))**(alpha-1 )*(((b-y)/(b-a))**(beta-1))/(b-a)/scs.beta(alpha,beta)
    return pdf

def cdfbeta(y,parms):
    a, b, alpha, beta = parms[:4]
    cdf = beta_forward(parms,y)
    return cdf

#define distribution function and quantile function for the Beta distribution
def beta_forward(parms,y):
    a, b, alpha, beta = parms[:4]
    if (y<a):
        return 0.0
    elif (y>b):
        return 1.0
    else:
        return(scs.betainc(alpha,beta,(y-a)/(b-a)))

def beta_reverse(parms,p):
    a, b, alpha, beta = parms[:4]
    return(scs.betaincinv(alpha,beta,p)*(b-a)+a)

def P_from_Beta(beta):
    return (0.5+0.5*scs.erf(beta/(2**0.5)))
def Q_from_Beta(beta):
    return (0.5+0.5*scs.erf(-beta/(2**0.5)))
def Beta_from_Q(Q):
    return(-scs.erfinv(2*Q - 1.0)*(2**0.5))
def Beta_from_P(Q):
    return(scs.erfinv(2*Q - 1.0)*(2**0.5))

def BetaMoments(parms):
    A, B, alpha, beta = parms[:4]
    mean = alpha/(alpha+beta)*(B-A)+A
    variance = alpha*beta/(alpha+beta)**2./(alpha+beta+1)*(B-A)
    skewness = 2.*(beta-alpha)*(alpha+beta+1)**0.5 \
             /((alpha+beta+2)*(alpha*beta)**0.5)
    exkurt = 6*((alpha-beta)**2.*(alpha+beta+1)-(alpha*beta)*(alpha+beta+2)) \
             /(alpha*beta*(alpha+beta+2)*(alpha+beta+3))
    return (mean, variance, skewness, exkurt)




# R.Leander
# 10-12-2020
# Jackknife estimate of beta parameters from a sample moments.
# Unbiased estimate is approximated by applying the jackknife principle (pseudo-values)
# fitBetaMomentsJack(<sample>,[alpha=...,beta=...,min=...,max=...)
# arguments: sample - numpy array of values (not necessarily sorted), compulsary
#            alpha  - fixed alpha parameter, optional
#            beta   - fixed beta parameter, optional
#            min    - fixed underbound, optional
#            max    - fixed upperbound, optional

def quantilesetbeta(n,a,b,alpha,beta):
    return [qfunbeta(alpha,beta,a,b,(i+1.0-0.3)/(n+0.4)) for i in range(n)]

def quantilesetbeta(n,a,b,alpha,beta):
    pvalues = [(i+1.0-0.3)/(n+0.4) for i in range(n)]
    qvalues = [qfunbeta(alpha,beta,a,b,p) for p in pvalues]
    return pvalues, qvalues

def qfunbeta(parms,P):
    return beta_reverse(parms,P)

def jackknife(estimator,sample,**kwargs):
    # obtain unbiased estimate using the jackknife approach
    estimate = estimator(sample)
    m = np.size(estimate)
    pseudo = np.array([])
    subsample = np.ma.array(sample)
    nn = len(sample)
    for i in range(nn):
        subsample.mask = False
        subsample.mask[i] = True
        pseudo_estimate = estimator(subsample,**kwargs)
        pseudo = np.append(pseudo, nn*estimate - (nn-1)*pseudo_estimate)
    pseudomean = np.array([])
    pseudovar  = np.array([])
    for j in range(m):
        select = [j+m*i for i in range(nn)]
        pseudomean = np.append(pseudomean, pseudo[select].mean())
        pseudovar  = np.append(pseudovar, pseudo[select].var())
    return pseudomean, pseudovar


def meanVarSkewKurt(sample):
    # return sample Mean, Variance, Skewness and Kurtosis
    s = [0,0,0,0,0]
    for i in range(5):
        s[i] = (sample**i).sum()
    n = s[0]
    mean = s[1]/s[0]
    variance = (s[2]-(s[1]*s[1]/n))/(n-1)
    skewness = (s[3] - 3.*s[1]*s[2]/n + 2*s[1]*s[1]*s[1]/(n*n)) / (variance**(3./2.)*n)
    kurtosis = (s[4] - 4.*s[1]*s[3]/n + 6*s[1]*s[1]*s[2]/(n*n) - 3*s[1]*s[1]*s[1]*s[1]/(n*n*n)) / (variance**2.*n)
    return (mean,variance,skewness,kurtosis)  


def fitBetaMoments(sample,**kwargs):
    alpha_fixed = False
    beta_fixed = False
    a_fixed = False
    b_fixed = False
    if 'alpha' in kwargs:
        p = float(kwargs['alpha'])
        alpha_fixed = True
    if 'beta' in kwargs:
        q = float(kwargs['beta'])
        beta_fixed = True
    if 'a' in kwargs:
        a = float(kwargs['a'])
        a_fixed = True
    if 'b' in kwargs:
        b = float(kwargs['b'])
        b_fixed = True
    
    M1, M2, A3, A4 = meanVarSkewKurt(sample)
    if (a_fixed and b_fixed):       # fit a 2-parameter beta on mean and variance
        positive = (M1>(a+b)/2.)
        M1 = (M1 - a)/(b - a)
        M2 = M2/(b - a)**2
        rho = (1.-M1)/M1
        q = (rho-M2*(1.+rho)**2.)/(M2*(1.+rho)**3.)
        p = rho * q
        if positive:		# positively skewed
            p,q = max(p,q),min(p,q)
        else:
            q,p = max(p,q),min(p,q)
    else:                           # fit a 4-parameter beta on mean, variance, skewness and kurtosis
        r = 6 * (A4-A3**2.0-1.0) / (6.0+3*A3**2.0-2*A4)
        p = r/2.0 * (1.0 - (1.0 - 24.0*(r+1)/(A4*(r+2)*(r+3) - 3*(r+1)*(r-6)))**0.5)
        q = r - p
        if A3>0.0:		# positively skewed
            q,p = max(p,q),min(p,q)
        else:
            p,q = max(p,q),min(p,q)

        range = (M2 * r**2.*(r+1)/(p*q))**0.5
        a = M1 - p/r*range
        b = a + range
    return np.array([a,b,p,q])


def fitBetaMomentsJack(sample,**kwargs):
    return jackknife(fitBetaMoments,sample,**kwargs)


def beta_quantiles_residuals(par,**kwargs):
    xs = np.sort(kwargs['qval'])
    fs, xx = quantilesetbeta(len(xs),par[0],par[1],par[2],par[3])
    return xs-xx
    

def fitBetaQuantiles(sample):
    AA = np.min(sample)
    BB = np.max(sample)
    a = 2.0             # start with parabola distribution covering the total interval
    b = 2.0
    est2 = np.array([AA,BB,a,b])
    op_result = op.least_squares(beta_quantiles_residuals,est2,kwargs={'qval':np.array(sample)})
    est2 = op_result.x
    return est2


def main():
    # testcode
    betaparstr = sys.argv[1] 	# alpha, beta, a, b
    alpha, beta, a, b = [float(s) for s in sys.argv[1].split()]
    samplesize = int(sys.argv[2])
#   mysample = np.random.beta(alpha,beta,size=samplesize)*(b-a)+a		# random sample for the given parameter
    mysample = np.array(quantilesetbeta(samplesize,alpha,beta,a,b))             # exact quantiles (order statistics) for the given parameters
    est,var = fitBetaMomentsJack(mysample)

    print ("alpha = %15.5e +/- %12.4e" % (est[0],var[0]**0.5))
    print ("beta  = %15.5e +/- %12.4e" % (est[1],var[1]**0.5))
    print ("min   = %15.5e +/- %12.4e" % (est[2],var[2]**0.5))
    print ("max   = %15.5e +/- %12.4e" % (est[3],var[3]**0.5))

if __name__ == '__main__':
    main()

    




    

 

