import sys
import numpy as np
import scipy.special as scs
import scipy.optimize as op

def qfunbeta(alpha,beta,a,b,P):
    return scs.betaincinv(alpha,beta,P)*(b-a)+a

def quantilesetbeta(n,alpha,beta,a,b):
    pvalues = [(i+1.0-0.3)/(n+0.4) for i in range(n)]
    qvalues = [qfunbeta(alpha,beta,a,b,p) for p in pvalues]
    return pvalues, qvalues

def quantilesetnormal(n,mu,sigma):
    pvalues = [(i+1.0-0.3)/(n+0.4) for i in range(n)]
    qvalues = [mu+sigma*(2.**0.5)*scs.erfinv(2.*p-1.) for p in pvalues]
    return pvalues, qvalues

def pdfbeta(n,alpha,beta,a,b,bnd):
    xvalues = [bnd[0]+(bnd[1]-bnd[0])/(n-1)*i for i in range(n)]
    qvalues = [((x-a)/(b-a))**(alpha-1 )*(((b-x)/(b-a))**(beta-1))/(b-a)/scs.beta(alpha,beta)   for x in xvalues]
    return xvalues, qvalues

#define distribution function and quantile function for the Beta distribution
def beta_forward(alpha,beta,a,b,y):
    if (y<a):
        return 0.0
    elif (y>b):
        return 1.0
    else:
        return(scs.betainc(alpha,beta,(y-a)/(b-a)))

def beta_reverse(alpha,beta,a,b,p):
    return(scs.betaincinv(alpha,beta,p)*(b-a)+a)

def P_from_Beta(beta):
    return (0.5+0.5*scs.erf(beta/(2**0.5)))
def Q_from_Beta(beta):
    return (0.5+0.5*scs.erf(-beta/(2**0.5)))
def Beta_from_Q(Q):
    return(-scs.erfinv(2*Q - 1.0)*(2**0.5))
def Beta_from_P(Q):
    return(scs.erfinv(2*Q - 1.0)*(2**0.5))

def BetaMoments(alpha,beta,A,B):
    mean = alpha/(alpha+beta)*(B-A)+A
    variance = alpha*beta/(alpha+beta)**2./(alpha+beta+1)*(B-A)
    skewness = 2.*(beta-alpha)*(alpha+beta+1)**0.5 \
             /((alpha+beta+2)*(alpha*beta)**0.5)
    exkurt = 6*((alpha-beta)**2.*(alpha+beta+1)-(alpha*beta)*(alpha+beta+2)) \
             /(alpha*beta*(alpha+beta+2)*(alpha+beta+3))
    return (mean, variance, skewness, exkurt)


def quantileplot(ensembles, results, lbl):
    q_obs = sorted(ensembles[lbl])
    n_obs = len(q_obs)
    result = results[lookup[lbl]]
    alpha = result['estimate'][0]
    beta = result['estimate'][1]
    a = result['estimate'][2]
    b = result['estimate'][3]
    pvalues, q_fit = quantilesetbeta(n_obs,alpha,beta,a,b)
    mysample = np.array(quantilesetbeta(n_obs,alpha,beta,a,b))             # exact quantiles (order statistics) for the given parameters
    plt.plot(pvalues,q_obs,'+',pvalues,q_fit,'r-')

def pdfplot_duo(ax, results_m, lbl):
    result = results_m[lbl]
    n_obs=1000
    alpha = result['estimate'][0]
    beta = result['estimate'][1]
    a = result['estimate'][2]
    b = result['estimate'][3]
    xvalues, q_pdf_m = pdfbeta(n_obs,alpha,beta,a,b)

    alpha = result['estimate2'][0]
    beta = result['estimate2'][1]
    a = result['estimate2'][2]
    b = result['estimate2'][3]
    xvalues, q_pdf_q = pdfbeta(n_obs,alpha,beta,a,b)
    ax.plot(xvalues,q_pdf_m,'r-',xvalues,q_pdf_q,'g-')        
    
def quantileplot_duo(ax, ensembles, results_m, lbl):
    q_obs = sorted(ensembles[lbl])
    n_obs = len(q_obs)
    result = results_m[lbl]
    alpha = result['estimate'][0]
    beta = result['estimate'][1]
    a = result['estimate'][2]
    b = result['estimate'][3]
    pvalues, q_fit_m = quantilesetbeta(n_obs,alpha,beta,a,b)

    alpha = result['estimate2'][0]
    beta = result['estimate2'][1]
    a = result['estimate2'][2]
    b = result['estimate2'][3]
    pvalues, q_fit_q = quantilesetbeta(n_obs,alpha,beta,a,b)
    
    mysample = np.array(quantilesetbeta(n_obs,alpha,beta,a,b))             # exact quantiles (order statistics) for the given parameters
    ax.plot(pvalues,q_obs,'b+',label="Sample")    
    ax.plot(pvalues,q_fit_m,'r-',label="Fit Moments")    
    ax.plot(pvalues,q_fit_q,'g-',label="Fit Quantile")    
    ax.set_xlabel('P-value [-]')
    ax.set_ylabel('Percentile [m3 s-1]')    
    ax.legend()
    ax.set_title("Fits, percentiles versus P-value")



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

def quantilesetbeta(n,alpha,beta,a,b):
    return [qfunbeta(alpha,beta,a,b,(i+1.0-0.3)/(n+0.4)) for i in range(n)]

def quantilesetbeta(n,alpha,beta,a,b):
    pvalues = [(i+1.0-0.3)/(n+0.4) for i in range(n)]
    qvalues = [qfunbeta(alpha,beta,a,b,p) for p in pvalues]
    return pvalues, qvalues


def qfunbeta(alpha,beta,a,b,P):
    return scs.betaincinv(alpha,beta,P)*(b-a)+a

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
    if 'beta' in kwargs:
        q = float(kwargs['beta'])
    if 'a' in kwargs:
        a = float(kwargs['min'])
    if 'b' in kwargs:
        b = float(kwargs['max'])
    
    M1, M2, A3, A4 = meanVarSkewKurt(sample)
    if (alpha_fixed and beta_fixed):
        r = p + q
        print(r, "<----")
    else:
        r = 6 * (A4-A3**2.0-1.0) / (6.0+3*A3**2.0-2*A4)
        if alpha_fixed:
            q = r - p
        elif beta_fixed:
            p = r - q
        else:
            p = r/2.0 * (1.0 - (1.0 - 24.0*(r+1)/(A4*(r+2)*(r+3) - 3*(r+1)*(r-6)))**0.5)
            q = r - p
            if A3>0.0:		# positively skewed
                q,p = max(p,q),min(p,q)
            else:
                p,q = max(p,q),min(p,q)

    if not(a_fixed and b_fixed):
        range = (M2 * r**2.*(r+1)/(p*q))**0.5
        if b_fixed:
            a = b - range
        else:
            if not a_fixed:
                a = M1 - p/r*range
            b = a + range
        return np.array([p,q,a,b])


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
    est2 = np.array([a,b,AA,BB])
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

    




    

 

