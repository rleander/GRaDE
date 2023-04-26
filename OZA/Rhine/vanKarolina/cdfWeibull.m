function res = cdfWeibull(x,a,k,c)

res = 1-exp(-(max((x-c)./a,0.0).^k));

