function res = cdfBeta(x,a,b,alpha,beta)
if ((x-a)/(b-a)<=0.0000)
   res = 0.0;
elseif ((x-a)/(b-a)>=1.0000)
   res = 1.0;
else
   res = betainc((x-a)/(b-a),alpha,beta);
endif
