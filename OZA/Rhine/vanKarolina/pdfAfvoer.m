function res = pdfAfvoer(val,prob,q)

delta     = 10;

prob_q       = exp(interp1(val,log(prob),q,'linear','extrap'));
prob_q_delta = exp(interp1(val,log(prob),q+delta,'linear','extrap'));

prob_q(prob_q>1)             = 1;
prob_q_delta(prob_q_delta>1) = 1;

res = (prob_q-prob_q_delta)./delta;

res(res<0) = 0;