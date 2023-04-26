clc
clear all
close all

inNum = dlmread('input/werklijn.txt',',',1,0)
Q         = inNum(:,3);
T         = inNum(:,2);
prob      = 1./T;

inNum = dlmread('input/beta_parameters_uniform.txt',',',1,0)
Q_par     = inNum(:,1);
a_par     = inNum(:,2);
b_par     = inNum(:,3);
alpha_par = inNum(:,4);
beta_par  = inNum(:,5);

delta_x  = 10;
delta_q  = 10;
x_grid = 200:delta_x:20000;
q_grid = 200:delta_q:20000;

a_par_res = interp1(Q_par,a_par,x_grid,'linear','extrap');
b_par_res = interp1(Q_par,b_par,x_grid,'linear','extrap');
alpha_par_res = max(interp1(Q_par,alpha_par,x_grid,'linear','extrap'),0.001);
beta_par_res = max(interp1(Q_par,beta_par,x_grid,'linear','extrap'),0.001);

for i = 1:length(q_grid)
    int(i) = sum((1-cdfBeta(x_grid(i)-q_grid,a_par_res,b_par_res,alpha_par_res,beta_par_res)).*pdfAfvoer(Q,prob,q_grid)*delta_x);
    

    
end

figure
semilogx(T,Q);
#semilogx(-1./(log(1.-1./T)),Q);
hold on
semilogx(1./int,x_grid,'+');
#semilogx(-1./(log(1.-int)),x_grid,'+');

# Laad het 'python' resultaat
inCheck = dlmread('check/uitgeintegreerd_beta.txt',',',1,0);
Tr_check = inCheck(:,2);
Q_check = inCheck(:,3);
semilogx(Tr_check,Q_check,'g--');
#semilogx(-1./(log(1.-1./Tr_check)),Q_check,'g--');

grid on
