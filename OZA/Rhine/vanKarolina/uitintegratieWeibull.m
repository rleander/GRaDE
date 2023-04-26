clc
clear all
close all

# [inNum,~] = xlsread('Input\weibull\weibull_new.xlsx','werklijn');
inNum = dlmread('input/werklijn.txt',',',1,0)
Q         = inNum(:,3);
T         = inNum(:,2);
prob      = 1./T;

#[inNum,~] = xlsread('Input\weibull\weibull_new.xlsx','par');
inNum = dlmread('input/weibull_parameters_uniform.txt',',',1,0)
Q_par     = inNum(:,1);
a_par     = inNum(:,2);
k_par     = inNum(:,3);
c_par     = inNum(:,4);

delta_x  = 10;
delta_q  = 10;
x_grid = 200:delta_x:20000;
q_grid = 200:delta_q:20000;

a_par_res = interp1(Q_par,a_par,q_grid,'linear','extrap');
k_par_res = interp1(Q_par,k_par,q_grid,'linear','extrap');
c_par_res = interp1(Q_par,c_par,q_grid,'linear','extrap');

for i = 1:length(q_grid)
    int(i) = sum((1-cdfWeibull(x_grid(i)-q_grid,a_par_res,k_par_res,c_par_res)).*pdfAfvoer(Q,prob,q_grid)*delta_x);
    

    
end

figure
semilogx(T,Q);
hold on
semilogx(1./int,x_grid,'+');

# Laad het 'python' resultaat
inCheck = dlmread('check/uitgeintegreerd.txt',',',1,0)
Tr_check = inCheck(:,2);
Q_check = inCheck(:,3);
semilogx(Tr_check,Q_check,'g--');

grid on
