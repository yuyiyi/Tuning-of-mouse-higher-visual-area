function f = logstic_fun(x1)

alpha = 5;
beta = 2;
x0 = 20;

f = alpha./(1+exp(-beta*(x1-x0)));
