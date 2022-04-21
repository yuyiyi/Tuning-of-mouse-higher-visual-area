function [f,df,ddf] = squarefun(x)
%  [f,df,ddf] = expfun(x)
%
%  replacement for 'exp' that returns 3 arguments (value, 1st & 2nd deriv)

f = x.^2;
df = 2*f;
ddf = 2;
