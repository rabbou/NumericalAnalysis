%%
clear all;
close all;
clc;

%% Euler's method
% Approximate solution for y' = t*exp(3t)-2y
w = euler(0,1,2,0,@(t, y) t*exp(3*t)-2*y);
fprintf("\n")

% Approximate solution for y' = exp(t-y)
w = euler(0,1,2,1,@(t, y) exp(t-y));
fprintf("\n")

%%
function [w, t] = euler(a, b, N, alpha, f)
h = (b-a)/N;
t = a;
w = alpha;
for i = 1:N
    w = w+h*f(t, w);
    fprintf("t = %f; w = %f\n", t, w)
    t = a+i*h;
end
end
