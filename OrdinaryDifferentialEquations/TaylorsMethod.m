%%
clear all;
close all;
clc;

%% Taylor's Method
% Use Taylor?s method of order two with h = 0.1 to approximate the solution
% for y' = 2/t*y+t^2*exp(t)
w = taylor(.1,1,2,0,@(t, y) 2/t*y+t^2*exp(t)+.1/2*(2*y/t^2+4*t*exp(t)+t^2*exp(t)));
fprintf("\n")

% Use Taylor?s method of order two with h = 0.1 to approximate the solution
% for y' = 2/t*y+t^2*exp(t)
w = taylor(.1,1,2,0,@(t, y) 2/t*y+t^2*exp(t) + .1*(y/t^2+2*t*exp(t)+t^2*exp(t)/2) ...
    + (0.1^2/6)*exp(t)*(6+6*t+t^2) + (0.1^2/24)*exp(t)*(12+8*t+t^2));
fprintf("\n")

%%
function w = taylor(h, a, b, alpha, T)
N = (b-a)/h;
t = a;
w = alpha;
for i = 1:N
    w = w + h * T(t, w);
    t = a+i*h;
    act = t^2*(exp(t)-exp(1));
    error = abs(act-w);
    fprintf("t = %.1f; w = %.4f; Actual y(%.1f) = %.4f; Error: %.4f\n", t, w, t, act, error)
end
end