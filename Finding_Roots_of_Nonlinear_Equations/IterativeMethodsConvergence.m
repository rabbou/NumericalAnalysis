%%
clear all;
close all;
clc;

%% Finding 3^(1/3) using Fixed point iterations, Newton's method and a Cubic Method
p0 = 1.3;
N0 = 20;
f = @(x) x^3-3;
df = @(x) 3*x^2;
ddf = @(x) 6*x;

g_fp = @(x) x-.1*f(x); % fixed point update
g_nm = @(x) x-f(x)/df(x); % Newton Method update
g_cm = @(p0) p0 - f(p0)/df(p0) - ddf(p0)/(2*df(p0))*(f(p0)/df(p0))^2; % Cubic Method Update

% set the true solution as the result of the cubic method run for 20 iterations
[~, true_p] = e(p0, N0, g_cm, 0);

fixed_point = e(p0, N0, g_fp, true_p);
newt_met = e(p0, N0, g_nm, true_p);
cub_met = e(p0, N0, g_cm, true_p);

fig = figure(1);
a1 = semilogy(1:N0,fixed_point(1:N0), 'linewidth', 1.5);
t1 = "Fixed Point Iteration";
hold on
a2 = semilogy(1:N0,newt_met(1:N0), 'linewidth', 1.5);
t2 = "Newton's Method";
hold on
a3 = semilogy(1:N0,cub_met(1:N0), 'linewidth', 1.5);
t3 = "Cubic Method";
legend([a1;a2;a3], t1,t2,t3,'FontSize',12,'interpreter','latex');
set(gca,'YMinorTick','off');

xlabel('Iteration number','interpreter','latex','FontSize',15);
ylabel('Absolute error','interpreter','latex','FontSize',15);
title('Error Convergence in approximating $\sqrt[3]{3}$','interpreter','latex','FontSize',15);

saveas(fig, 'CubeRootThreeNewton.jpg');

%% Function for iterative methods for finding roots.
% It takes a given update (Fixed point, Newton, etc.) and a true solution
% and returns the approximation and an array of error for each iteration
function [error_terms, p] = e(p0, N0, g, true_p)
error_terms = zeros(1,N0);
Iter = 1;
p = p0;
error_terms(Iter) = abs(p - true_p);
while Iter < N0
    Iter = Iter+1;
    p = g(p);
    error_terms(Iter) = abs(p-true_p);
end
end