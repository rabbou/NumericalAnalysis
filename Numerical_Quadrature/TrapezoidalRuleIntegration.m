%%
clear all;
close all;
clc;

%% Approximate erf function using Trapezoidal Rule
xx = linspace(0,5,1000);

fig = figure('Position', [0,0,900,400]);
subplot(1,2,1);
approx = arrayfun(@(z) erf0(z, 20), xx);
plot(xx, approx, 'color', '#00CED1', 'Linewidth', 1.1);
xlabel('$x$','interpreter','latex','FontSize',15)
ylabel('$erf(x)$','interpreter','latex','FontSize',15)
title('Approximation of ${erf}(x)$', ...
    'interpreter','latex','FontSize',15)

subplot(1,2,2);
plot(xx, abs(erf(xx) - approx), 'color', '#ffa500', 'Linewidth', 1.1);
xlabel('$x$','interpreter','latex','FontSize',15)
ylabel('$|erf_{0}(x)-{erf}(x)|$','interpreter','latex','FontSize',15)
title('Pointwise error for $erf(x)$', ...
    'interpreter','latex','FontSize',15)
saveas(fig, 'ErfTrapezoidalApprox.jpg')

%%
function p = erf0(x, n)
f = @(x) exp(-x.^2)*2/sqrt(pi);
h = x/n;
xj = h:h:x-h;
p = h/2*(f(0)+ 2*sum(f(xj))+f(x));
end