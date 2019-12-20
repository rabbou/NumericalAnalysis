%%
clear all;
close all;
clc;

%% Natural Cubic Spline approximation for cos(5*cos(5*x))

f = @(x) cos(5*cos(5*x));
xj1 = 2*pi.*(0:1:10)/10;
xj2 = 2*pi.*(0:1:100)/100;
xx = linspace(0, 2*pi, 1000);
yy1 = spline(xj1,f(xj1),xx);
yy2 = spline(xj2,f(xj2),xx);
fig = figure('Position', [0,0,800,400]);
plot(xx, f(xx), 'r', 'Linewidth', 1.3)
hold on
plot(xx,yy1,'color', '#136207', 'Linewidth', 1.1)
hold on
plot(xx, yy2,'b')
hold off
legend('$f(x)$','$P(x)$ w/ step size $\frac{2\pi}{10}$', ...
        '$P(x)$ w/ step size $\frac{2\pi}{100}$', ...
        'interpreter','latex','FontSize',12,'Location','northwest')
title('Cubic Spline Approximation of $\cos(5\cos(5x))$','interpreter','latex','FontSize',15)
xlabel('$x$','interpreter','latex','FontSize',15)
ylabel('$y$','interpreter','latex','FontSize',15)
xlim([0 2*pi])
ylim([-1 1.6])

saveas(fig, 'CubicSplineInterpolation.jpg');