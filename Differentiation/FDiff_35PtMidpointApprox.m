%%
clear all;
close all;
clc;

%% Approximate f' at a point using the forward difference formula,
% the three-point midpoint formula, and the five-point midpoint formula
% for f(x) = 2sin(5x) + cos(3x) (f'(x) = 10cos(5x) - 3sin(3x)

derivForwardDiff = @(f, h, x) (f(x+h)-f(x))./h;
deriv3PointMid = @(f, h, x) (f(x+h)-f(x-h))./(2*h);
deriv5PointMid = @(f, h, x) (f(x-2*h)-8*f(x-h)+8*f(x+h)-f(x+2*h))./(12*h);

x = .5;
f = @(x) 2*sin(5*x)+cos(3*x);
df = @(x) 10*cos(5*x)-3*sin(3*x);
h = 10.^(-(1:1:20)/2);

err = @(P, h) abs(df(x) - P(f, h, x));
fig = figure(1);
loglog(h, err(derivForwardDiff,h), 'b')
hold on
loglog(h, err(deriv3PointMid, h), 'color', '#136207')
loglog(h, err(deriv5PointMid, h), 'r')
loglog(h, h.^1, 'b--')
loglog(h, h.^2, '--', 'color', '#136207')
hold off
yticks(10.^(-20:2:0))
ylim([1e-20 1e1])
xlim([h(20) h(1)])
set(gca,'TickLabelInterpreter', 'latex','XMinorTick','off','YMinorTick','off')
xlabel('$h$','interpreter','latex','FontSize',15)
ylabel('$\left|D_{h} f(x)-f^{\prime}(x)\right|$','interpreter','latex','FontSize',15)
legend('Forward difference (Actual)','3 point mid-point (Actual)', ...
        '5 point mid-point (Actual)', 'Forward Diff. (Expected)', '3 point mid-point (Expected)', ...
        'interpreter','latex','FontSize',10, 'Location', 'northwest')
title('Expected error and Actual error rate in the approximation of $10 \cos (5 x)-3\sin (3 x)$ at $x=.5$', ...
    'interpreter','latex','FontSize',12)

saveas(fig, 'FDiff_35PtMidpointApprox.jpg')