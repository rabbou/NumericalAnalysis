%%
clear all;
close all;
clc;

%% Newton's method
% Find solution for x^3 ? 2x^2 ? 5 = 0 on [1,4]
newton(2.5, 1e-4, 20, @(x) x^3-2*x^2-5, @(x) 3*x^2 - 4*x)
% Find solution for x^3 + 3x^2 ? 1 = 0 on [?3,?2]
newton(-2.5, 1e-4, 20, @(x) x^3+3*x^2-1, @(x) 3*x^2 + 6*x)
% Find solution for 1/2+x^2/4-x*sin(x)-1/2*cos(2*x) = 0 with p0 = pi/2,
% 5pi, 10pi
newton(pi/2, 1e-5, 50, @(x) 1/2+x^2/4-x*sin(x)-1/2*cos(2*x), @(x) x/2-sin(x)-x*cos(x)+sin(2*x))
newton(5*pi, 1e-5, 50, @(x) 1/2+x^2/4-x*sin(x)-1/2*cos(2*x), @(x) x/2-sin(x)-x*cos(x)+sin(2*x))
newton(10*pi, 1e-5, 100000, @(x) 1/2+x^2/4-x*sin(x)-1/2*cos(2*x), @(x) x/2-sin(x)-x*cos(x)+sin(2*x))

%% Find three positive zeros for (cos(x)+sin(sqrt(2)*x)) * exp(-x)

f1 = @(x) (cos(x)+sin(sqrt(2)*x)) .* exp(-x);
df1 = @(x) -exp(-x)*(sin(x)+sin(sqrt(2)*x)+cos(x)-sqrt(2)*cos(sqrt(2)*x));

x = linspace(0,3*pi);
y = f1(x);
fig = figure(1);
plot(x, y, 'b', 'linewidth', 1.5);
hold on
plot(x, zeros(1,size(x,2)),'--r','linewidth', 1);
hold off
xlim([0 3*pi]);
xlabel('$x$','interpreter','latex','FontSize',15);
ylabel('$f(x)$','interpreter','latex','FontSize',15);
title('Plot of $f(x) = \left(\cos(x)+\sin\left(\sqrt2x\right)\right) \cdot e^{-x}$','interpreter','latex','FontSize',15);
saveas(fig, 'TrigFcnZerosNewton.jpg');

% according to the plot, we can use the following starting points:
p01 = 2;
p02 = 4.5;
p03 = 7;
TOL = 1e-10;
newton(p01, TOL, 20, f1, df1)
newton(p02, TOL, 20, f1, df1)
newton(p03, TOL, 20, f1, df1)

%%
function p = newton(p0, TOL, N0, f, df)
i = 1;
while i <= N0
   p = p0 - f(p0)/df(p0);
   if abs(p - p0) < TOL
       fprintf('Convergence obtained using Newton after %d itetations', i)
       return
   end
   i = i+1;
   p0 = p;
end
fprintf('The method failed after %d iterations\n', N0);
end