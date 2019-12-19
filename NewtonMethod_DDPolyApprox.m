%%
clear all;
close all;
clf;

%% Interpolate 3^x on [0,2.5] using Newton approximations
% Evaluate at x = 0.5 ==> x = 3^(0.5) ==> x^2 = 3
newton(2, 1e-8, 50, @(x) x^2 - 3, @(x) 2*x)

% Evaluate at x = 1 ==> x = 3
newton(2, 1e-8, 50, @(x) x-3, @(x) 1)

% Evaluate at x = 1.5 ==> x = 3^(1.5) ==> x^2 = 3^3 = 27
newton(2.5, 1e-8, 50, @(x) x^2-27, @(x) 2*x)

%% Divided Differences Polynomial approximation of 3^x using the three query points
xi = [.5, 1, 1.5];
f = @(x) 3.^x;
fi = f(xi);

x = 0:.01:2.5;
fx = f(x);
n = max(size(xi));
DD = DividedDifferences(xi, fi);
p_x = DD_Eval(DD, x, xi);

fig = figure('Position',[0,0,1000,400]);
subplot(1,2,1)
plot(x,p_x,'r', x,fx, 'b', 'Linewidth', 1.1);
legend('Polynomial Approx','$f(x)$','interpreter','latex','FontSize',13,'Location','northwest')
title('Polynomial Approximation of $3^x$ using Divided Differences','interpreter','latex','FontSize',13)
xlabel('x','interpreter','latex','FontSize',13)
ylabel('y','interpreter','latex','FontSize',13)

subplot(1,2,2)
E = abs(p_x-fx)./abs(fx);
plot(x,E, 'color', '#77AC30', 'Linewidth', 1.1);
title('Point-wise relative error of the polynomial approximation of $3^x$','interpreter','latex','FontSize',13)
xlabel('x','interpreter','latex','FontSize',13)
ylabel('Relative Error','interpreter','latex','FontSize',13)
saveas(fig, '3PowerxDDPolyApprox.jpg');


%%
function f_x = DD_Eval(div_diff, x, xi)
n = max(size(div_diff));
I = ones(size(x));
f_x = div_diff(n)*I;
for i = fliplr(1:n-1)
    f_x = f_x .* (x-xi(i)) + div_diff(i)*I;
end
end

function [div_diff,full_div] = DividedDifferences(xi,f_xi)
n = max(size(xi));
full_div = zeros(n,n);
full_div(:,1) = f_xi(:);
for i = 2:n
   full_div(i:n,i) = (full_div( i:n,i-1) - full_div((i-1):(n-1),i-1)) ./ (xi(i:n)'-xi(1:(n-i+1))');
end
div_diff = diag(full_div);
end

function p = newton(p0, TOL, N0, f, ff)
i = 1;
while i <= N0
   p = p0 - f(p0)/ff(p0);
   if abs(p - p0) < TOL
       fprintf('Convergence obtained using Newton after %d itetations \n', i)
       return
   end
   i = i+1;
   p0 = p;
end
fprintf('The method failed after %d iterations\n', N0);
end
