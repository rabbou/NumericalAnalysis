%%
clear all;
close all;
clc;

%% Runge-Kutta
% We have 10 particles living in the Reals with initial conditions x_j(0) = j 
% and x'j(0) = 0 except x_1(0) = .5
% u_j = x_j ? j is the displacement.
% The force equations are u_1''(t)= ?10(u_(t)?u_2(t))
% u_n''(t)=?10(u_n(t)?u_{n-1}(t))
% u_j''(t) = ?10 (2u_j(t) ? u_{j-1}(t) ? u_{j+1}(t))
% the initial value problem which corresponds is y(0)=0, y' = y?t
% Approximate the solutions for this IVP using 4 th order Runge- Kutta
N = 10;
one = ones(N, 1) ;
A = diag(2 * one, 0) - diag(one(1:N-1), -1) - diag(one(1:N-1), 1);
A(1,1) = 1; 
A(N,N) = 1;
f = [zeros(N, N) eye(N); -10*A zeros(N, N)];
y0 = zeros(2*N,1);
y0(1) = -.5;
h = .05;
W = runge_kutta_4(f, y0, N, h, 8);
X = W(1:10,:) + (1:10)';

plot1 = scatter(X(:,1), zeros(N,1), 30, 'fill');
xlim([0,11]);
ylim([-1,1]);
xlabel('x')
ylabel('y')
title('Runge-Kutta', 'interpreter', 'latex', 'FontSize',15)
for k=2:8/h
    plot1.XData = X(:,k);
    pause(h)
end


%% approximate y(8) for the following time-steps: h = 2?s for s = 4,5,...,13
Y = zeros(10,10);
for j=1:10
    h = 2^(-1*(j+3));
    W = runge_kutta_4(f, y0, N, h, 8);
    Y(:,j) = W(1:10,length(W)) + (1:10)';
end

H = zeros(9);
Error = zeros(9);
for j=1:9
    Error(j) = norm((Y(:,j) - Y(:,10)), 1);
    H(j) = 2^(-1*(j+3));
end

loglog(H, Error, H, H.^4)
title('Error vs. step size', 'interpreter', 'latex', 'FontSize',15)
xlabel('$\log h$', 'interpreter', 'latex', 'FontSize',13)
ylabel('Error', 'interpreter', 'latex', 'FontSize',13)


%% functions

function r = runge_kutta_4(f, y0, N, h, T)
p = T / h;
r = zeros(2*N,p);
r(:,1) = y0;

for i=2:p
    w = r(:,i-1);
    k1 = h * f * w;
    k2 = h * f * (w + k1/2);
    k3 = h * f * (w + k2/2);
    k4 = h * f * (w + k3);
    r(:,i) = w + (k1 + 2*k2 + 2*k3 + k4)/6;
end
end