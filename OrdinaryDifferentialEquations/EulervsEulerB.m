%%
clear all;
close all;
clc;

%% We have a Harmonic Oscillator q''(t) = -kq(t),q(0) = 0,q'(0) = 1
% let p(t) = x?(t). Let's evolve this system for t in [0, 20]:

k = 2;
h = 1e-2;
a = 0;
b = 20;
q0 = 0;
p0 = 1;

[p, q] = euler(a, b, h, p0, q0, k);
[p_B, q_B] = euler_B(a, b, h, p0, q0, k);
xx = a:h:b;
fig1 = figure(1);
plot(xx, q, xx, q_B)
legend('Euler', 'Euler-B','interpreter','latex','FontSize',12,'location','northwest');
title('Position vs. Time for a Harmonic Oscillator using Euler''s Methods',...
    'interpreter','latex','FontSize',15);
xlabel('Time $t$','interpreter','latex','FontSize',15);
ylabel('Position $q(t)$','interpreter','latex','FontSize',15);
saveas(fig1, 'EulervsEulerB.jpg')

%% Energy Loss
E = @(p, q) .5 * p.^2 + .5 * k * q.^2;

fig2 = figure(2);
plot(xx, E(p, q), xx, E(p_B, q_B))
legend('Euler', 'Euler-B','interpreter','latex','FontSize',15,'location','northwest');
title('Energy vs. Time for a Harmonic Oscillator using Euler''s Methods',...
    'interpreter','latex','FontSize',15);
xlabel('Time $t$','interpreter','latex','FontSize',15);
ylabel('Energy $E(t)$','interpreter','latex','FontSize',15);
saveas(fig2, 'EulervsEulerBEnergy.jpg')

% We can see clearly that energy is conserved much better in the Euler B method,
% while Euler?s method is dependent of time.

%%
function [p, q] = euler(a, b, h, p0, q0, k)
N = (b-a)/h;
p = zeros(1, N);
q = zeros(1, N);
p(:,1) = p0;
q(:,1) = q0;
for n = 1:N
    q(n+1) = q(n) + h * p(n);
    p(n+1) = p(n) - h * k*q(n);
end
end

function [p, q] = euler_B(a, b, h, p0, q0, k)
N = (b-a)/h;
p = zeros(1, N);
q = zeros(1, N);
p(1) = p0;
q(1) = q0;
for n = 1:N
    p(n+1) = p(n) - h * k * q(n);
    q(n+1) = q(n) + h * p(n+1);
end
end
