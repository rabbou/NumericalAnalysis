%%
clear all;
close all;
clf;

%% Find fixed point of g(x) = pi + 0.5 sin(x/2) on [0, 2pi]
fixed_p_it(pi, 1e-2, 20, @(x) pi + 0.5*sin(x/2))

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

%% Finding 3^(1/3) using Fixed point iterations, Newton's method and a Cubic Method

p0 = 1.3;
N0 = 20;
g_fp = @(x) x-.1*(x^3-3);
f = @(x) x^3-3;
df = @(x) 3*x^2;
ddf = @(x) 6*x;
g_nm = @(x) x-f(x)/df(x);
g_cm = @(p0) p0 - f(p0)/df(p0) - ddf(p0)/(2*df(p0))*(f(p0)/df(p0))^2;

true_p = cubic(p0, N0, f, df, ddf);

fixed_point = e(p0, N0, g_fp, true_p);
newt_met = e(p0, N0, g_nm, true_p);
cub_met = e(p0, N0, g_cm, true_p);

fig = figure;
a1 = semilogy(1:N0,fixed_point(1:N0), 'linewidth', 1.5);
t1 = "Fixed Point Iteration";
hold on
a2 = semilogy(1:N0,newt_met(1:N0), 'linewidth', 1.5);
t2 = "Newton's Method";
hold on
a3 = semilogy(1:N0,cub_met(1:N0), 'linewidth', 1.5);
t3 = "Cubic Method";
legend([a1;a2;a3], t1,t2,t3,'FontSize',12,'interpreter','latex');

xlabel('Iteration number','interpreter','latex','FontSize',15);
ylabel('Absolute error','interpreter','latex','FontSize',15);
title('Error Convergence in approximating $\sqrt[3]{3}$','interpreter','latex','FontSize',15);

saveas(fig, 'CubeRootThree.jpg');

%% Find three positive zeros for (cos(x)+sin(sqrt(2)*x)) * exp(-x)

f1 = @(x) (cos(x)+sin(sqrt(2)*x)) .* exp(-x);
df1 = @(x) -exp(-x)*(sin(x)+sin(sqrt(2)*x)+cos(x)-sqrt(2)*cos(sqrt(2)*x));

x = linspace(1,3*pi);
y = f1(x);
figure
line([0 0], ylim);
line(xlim, [0 0]);
plot(x, y, 'linewidth', 1.5);
hold on
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
hold off

% according to the plot, we can use the following starting points:
p01 = 2;
p02 = 4.5;
p03 = 7;
TOL = 1e-10;
newton(p01, TOL, 20, f1, df1)
newton(p02, TOL, 20, f1, df1)
newton(p03, TOL, 20, f1, df1)

%%
function error_terms = e(p0, N0, g, true_p)
error_terms = zeros(1,N0);
Iter = 0;
while Iter < N0
    Iter = Iter+1;
    xn = p0;    
    p0 = g(xn);    
    error_terms(Iter) = abs(xn-true_p);
end
end

function p = fixed_p_it(p0, TOL, N0, g)
i = 1;
while i <= N0
   p =  g(p0);
   if abs(p - p0) < TOL
       fprintf('Fixed point obtained after %d itetations:', i)
       return
   end
   i = i+1;
   p0 = p;
end
fprintf('The method failed after %d iterations with p = %f\n', N0, p);
end

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

function p = cubic(p0, N0, f, df, ddf)
i = 1;
while i <= N0
   p =  p0 - f(p0)/df(p0) - ddf(p0)/(2*df(p0))*(f(p0)/df(p0))^2;
   i = i+1;
   p0 = p;
end
end