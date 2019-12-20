%%
clear all;
close all;
clc;

%% Bisection Method estimations
% Solution for 3x-e^x=0 on [1,2]
bisection(1, 2, 1e-5, 20, @(x) 3*x-exp(x))
% Solution for x^2-4x+4-log(x)=0 on [1,2]
bisection(1, 2, 1e-5, 20, @(x) x^2-4*x+4-log(x))
% Solution for x^2-4x+4-log(x)=0 on [2,4]
bisection(2, 4, 1e-5, 20, @(x) x^2-4*x+4-log(x))
% Solution for x^3-25=0 on [2,3]
bisection(2, 3, 1e-4, 20, @(x) x^3-25)

% approximation with accuracy 10^-4 to the solution of x^3 ? x ? 1 = 0
f15 = @(x) x^3-x-1;
a = 1;
b = 2;
FA = f15(a);
for i = 1:14
    p = a + (b - a)/2;
    FP = f15(p);
    if FA*FP > 0
        a = p;
        FA = FP;
    else
        b = p;
    end
end
fprintf('Result for 14 iterations: %f \n', p)

%% Evaluate P(x_0) for P(x) = a_n*x^n + a_{n-1}*x_{n-1} +···+ a_1*x + a_0
% using nested multiplications, assuming a_i = 1 for all i.

n = 10;
x_0 = .5;
% raise everything to n+1 since matlab does not include 0 indexes
a = ones([n+1,1]);
Px_0 = a(n+1);
for i = n:-1:1
    Px_0 = a(i) + x_0 * Px_0;
end
fprintf('Value for P(x): %f \n', Px_0)

function p = bisection(a, b, TOL, N0, f)
i = 1;
FA = f(a);
while i <= N0
    p = a + (b - a)/2;
    FP = f(p);
    if FP == 0 || (b - a)/2 <  TOL
        fprintf("number of iterations: %d", i)
        return
    end
    i = i+1;
    if FA*FP > 0
        a = p;
        FA = FP;
    else
        b = p;
    end
end
fprintf('Method failed after N0 iterations = %d\n', N0);
end
