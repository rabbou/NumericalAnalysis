%%
clear all;
close all;
clc;

%% Gauss-Seidel Iteration
% Movement of an object right or left with probability alpha = 1/2
% P_i = 1/2* P_{i?1} + 1/2 P_{i+1} for i = 1,2,...,n?1
make_A = @(n, alpha)  eye(n) - diag(alpha*ones(n-1, 1), -1) - diag((1-alpha)*ones(n-1, 1), 1);

alpha = 1/2;
% n = 10
n = 10;
b = zeros(n-1, 1);
b(1) = alpha;
x = Gauss_Seidel(make_A(n-1, alpha), b, zeros(n-1, 1), 10000, 1e-4);
% disp(x.');
% disp((make_A(n-1, alpha)*x).');

% n = 50
n = 50;
b = zeros(n-1, 1);
b(1) = alpha;
x = Gauss_Seidel(make_A(n-1, alpha), b, zeros(n-1, 1), 10000, 1e-4);
% disp(x.');
% disp((make_A(n-1, alpha)*x).');

% n = 100
n = 100;
b = zeros(n-1, 1);
b(1) = alpha;
x = Gauss_Seidel(make_A(n-1, alpha), b, zeros(n-1, 1), 10000, 1e-4);
% disp(x.');
% disp((make_A(n-1, alpha)*x).');

%% alpha = 1/3
alpha = 1/3;
% 10
n = 10;
b = zeros(n-1, 1);
b(1) = alpha;
x = Gauss_Seidel(make_A(n-1, alpha), b, zeros(n-1, 1), 10000, 1e-4);
% disp(x.');
% disp((make_A(n-1, alpha)*x).');

% 50
n = 50;
b = zeros(n-1, 1);
b(1) = alpha;
x = Gauss_Seidel(make_A(n-1, alpha), b, zeros(n-1, 1), 10000, 1e-4);
% disp(x.');
% disp((make_A(n-1, alpha)*x).');

% 100
n = 100;
b = zeros(n-1, 1);
b(1) = alpha;
x = Gauss_Seidel(make_A(n-1, alpha), b, zeros(n-1, 1), 10000, 1e-4);
% disp(x.');
% disp((make_A(n-1, alpha)*x).');

%%
function x = Gauss_Seidel(A, b, XO, N, TOL)
k = 1;
n = size(b, 1);
x = zeros(n, 1);
while k <= N
    for i = 1:n
        x(i) = (b(i) - A(i, 1:i-1)*x(1:i-1) - A(i, i+1:n)*XO(i+1:n))/A(i, i);
    end
    if norm(x - XO) < TOL
        fprintf('Success after %d iterations\n', k);
        return
    end
    k = k + 1;
    XO = x;
end
fprintf('Failed for Max iterations: %d\n', N);
return
end