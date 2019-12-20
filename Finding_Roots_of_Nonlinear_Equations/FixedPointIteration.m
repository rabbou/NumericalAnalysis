%%
clear all;
close all;
clc;

%% Find fixed point of g(x) = pi + 0.5 sin(x/2) on [0, 2pi]
fixed_p_it(pi, 1e-2, 20, @(x) pi + 0.5*sin(x/2))

%%
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