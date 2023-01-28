% 
% 
% This is to solve the diffusion in 1D/2D In a 1D/2D propagation, the operator is split into
% the ODE solution, and the PDE solution.
% 
% 
function x = TDMAsolver(gam,c,d)
% a, b, c are the column vectors for the compressed tridiagonal matrix, d is the right vector
n = length(c); % n is the number of rows

a = -gam;
b = 1.0+2*gam;
an = -2*gam;

% Modify the first-row coefficients
c(1) = c(1) / b;    % b is always 1+2*gam and b is always on RHS - so that is also a constant
d(1) = d(1) / b;
 
for i = 2:n-1
    temp = b - a * c(i-1); % a = -gam, and it is on RHS, b is 1+2*gam and only on RHS. c occurs recursively!
    c(i) = c(i) / temp;
    d(i) = (d(i) - a * d(i-1))/temp; % a is only on RHS, here it is -gam
end;
 
d(n) = (d(n) - an * d(n-1))/(b - an * c(n-1)); % a is only on RHS, here it is -2*gam
 
% Now back substitute.
x(n) = d(n);
for i = n-1:-1:1
    x(i) = d(i) - c(i) * x(i + 1);
end;
end
