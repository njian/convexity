%----------------********************************--------------------------
% findStepSize.m
% Description:
%   Find the minimum or maximum stepsize along a given vector so that the 
%   resulted vector is still convex.
% Inputs:
%   xv: d-by-r, the design points to sample on.
%   mu_n: 1-by-r, the center of the posterior distribution.
%   D: the direction (Lambdahalf*Z) that is uniformly distributed on the
%   surface of the ellipsoid defined by the posterior distribution.
%   solveFor: 'min' or 'max' stepsize to solve for.
% Outputs:
%   t: the min or max stepsize.
%----------------********************************--------------------------
function [ t ] = findStepSize( xv, mu_n, D, solveFor)
[d, r] = size(xv);

% Decision variables: a1, a2, ..., ar in R^(d by 1), b1, b2, ..., br in
% R, and t in R.
if strcmp(solveFor, 'min')
    % minimize t
    f = [zeros(d*r+r, 1); 1]; 
elseif strcmp(solveFor, 'max')
    % minimize -t
    f = [zeros(d*r+r, 1); -1];
end

% Construct the equality constraints
Aeq = [];
% x1', x2', ..., xr' on the diagonal, and 0 off-diagonal
for i = 1:r
    Aeq = blkdiag(Aeq, xv(:,i)');
end
Aeq = [Aeq, eye(r), -D];
beq = mu_n';

% Construct the inequality constraints
% First do all th r^2 constraints ignoring i \neq j
A = [];
b = [];
for i = 1:r
    A = [A; zeros(r,d*(i-1)), xv', zeros(r,d*(r-i)), zeros(r,i-1), ones(r,1), zeros(r,r-i), -D];
    b = [b; mu_n];
end
% Then delete the i=j rows
for i = 1:r
    A = removerows(A, 'ind', r*(i-1)+1);
    b = removerows(b, 'ind', r*(i-1)+1);
end

% Solve for t!
options = optimset('Display', 'off', 'Algorithm', 'simplex', 'MaxIter', 1000); % turn off screen output for linprog
[~, t, exitflag, output] = linprog(f, A, b, Aeq, beq, [], [], [], options);

if strcmp(solveFor, 'max')
    if exitflag == -2 || exitflag == -5 % primal infeasible, or both primal and dual infeasible
        t = -Inf;
    elseif exitflag == -3 % unbounded
        t = +Inf;
    else
        % minimize -t
        t = -t;
    end
    
elseif strcmp(solveFor, 'min')
    if exitflag == -2 || exitflag == -5 % primal infeasible, or both primal and dual infeasible
        t = Inf;
    elseif exitflag == -3
        t = -Inf;
    end
end

% fprintf('linprog, t = %d, exitflag = %d, output message = %s. \n', t, exitflag, output.message);


end

