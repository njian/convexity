%----------------********************************--------------------------
% findStepSize.m
% Description:
%   Find the minimum or maximum stepsize along a given vector so that the 
%   resulted vector is still convex. This is a version that requires the
%   installation of cvx and Gurobi.
%
% Inputs:
%   xv: the design points to sample on.
%   mu_n: the center of the posterior distribution.
%   D: the direction (Lambdahalf*Z) that is uniformly distributed on the
%   surface of the ellipsoid defined by the posterior distribution.
%   solveFor: 'min' or 'max' stepsize to solve for.
%
% Outputs:
%   t: the min or max stepsize.
%
% Author: Nanjing Jian, Cornell University. 
% Reference:
%   Jian and Henderson, 2016, "Estimating the Probability that a Function
%   Observed with Noise is Convex", submitted to INFORMS Journal of
%   Computing.
%----------------********************************--------------------------
function [ t ] = findStepSize_cvx( xv, mu_n, D, solveFor)
[~, r] = size(xv);

if strcmp(solveFor, 'min')
    % Full LP
    cvx_begin quiet
        variable t
        variable a(d,r)
        variable b(r)
        expression mul(r,1)
        minimize( t );
        subject to
           for i = 1:r
              mul = ones(r,1); mul(i) = 0;
              mu_n(i) + D(i)*t == (a(:,i)'*xv(:,i))' + b(i);
              mu_n.*mul + D*t.*mul >= (a(:,i)'*xv)'.*mul + b(i)*mul;
           end
    cvx_end

elseif strcmp(solveFor, 'max')
    % Full LP
    cvx_begin quiet
        variable t
        variable a(d,r)
        variable b(r)
        expression mul(r,1)
        maximize( t );
        subject to
           for i = 1:r
              mul = ones(r,1); mul(i) = 0;
              mu_n(i) + D(i)*t == (a(:,i)'*xv(:,i))' + b(i);
              mu_n.*mul + D*t.*mul >= (a(:,i)'*xv)'.*mul + b(i)*mul;
           end
    cvx_end
    
end

if cvx_optval == -Inf %infeasible
    t = -Inf;
elseif cvx_optval == Inf %unbounded
    t = Inf;
end

% fprintf('gurobi, t = %d, cvx_status = %s. \n', t, cvx_status);

end

