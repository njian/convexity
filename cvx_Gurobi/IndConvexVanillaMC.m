%----------------********************************--------------------------
% IndConvexVanillaMC.m
% Description:
%   Gives a vector of the indicators of whether each row of sample is
%   convex. Used for estimating the probability that the samples on xv is
%   convex using vanilla Monte Carlo method. This is a version that 
%   requires the installation of cvx and Gurobi.
%
% Inputs:
%   sample: the samples on xv.
%   xv: d-by-r, each column is a d-dimensional design point to sample on.
%
% Outputs:
%   Ind: K-by-1 indicators of whether the input is convex, where K is the
%   number of Monte Carlo samples.
%
% Author: Nanjing Jian, Cornell University. 
% Reference:
%   Jian and Henderson, 2016, "Estimating the Probability that a Function
%   Observed with Noise is Convex", submitted to INFORMS Journal of
%   Computing.
%----------------********************************--------------------------
function [ Ind ] = IndConvexVanillaMC(sample, xv)
% Algorithm 2: For when Gamma is unknown (Normal-inv-Wishart case)
[~, r] = size(xv);
[K, ~] = size(sample);

Ind = ones(K,1);
for k = 1:K
    y = sample(k,:)';
    
    for i = 1:r
        % Solve for LS(i)
        cvx_begin quiet
            variable a(d)
            variable b
            expression mul(r,1)
            subject to
              mul = ones(r,1); mul(i) = 0;
              y(i) == (a'*xv(:,i))' + b;
              y.*mul >= (a'*xv)'.*mul + b*mul;
        cvx_end
        
        if strcmp(cvx_status,'Infeasible') % if there is NO feasible soln for a and b
            Ind(k) = 0;
            break; % infeasible LP(i) found, proceed to next k
        end
        
    end

end
end