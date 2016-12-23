%----------------********************************--------------------------
% IndConvexVanillaMC.m
% Description:
%   Gives a vector of the indicators of whether each row of sample is
%   convex. Used for estimating the probability that the samples on xv is
%   convex using vanilla Monte Carlo method.
% Inputs:
%   sample: the samples on xv.
%   xv: d-by-r, each column is a d-dimensional design point to sample on.
% Outputs:
%   Ind: K-by-1 indicators of whether the input is convex, where K is the
%   number of Monte Carlo samples.
% Author: Nanjing Jian, Cornell University. 
% Reference:
%   Jian and Henderson, 2016, "Estimating the Probability that a Function
%   Observed with Noise is Convex", submitted to INFORMS Journal of
%   Computing.
%----------------********************************--------------------------
function [ Ind ] = IndConvexVanillaMC(sample, xv)
% Algorithm 2: For when Gamma is unknown (Normal-inv-Wishart case)
[d, r] = size(xv);
[K, ~] = size(sample);

options = optimset('Display', 'off'); % turn off screen output for linprog
Ind = ones(K,1);
for k = 1:K
    y = sample(k,:)';
    
        for i = 1:r
            % Solve for LS(i)
            
            % set up equality constraints x_i*a_i + b_i = g_i
            Aeq = [xv(:,i)', 1];
            beq = y(i);
    
            % set up inequality constraints : x_j*a_i + b_i <= g_j
            % first ignore j ~= i
            A = [xv', ones(r,1)];
            b = y;
            % then delete ith row
            A(i,:) = [];
            b(i) = [];
    
            % check feasibility
            c = zeros(d+1,1);
            [~, ~, exitflag] = linprog(c,A,b,Aeq,beq,[],[],[],options);
            if exitflag ~= 1 % if there is NO feasible soln for a and b
                Ind(k) = 0;
                break; % infeasible LP(i) found, proceed to next k
            end
        
        end 

end
end

