%----------------********************************--------------------------
% LRstudentt.m
% Description:
%   Calculate the likelihood ratio of two Normal densities evaluated at
%   a given sample.
% Inputs:
%   sample: the value that the two Normal densities are evaluated at.
%   mu0, Lambda0, iLambda0: mean and variance of the Normal density on the
%   denominator (iLambda0 is the inverse of Lambda0).
%   mu, Lambda, iLambda: mean and variance of the Normal density on the
%   numerator (iLambda is the inverse of Lambda).
% Outputs:
%   likelihood: the likelihood ratio of the Normal densities.
%----------------********************************--------------------------
function [ likelihood ] = ...
         LRnormal( sample, mu0, Lambda0, iLambda0, mu, Lambda, iLambda )
     
[K, ~] = size(sample);
likelihood = zeros(K,1);
for k = 1:K
    % This is more robust than LR(k) = mvnpdf(Y,mu,Lambda)/mvnpdf(Y,mu0,Lambda0);
    Y = sample(k,:);
    diff1 = Y-mu;
    diff0 = Y-mu0;
    expo = -1/2*((diff1*iLambda)*(diff1')-(diff0*iLambda0)*(diff0'));
    rateDet = det(Lambda0)/det(Lambda);
    likelihood(k) = (rateDet)^(1/2) * exp(expo);
end

if min(likelihood) < 0
   fprintf('When min(LR) = %d, Y - mu is %d. \n', min(likelihood), diff1) 
end

end

