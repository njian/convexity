%----------------********************************--------------------------
% ARParam_Normal.m
% Description:
%   Calculate the maximum likelihood ratio of two Normal densities.
% Inputs:
%   l: number of iterations since a brand new set of sample was obtained.
%   sample: K-by-r, the value that the likelihood ratio is calculated for.
%   mu0, Lambda0, iLambda0: mean and variance of the Normal density on the
%   denominator (iLambda0 is the inverse of Lambda0).
%   mu, Lambda, iLambda: mean and variance of the Normal density on the
%   numerator (iLambda is the inverse of Lambda).
% Outputs:
%   ratio: a K-by-1 vector, each component is the likelihood ratio of the
%   given T densities evaluated at sample.
%   C: the maximum likelihood ratio of the given T densities.
%----------------********************************--------------------------
function [ ratio, C ] = ARParam_Normal( l, sample, mu0, iLambda0, mu, iLambda, Lambda0, Lambda, Gamma )
% Calculate the likelihood ratio of normal(mu, Sigma)/normal(mu0, Sigma0)
[K, ~] = size(sample);
ratio = zeros(1,K);

temp0 = iLambda0*mu0';
temp = iLambda*mu';

% Also tried ratio of mvnpdf at Ymin = 1/l*Gamma*(temp-temp0), same answer
% but more numerical error later on (divide by 0)
C = sqrt((det(Lambda0)/det(Lambda))*exp(1/l*(temp-temp0)'*Gamma*(temp-temp0) + mu0*temp0 - mu*temp));

for k = 1:K
    Y = sample(k,:);
    ratio(k) = LRnormal( Y, mu0, Lambda0, iLambda0, mu, Lambda, iLambda );
end
end

