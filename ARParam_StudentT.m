prin%----------------********************************--------------------------
% ARParam_StudentT.m
% Description:
%   Calculate the maximum likelihood ratio of two T densities.
% Inputs:
%   l: number of iterations since a brand new set of sample was obtained.
%   sample: K-by-r, the value that the likelihood ratio is calculated for.
%   mu0, Sigma0, kappa0, nu0: parameters of the T density on the
%   denominator.
%   mu, Sigma, kappa, nu: parameters of the T density on the numerator.
% Outputs:
%   ratio: a K-by-1 vector, each component is the likelihood ratio of the
%   given T densities evaluated at sample.
%   C: the maximum likelihood ratio of the given T densities.
%   heavytail: if the heavytail of the likelihood ratio makes it impossible
%   to find C.
%----------------********************************--------------------------
function [ ratio, C, heavytail ] = ARParam_StudentT( l, sample, mu0, kappa0, nu0, Sigma0, mu, kappa, nu, Sigma )
[K, ~] = size(sample);
ratio = zeros(1,K);

% Numerically maximize LRstudentt (or minimize its negative)
LRratio = @(y) - LRstudentt(y, mu0, Sigma0, kappa0, nu0, mu, Sigma, kappa, nu);
options = optimset('MaxFunEvals', 1e10, 'MaxIter', 1e10, 'Display', 'notify');
[~, C, exitflag, output] = fminsearch(LRratio, 1, options);
C = -C;

heavytail = 0;
if exitflag == 0 % Maximum number of iterations has been exceeded
    heavytail = 1;
end

for k = 1:K
    Y = sample(k,:);
    ratio(k) = LRstudentt(Y, mu0, Sigma0, kappa0, nu0, mu, Sigma, kappa, nu);
end

end

