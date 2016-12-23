%----------------********************************--------------------------
% LRstudentt.m
% Description:
%   Calculate the likelihood ratio of two T densities evaluated at a
%   given sample.
% Inputs:
%   sample: the value that the two T densities are evaluated at.
%   mu0, Sigma0, kappa0, nu0: parameters of the T density on the
%   denominator.
%   mu, Sigma, kappa, nu: parameters of the T density on the numerator.
% Outputs:
%   likelihood: the likelihood ratio of the T densities.
%----------------********************************--------------------------
function [ likelihood ] = ...
         LRstudentt(sample, mu0, Sigma0, kappa0, nu0, mu, Sigma, kappa, nu)
     
[K, r] = size(sample);
likelihood = zeros(K,1);

dof0 = nu0-r+1;
dof = nu-r+1;

scale0 = Sigma0/(kappa0*dof0);
scale = Sigma/(kappa*dof);
for k = 1:K
    Y = sample(k,:);
    pdf = mvtpdf(Y - mu, scale, dof);
    pdf0 = mvtpdf(Y - mu0, scale0, dof0);
    likelihood(k) = pdf / pdf0;
end

end

