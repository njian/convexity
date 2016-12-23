%----------------********************************--------------------------
% pnConvexCondMC.m
% Description:
%   Uses conditional Monte Carlo method to estimate the probability of
%   convex of the given posterior information.
%
% Inputs:
%   Nrand: the standard normal or T samples used to generate a random
%   direction on the unit shell.
%   xv: d-by-r, each column is a d-dimensional design point to sample on.
%   mu: posterior mean.
%   Lambdahalf: square root of the posterior sampling variance.
%   K: number of samples used in estimating pnConvex by Monte Carlo, i.e.
%      the dimension of the output.
%   nu: posterior hyperparameter (only nonempty if Gammaknown is 1).
%   GammaKnown: whether assuming the sampling variance is known.
%
% Outputs:
%   p: K-by-1 vector, each component is an estimate to pConvex.
%   I: K-by-1 vector, each component is a "free" vanilla Monte Carlo 
%      indicator of convexity.
%
% Author: Nanjing Jian, Cornell University. 
% Reference:
%   Jian and Henderson, 2016, "Estimating the Probability that a Function
%   Observed with Noise is Convex", submitted to INFORMS Journal of
%   Computing.
%----------------********************************--------------------------
function [ p, I ] = pnConvexCondMC( Nrand, x, mu, Lambdahalf, nu, GammaKnown )
[~, r] = size(x);
[K, ~] = size(Nrand);
p = zeros(K,1);
I = zeros(K,1);

for k = 1:K
    %% Conditional MC estimator (p)
    % Generate a random direction on the unit shell
    U = (Nrand(k,:) ./ norm(Nrand(k,:), 2))';
    D = Lambdahalf*U; % a random direction on the unit ellipsoid defined by Lambdahalf
    
    % Find the minimum and maximum stepsize along direction D so that the
    % vector is still convex.
    tmin = findStepSize(x, mu', D, 'min');
    tmax = findStepSize(x, mu', D, 'max');
    
    if tmin == Inf || tmax == -Inf || isnan(tmin) || isnan(tmax) || tmin > tmax % infeasible
        p(k) = 0;
    elseif tmin == -Inf && tmax == Inf % unbounded
        p(k) = 1; 
    else
        if GammaKnown == 0
            c = nu - r + 1;
            p(k) = sign(tmax) * 0.5 * fcdf(tmax^2 / r, r, c) - ...
                   sign(tmin) * 0.5 * fcdf(tmin^2 / r, r, c);
        elseif GammaKnown == 1
            p(k) = sign(tmax) * 0.5 * chi2cdf(tmax^2, r) - ...
                   sign(tmin) * 0.5 * chi2cdf(tmin^2, r);
        end
    end
    
    %% Free vanilla MC estimator (I)
    % Take a sample of the step size T. If it is between a and b, it means
    % the current sample is convex.
    temp = rand();
    if GammaKnown == 0
        if temp >= 1/2
            % take a positive sample
            sampleT = sqrt(finv(2 * (temp - 1/2), r, nu));
        else
            % take a negative smaple
            sampleT = -sqrt(finv(2 * (-1) * (temp - 1/2), r, nu));
        end
        
        if sampleT >= tmin && sampleT <= tmax
            I(k) = 1;
        end
    
    elseif GammaKnown == 1
        if rand() >= 1/2
            % take a positive sample
            sampleT = sqrt(chi2inv(2 * (rand() - 1/2), r));
        else
            % take a negative smaple
            sampleT = -sqrt(chi2inv(2 * (-1) * (rand() - 1/2), r));
        end  
        
        if sampleT >= tmin && sampleT <= tmax
            I(k) = 1;
        end
    end
     
end

end

