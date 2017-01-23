%----------------********************************--------------------------
% algorithm.m
% Description: 
%   Performs sequential update to the posterior and calls Monte Carlo
%   methods with different variance reduction tricks to estimate the
%   posterior probability that the sample function is convex.
%
% Calls:
%   sampleFcn.m, MonteCarlo.m, (for all)
%   IndConvexVanillaMC(_cvx).m, (for vanilla Monte Carlo)
%   ARParam_StudentT.m, ARParam_Normal.m, LRstudentt.m, LRnormal.m, (for
%   change of measure and generalized acceptance/rejection methods)
%   pnConvexCondMC.m, findStepSize(_cvx).m, (for conditional Monte Carlo)
%
% Inputs:
%   x: the design points to sample on.
%   mu0, Lambda0: prior (Normal) mean and covariance.
%   kappa0, nu0, Sigma0: prior (inv-Wishart) parameters.
%   Gamma: sampling covariance matrix.
%   maxIte: max number of iterations allowed for Bayesian updates.
%   K: number of samples used in MC estimation of p(convex).
%   GammaKnown: 1 if assume Gamma is known, 0 otherwise.
%   method: 0 for 'vanillaMC', 1 for 'conditionalMC', 2 for 'changeMeasure',
%   or 3 for 'acceptanceRejection'. See reference for details.
%
% Outputs:
%   phat: estimated posterior probability that the sampleFcn is convex.
%   half: half width of phat with 95% confidence.
%   mu: posterior mean of the function values when algorithm stops.
%   Lambda: posterior covariance of the function values when algorithm 
%           stops.
%   time: time spent in estimating phat via Monte Carlo in each iteration.
%   efficiency: 1/(variance*time) of the Monte Carlo estimator for phat.
%
% Author: Nanjing Jian, Cornell University.
% Reference:
%   Jian and Henderson, 2016, "Estimating the Probability that a Function
%   Observed with Noise is Convex", submitted to INFORMS Journal of
%   Computing.
%----------------********************************--------------------------
function [ phat, half, mu, Lambda, time, efficiency ] = ...
            algorithm( x, mu0, Lambda0, kappa0, nu0, Sigma0, Gamma, ...
                        maxIte, K, GammaKnown, method )
% Initialize
[~, r] = size(x);
phat = zeros(maxIte,1);
half = zeros(maxIte,1);
stdev = zeros(maxIte,1);
time = zeros(maxIte,1);
efficiency = zeros(maxIte,1);
fprintf('==========================SETUP================================\n')
prompt = 'Please enter the number of samples in each iteration of posterior update:\n (integer >= 1, default = 1)  \n';
s = input(prompt); % #samples in each iteration
checkPrompt(s, '>=1Integer')

% Prompt for additional inputs if "skipping" every few iterations.
initialItn = 0;
probUpdate = 1;
if method ~= 0
    prompt = 'Please enter the initial number of iterations WITHOUT reusing samples:\n (integer >= 1, default = 10) \n';
    initialItn = input(prompt);
    checkPrompt(initialItn, 'positiveInteger');
    prompt = 'And the samples from an earlier iteration will be reused for this many iterations:\n (integer >=1, 1 means to generate new samples every iteration) \n';
    probUpdate = input(prompt);
    checkPrompt(probUpdate, 'positiveInteger');
end
fprintf('========================SETUP INFO=============================\n')
fprintf('maxIte = %d, K = %d, method = %d. \n', maxIte, K, method);
fprintf('initialItn = %d, probUpdate = %d. \n', initialItn, probUpdate);
fprintf('=======================START ITERATIONS========================\n')

%% When the sampling covariance Gamma is known: Normal-Normal case.
if GammaKnown == 1
    iGamma = inv(Gamma);
    last_mu = zeros(r,probUpdate)';
    last_Lambda = zeros(r,r*probUpdate)';
    last_iLambda = zeros(r,r*probUpdate)';
    
    n=1;
    while (n <= maxIte)
        fprintf('Iteration %d: \n', n);
        
        Yn = sampleFcn(x, s);  % new samples
        if s ~= 1
            Yn = mean(Yn);
        end
        
        % WITH ORIGINAL FORMULA
        if n == 1
            iLambda0 = inv(Lambda0);
        end
        iLambda = iLambda0 + s * iGamma;
        Lambda = inv(iLambda);
        mu = (iLambda \ (iLambda0 * mu0' + s * iGamma * Yn'))';
        
        % WITH FORMULA BY SMW, useful only when Gamma is diagonal
        % temp1 = Lambda0*s/Gamma;
        % temp = temp1/(eye(r) + temp1);
        % Lambda = Lambda0 - temp*Lambda0;
        % mu = (mu0' + temp*(Yn-mu0)')';
        
        % In case Lambda becomes singular due to numerical error.
        minEig = min(eig(Lambda));
        fprintf('Min eigenvalue of Lambda before adjuested: %d \n', minEig);
        while minEig < 0
            Lambda = perturbPSD(Lambda, minEig);
            minEig = min(eig(Lambda));
            fprintf('Min eigenvalue of Lambda after adjuested: %d \n', minEig);
        end
        Lambdahalf = sqrtm(Lambda);
        Y = Nrand*Lambdahalf + ones(K,1)*mu;
        
        % DIFFERENT METHODS
        switch method
            % Vanilla MC --------------------------------------------------
            case 0
                tic;
                Ind = IndConvexVanillaMC(Y, x);
                [ phat(n), stdev(n), half(n) ] = MonteCarlo(Ind);
                time(n) = toc;
                efficiency(n) = 1/(time(n) * stdev(n)^2);
                fprintf('CI for Pconvex: %.2e +/- %.2e. \n', phat(n), half(n));
                fprintf('time of the iteration: %.2f. \n', time(n));
                fprintf('Efficiency (1/sigma^2t) of the iteration: %.2f. \n', efficiency(n));
      
            % Conditional MC ----------------------------------------------
            case 1
                tic;
                % the first initialItn iterations
                if (n <= initialItn) || (mod((n - initialItn), probUpdate) == 0)
                    [ p, I ] = pnConvexCondMC(Nrand, x, mu, Lambdahalf, [], GammaKnown);
                    last_mu = mu;
                    last_Lambda = Lambda;
                    last_iLambda = iLambda;
                    last_I = I;
                else
                    LR = LRnormal( Y, last_mu, last_Lambda, last_iLambda, mu, Lambda, iLambda );
                    I = LR.*last_I;
                    p = mean(I);
                end
                [ phat(n), stdev(n), half(n) ] = MonteCarlo(p);
                time(n) = toc;
                efficiency(n) = 1/(time(n) * stdev(n)^2);
                fprintf('CI for Pconvex from p: %.2e +/- %.2e. \n', phat(n), half(n));
                fprintf('CI for Pconvex from 1 - p: %.2e +/- %.2e. \n', 1 - mphat(n), mhalf(n));
                fprintf('time of the iteration: %.2f. \n', time(n));
                fprintf('Efficiency (1/sigma^2t) of the iteration: %.2f. \n', efficiency(n));

            % Change of Measure -------------------------------------------
            case 2
                tic;

                % For the first initialItn and every probUpdate iterations
                % after that, get brand new set of samples.
                if (n <= initialItn) || ...
                        (mod((n - initialItn), probUpdate) == 0)
                    I = IndConvexVanillaMC(Y, x);
                    last_mu = mu;
                    last_Lambda = Lambda;
                    last_iLambda = iLambda;
                    last_I = I;
                else
                    LR = LRnormal( Y, last_mu, last_Lambda, ...
                                   last_iLambda, mu, Lambda, iLambda );
                    I = LR.*last_I;
                    LRall(:, n) = LR';
                end
                [phat(n), stdev(n), half(n)] = MonteCarlo(I);
                time(n) = toc;
                efficiency(n) = 1/(time(n) * stdev(n)^2);
                fprintf('CI for Pconvex: %.2e +/- %.2e. \n', phat(n), ...
                        half(n));
                fprintf('time of the iteration: %.2f. \n', time(n));
                fprintf('Efficiency (1/sigma^2t) of the iteration: %.2f. \n', efficiency(n));
              
            % Acceptance / Rejection --------------------------------------
            case 3
                tic;
                rej = 0;
                est = zeros(K,1);
                % For the first initialItn and every probUpdate iterations
                % after that, get brand new set of samples.
                l = mod((n - initialItn), probUpdate);
                if (n <= initialItn) || (l == 0)
                    % Record the last estimators when a complete new set of
                    % samples was obtained.
                    est = IndConvexVanillaMC(Y, x);
                    last_est = est;
                    last_Y = Y;             
                    last_mu = mu;
                    last_Lambda = Lambda;
                    last_iLambda = iLambda;
                else
                    [ratio, c] = ARParam_Normal(l, last_Y, last_mu, last_iLambda, mu, iLambda, last_Lambda, Lambda, Gamma);
                    rej = 0;
                    for k = 1:K
                        if rand() > ratio(k) / c
                            % Reject: generate a new sample and calculate the
                            % regarding indicator.
                            new_sample = mvnrnd(zeros(1,r),eye(r),1)*Lambdahalf + mu;
                            est(k) = IndConvexVanillaMC(new_sample, x);
                            rej = rej + 1;
                        else
                            % Accept: reuse the old indicator.
                            est(k) = last_est(k);
                        end
                    end
                end
                [phat(n), stdev(n), half(n)] = MonteCarlo(est);
                time(n) = toc;
                efficiency(n) = 1/(time(n) * stdev(n)^2);
                fprintf('Number of samples rejected: %d \n', rej);
                fprintf('CI for Pconvex: %.2e +/- %.2e. \n', phat(n), half(n));
                fprintf('time of the iteration: %.2f. \n', time(n));
                fprintf('Efficiency (1/sigma^2t) of the iteration: %.2f. \n', efficiency(n));        
            
        end % end of switch cases for method
        
        % Use the posterior as the prior for next itn
        mu0 = mu;
        Lambda0 = Lambda;
        iLambda0 = iLambda;
        n = n+1;
    end % end of while loop for iterations
    
%% When the sampling covariance Gamma is unknown: Normal-inv-Wishart case.
elseif GammaKnown == 0  
    last_mu = zeros(r,probUpdate)';
    last_Sigma = zeros(r,r*probUpdate)';
    last_kappa = zeros(1,probUpdate)';
    last_nu = zeros(1,probUpdate)';
    
    n = 1;
    while (n <= maxIte)
        fprintf('Iteration %d: \n', n);
        
        Yn = sampleFcn(x, s);  % new samples
        
        % update posterior
        S = zeros(r,r);
        if s~= 1
            for i = 1:s
                S = S + (Yn(i,:) - mean(Yn))'*(Yn(i,:) - mean(Yn));
            end
        end
        Yn = mean(Yn);
        kappa = kappa0 + s;
        nu = nu0 + s;
        mu = (kappa0/(kappa0+s))*mu0 + (s/(kappa0+s))*Yn;
        temp = (Yn-mu0)'*(Yn-mu0);
        Sigma = Sigma0 + S + (kappa0*s/(kappa0+s))*temp;
        % In case Sigma becomes singular due to numerical error.
        minEig = min(eig(Sigma));
        fprintf('Min eigenvalue of Sigma before adjuested: %d \n', minEig);
        while minEig < 0
            Sigma = perturbPSD(Sigma, minEig);
            minEig = min(eig(Sigma));
            fprintf('Min eigenvalue of Sigma after adjuested: %d \n', minEig);
        end
        
        % GENERATE NEW SAMPLES FROM POSTERIOR
        c = nu-r+1;
        Lambda = Sigma./(kappa*c);
        Lambdahalf = sqrtm(Lambda);
        Trand =  mvtrnd(eye(r), c, K); % generate new samples, size K-by-r
        Y = Trand*Lambdahalf + ones(K,1)*mu;
        
        % DIFFERENT METHODS
        switch method
            % Vanilla MC --------------------------------------------------
            case 0
                tic;
                Ind = IndConvexVanillaMC(Y, x);
                [ phat(n), stdev(n), half(n) ] = MonteCarlo(Ind);
                time(n) = toc;
                efficiency(n) = 1/(time(n) * stdev(n)^2);
                fprintf('CI for Pconvex: %.2e +/- %.2e. \n', phat(n), half(n));
                fprintf('time of the iteration: %.2f. \n', time(n));
                fprintf('Efficiency (1/sigma^2t) of the iteration: %.2f. \n', ...
                        efficiency(n));
        
            % Conditional MC --------------------------------------------------
            case 1
                tic;
                % the first initialItn iterations
                if (n <= initialItn) || ...
                        (mod((n - initialItn), probUpdate) == 0)
                    [ p, I ] = pnConvexCondMC( Trand, x, mu, Lambdahalf, ...
                                               nu, GammaKnown );
                    last_mu = mu;
                    last_Sigma = Sigma;
                    last_kappa = kappa;
                    last_nu = nu;
                    last_I = I;
                    [ phat(n), stdev(n), half(n)] = MonteCarlo( p );
                else
                    LR = LRstudentt( Y, last_mu, last_Sigma, last_kappa, ...
                                     last_nu, mu, Sigma, kappa, nu );
                    I = LR .* last_I;
                    [ phat(n), stdev(n), half(n)] = MonteCarlo( I );
                end

                time(n) = toc;
                efficiency(n) = 1/(time(n) * stdev(n)^2);
                fprintf('CI for Pconvex: %.2e +/- %.2e. \n', phat(n), half(n));
                fprintf('time of the iteration: %.2f. \n', time(n));
                fprintf('Efficiency (1/sigma^2t) of the iteration: %.2f. \n', efficiency(n));
        
            % Change of Measure -------------------------------------------
            case 2
                tic;
                if (n <= initialItn) || (mod((n-initialItn),probUpdate)==0)
                    I =  IndConvexVanillaMC(Y, x);
                    last_mu = mu;
                    last_Sigma = Sigma;
                    last_kappa = kappa;
                    last_nu = nu;
                    last_I = I;
                    [ phat(n), stdev(n), half(n)] = MonteCarlo( I );
                else
                    LR = LRstudentt( Y, last_mu, last_Sigma, ...
                                     last_kappa, last_nu, mu, Sigma, ...
                                     kappa, nu );
                    I = LR .* last_I;
                    [ phat(n), stdev(n), half(n)] = MonteCarlo( I );
                end
                time(n) = toc;
                efficiency(n) = 1/(time(n) * stdev(n)^2);
                fprintf('CI for Pconvex: %.2e +/- %.2e. \n', phat(n), half(n));
                fprintf('time of the iteration: %.2f. \n', time(n));
                fprintf('Efficiency (1/sigma^2t) of the iteration: %.2f. \n', efficiency(n));

            % Acceptance / Rejection --------------------------------------
            case 3
                tic;
                rej = 0;
                est = zeros(K,1);
                l = mod((n - initialItn), probUpdate);
                if (n <= initialItn) || (l == 0)
                    % Record the last estimators when a complete new set of
                    % samples was obtained.
                    est = IndConvexVanillaMC(Y, x);
                    last_est = est;
                    last_Y = Y;
                    last_mu = mu;
                    last_kappa = kappa;
                    last_nu = nu;
                    last_Sigma = Sigma;
                else                
                    [ratio, c] = ARParam_StudentT(l, last_Y, last_mu, last_kappa, last_nu, last_Sigma, mu, kappa, nu, Sigma);   
                    rej = 0;
                    for k = 1:K
                        if rand() > ratio(k) / c
                            % Reject: generate a new sample and calculate the
                            % regarding indicator.
                            new_sample = mvtrnd(eye(r), c, 1)*Lambdahalf + mu;
                            est(k) = IndConvexVanillaMC(new_sample, x);
                            rej = rej + 1;
                        else
                            % Accept: reuse the old indicator.
                            est(k) = last_est(k);
                        end
                    end
                end
                [ phat(n), stdev(n), half(n)] = MonteCarlo(est);
                time(n) = toc;
                efficiency(n) = 1/(time(n) * stdev(n)^2);
                fprintf('Number of samples rejected: %d \n', rej);
                fprintf('CI for Pconvex: %.2e +/- %.2e. \n', phat(n), half(n));
                fprintf('Time of the iteration: %.2f. \n', time(n));
                fprintf('Efficiency (1/sigma^2t) of the iteration: %.2f. \n', efficiency(n));

        end % end of switch cases for method
      
        % Use the posterior as the prior for next itn
        mu0 = mu;
        kappa0 = kappa;
        nu0 = nu;
        Sigma0 = Sigma;
        n = n+1;   
    end % end of while loop for iterations
    
end % end of if for GammaKnown

end

