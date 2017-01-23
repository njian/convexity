%----------------********************************--------------------------
% SequantialConvex.m
% Description: 
%   Gives an estimation to the probability that the input function is 
%   convex.
%
% Calls:
%   algorithm1.m, sampleFcn.m, perturbPSD.m, plotOutputs.m
%
% Inputs:
%   sampleFcn: the function to evaluate convexity of.
%   nSampleLines: the number of lines in space to sample the design 
%                 points on.
%   maxIte: a preset number of iterations of the posterior update on the
%           function values.
%   K: the number of Monte Carlo simulation iterations for probability
%      estimation.
%   seed: an integer that's  used as the random seed.
%   method: 'vanillaMC' or 'conditionalMC' or 'changeMeasure' or
%   'acceptanceRejection'. See reference for details.
%   fileName: save the data, plot, and console prints to this fileName.
%
% Outputs:
%   phat: the estimated probability of convex along maxIte iterations.
%   half: the estimated CI halfwidth for phat.
%   x: the sampled design points.
%   mu: the posterior mean of the function.
%   Lambda: the posterior covariance matrix of the function.
%   time: time spent in estimating phat via Monte Carlo in each iteration.
%   efficiency: 1/(variance*time) of the Monte Carlo estimator for phat.
%
% Author: Nanjing Jian, Cornell University.
% TODO: This implementation in Matlab is not fully "sequential" as it does 
% not allow user to stop (ctrl+c) at any time and plot the current estimate 
% and resume afterwards. Instead it uses a preset max number of iterations.
% It would be ideal if the future version of Matlab allows to catch ctrl+c
% as an error and output the current estimate.
%----------------********************************--------------------------
function [ phat, half, x, mu, Lambda, time, efficiency ] = ...
            sequentialConvex( sampleFcn, maxIte, K, seed, method, fileName)
%% -------------------------PROCESS INPUTS---------------------------------
% Information of the sample function.
[ ~, ~, d, xlb, xub, GammaKnown ] = feval(sampleFcn, [], []);

% Random streams: using Pierre L'Ecuyer's RngStream 
% Reference: https://www.iro.umontreal.ca/~lecuyer/
[Stream1, Stream2] = RandStream.create('mrg32k3a', 'NumStreams', 2);
Stream1.Substream = seed;
Stream2.Substream = seed;
OldStream = RandStream.setGlobalStream(Stream1);

%% -------------------------SAMPLE POINTS----------------------------------
% Sample 3 points along each random line in space.
nSampleLines = d + 1; % Number of such lines = dimension of x space + 1.
r = 3 * nSampleLines;
x = zeros(d, r);
for i = 1:nSampleLines
    D = mvnrnd(zeros(d, 1),eye(d))';
    D = D ./ norm(D, 2);
    center = xlb + rand(d,1) .* (xub - xlb);
    maxt = zeros(d, 1);
    mint = zeros(d, 1);
    for j = 1:d
        if D(j)>=0
            maxt(j) = (xub(j) - center(j)) / D(j);
            mint(j) = (xlb(j) - center(j)) / D(j);
        else
            maxt(j) = (xlb(j) - center(j)) / D(j);
            mint(j) = (xub(j) - center(j)) / D(j);
        end
        tmax = min(maxt);
        tmin = max(mint);
    end
    x(:, 3 * i - 0) = center + (tmin + rand() * (tmax-tmin)) * D;
    x(:, 3 * i - 2) = center + (tmin + rand() * (tmax-tmin)) * D;
    x(:, 3 * i - 1) = center + (tmin + rand() * (tmax-tmin)) * D;
end

% Uniform random sampling
% r = 100;
% for i = 1:r
%     x(:, i) = xlb + rand(d,1) .* (xub - xlb);
% end

% Sampling on a grid
% x = zeros(d, r);
% k = 1;
% for i = 1:sqrt(r)
%     for j = 1:sqrt(r)
%         x(:,k) = xlb + [i/sqrt(r); j/sqrt(r)];
%         k = k + 1;
%     end
% end

%% -------------------------PRIOR DISTRIBUTION-----------------------------
Gamma = [];
if GammaKnown == 1
    % Uninformative: constant mean, huge variance
    [~, Gamma] = feval(sampleFcn, x, 1);
    mu0 = zeros(1, r);
    Lambda0 = 1e4 * eye(r);
    kappa0 = r + 1;
    nu0 = r;
    Sigma0 = [];
else
    % Jeffery's prior
    m0 = r + 1;
    Y = feval(sampleFcn, x, m0);
    kappa0 = m0;
    nu0 = m0 - 1;
    mu0 = mean(Y);
    Sigma0 = zeros(r,r);
    for i = 1:m0
        Sigma0 = Sigma0 + (Y(i,:)-mu0)' * (Y(i,:)-mu0);
    end
    % If Sigma0 is close to singular, add a small number of the diagonal.
    Sigma0 = perturbPSD(Sigma0, min(eig(Sigma0)));
    Lambda0 = [];
end


%% -------------------------RUN ALGORITHM----------------------------------
% Run algorithm
[ phat, half, mu, Lambda, time, efficiency ] = algorithm( ...
    x, mu0, Lambda0, kappa0, nu0, Sigma0, Gamma, ...
    maxIte, K, GammaKnown, method );
% Save a snapshot of the current workspace
save([pwd, '\outputs\', fileName]);
fprintf('========================SUMMARY================================\n')
fprintf('After %d iterations (samples) of posterior updates, \n using %d samples per iteration in the Monte Carlo simulation to estimate Pconvex, \n the estimated CI for Pconvex is %.2e +/- %.2e. \n',maxIte,K,phat(maxIte),(half(maxIte)));
fprintf('Average time per iteration: %d. \n', mean(time));
fprintf('Average efficiency per iteration (over finite values): %d. \n', mean(efficiency(isfinite(efficiency))));

% Plot output
plotOutputs( phat, half, time, efficiency, fileName );
fprintf('Figure saved to %s_phat/time/efficiency.fig in the current directory.\n', fileName);

% Restore old stream
RandStream.setGlobalStream(OldStream);
end




