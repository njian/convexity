%----------------********************************--------------------------
% MAIN.m Description:
%   Setting the file name of the test function and the parameters of the
%   algorithm.
% Calls:
%   sequentialConvex.m
% Author: Nanjing Jian, Cornell University. 
% Reference:
%   Jian and Henderson, 2016, "Estimating the Probability that a Function
%   Observed with Noise is Convex", submitted to INFORMS Journal of
%   Computing.
%----------------********************************--------------------------
clc;clear;close all;

%% INPUTS
% The .m file name that contains the test function.
sampleFcn = 'sampleFcn'; 

% Method to use: choose from 'vanillaMC' or 'conditionalMC' or 
% 'changeMeasure' or 'acceptanceRejection'. See reference for details.
prompt = 'Please choose a method by entering a number: \nAvailable options are (default = 0)\n 0 - vanillaMC \n 1 - conditionalMC \n 2 - changeMeasure \n 3 - acceptanceRejection \n';
method = input(prompt);
checkPrompt(method, 'method');

% Number of posterior updates.
prompt = 'Please enter the number of iterations of posterior updates: \n (integer >= 1, default = 100) \n';
maxIte = input(prompt);
checkPrompt(maxIte, 'positiveInteger');

% Number of iterations of Monte Carlo simulation to estimate the posterior
% probability of convexity.
prompt = 'Please enter the number of Monte Carlo samples from the posterior distribution to estimate the probability: \n (integer >= 1, default = 100) \n';
K = input(prompt);
checkPrompt(K, '>=1Integer');

% random seed for sampling the function
prompt = 'Please choose a random seed: \n (integer >= 1, default = 1) \n';
seed = input(prompt);
checkPrompt(seed, '>=1Integer');

%% RUN ALGORITHM
fileName = datestr(clock, 'mm_dd_yyyy_HH_MM');
diary([pwd, '\outputs\', fileName, '.txt']);
[ phat, half, x, mu, Lambda, time, efficiency ] = sequentialConvex( ...
    sampleFcn, maxIte, K, seed, method, fileName);
diary; % turn off diary mode