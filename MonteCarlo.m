function [ avg, stdev, half ] = MonteCarlo( y )
% Standard Monte Carlo function that returns the sample average, standard
% deviation, and the half width of a 95% CI
n = length(y);
avg = mean(y);
stdev = std(y);
half = 1.96*std(y)/sqrt(n);
end

