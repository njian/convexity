%----------------********************************--------------------------
% sampleFcn.m
% Description:
%   Put the black box function that you want to sample here. Make sure it
%   returns all the required outputs.
% Inputs:
%   x: the design points to sample on.
%   num: the number of samples to obtain.
% Outputs:
%   y: num-by-r matrix, each row is a set of samples on x.
%   Gamma: the r-by-r sampling covariance matrix.
%   d: the dimension of x space.
%   xlb: d-by-1 lower bound on the x space.
%   xub: d-by-1 upper bound on the x space.
%   GammaKnown: 0 for assuming unknown sampling variance, and 1 otherwise.
%----------------********************************--------------------------
function [ y, Gamma, d, xlb, xub, GammaKnown ] = sampleFcn( x, num )
%% Sample (x) space:
y = []; Gamma = [];
d = 1; % dimension of the X space
bound = 1;
xub = bound * ones(d,1); % upper bounds for x
xlb = -bound * ones(d,1); % lower bounds for x
GammaKnown = 0; % assume known (1) or unknown (0) variance
% if no design points is provided, only return the basic information.
if isempty(x)
    fprintf('Sample function information: d = %d, GammaKnown = %d. \n', d, GammaKnown)
    return
end

%% Now construct the sampling covariance matrix
[ d, r ] = size(x);
% B controls how large the correlation between points is (compared to the
% ones on the diagonal)
B = zeros(r);
for i = 1:r
    for j = 1:i
        % constructs correlation using Gaussian radial basis function
        B(i,j) = exp(- norm(x(:,i) - x(:,j))^2 / 2);
    end
end
B = 0.01 * (tril(B, -1) + tril(B, -1)');
% Gamma constructed by ones on the diagonal and B off-diagonal, scaled by
% scaleGamma
scaleGamma = 1e-2;
Gamma = scaleGamma * (eye(r) + B);

%% Test functions
%% Strictly convex: sphere (square norm)
y = zeros(num, r);
for i = 1:r
    y(:, i) = sum(x(:, i) .^ 2) * ones(num, 1);
end
% Gamma constructed by ones on the diagonal and B off-diagonal, with (i,j)
% th component scaled by the size of y(i)y(j) times a scaling factor
scaleGamma = 0.04*sqrt(y' * y); % like y
Gamma = scaleGamma .* (eye(r) + B);
y = y + mvnrnd(zeros(num, r), Gamma);

%% Strictly concave: -sphere
% y = zeros(num, r);
% for i = 1:r
%     y(:, i) = - (1/(bound^2*d)) * sum(x(:, i) .^ 2) * ones(num, 1);
% end
% y = y + mvnrnd(zeros(num, r), Gamma);

%% Slightly convex: -logarithm
% y = zeros(num, r);
% for i = 1:r
%     y(:,i) = -sum(log(x(:,i))) * ones(num, 1);
% end
% y = y + mvnrnd(zeros(num,r),Gamma);
%
%% Slightly concave: logarithm
% y = zeros(num, r);
% for i = 1:r
%     y(:,i) = sum(log(x(:,i))) * ones(num, 1);
% end
% y = y + mvnrnd(zeros(num,r),Gamma);

%% In between: plane
% y = zeros(num,r);
% y = y + mvnrnd(zeros(num, r),Gamma);

%% Special non-convex: "w"-shaped
% y = zeros(num, r);
% for i = 1:r
%     y(:,i) = log(14 * sum(x(:,i).^2,1) / d + 10 * exp( - 2 * sum(x(:,i).^2,1) / d ));
% end
% y = y + mvnrnd(zeros(num, r), Gamma);

%% Ambulance
% y = zeros(num,r);
% i = 1;
% for n = 1:num
%     while i <= r
%         % seed = ceil(1000*rand());
%         seed = 1;
%         y(n,i) = Ambulance(x(:,i)', 720, seed, []);
%         if isnan(y(n,i))~=1
%             i = i+1;
%         end
%     end
%     i=1;
% end

end

