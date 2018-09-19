function [likelihood, window] = likelihoodFcnMeng(window, observation, W, Vth, delta)
% Calculates probabilities of each particle given observation
% Inputs:
%	window ...					m x n x k array where m is the number of
%								variables, n is the number of particles and
%								k is the time span of the window (in
%								samples).
%	observation ...				Indicator measurement (is there a spike in the bin)
%	t ...						k x 1 vector of times
%	Vth ...						scalar voltage threshold [mV]

%% Set parameters
h = 1/W; % weight
b = h/100; % allowance
N = size(window, 2); % number of particles

v = squeeze(window(1, :, :))'; % Keep only the voltage of each particle
% u = squeeze(window(2, :, :))';
crossings = sum(diff(v > Vth) > 0);
% crossings = any(window > Vth);
lambda = h * 1./(abs(crossings - observation).^2 + 1);

% lambda = h * ones(1, N);	% Compute lambda_k of each particle
% lambda(window(1,:) >= Vth) = b;
lambda(all(v <= Vth)) = b;

% lambda = h * ones(1, N);
% lambda = h * (crossings == observation);
% if ~observation
% 	lambda(crossings == 0) = b;
% end
% lambda = b * ones(1, N);
% lambda(abs(crossings - observation) < 2) = h;

likelihood = exp(observation * log(lambda*delta) - lambda*delta);
% likelihood(any(v < -200)) = 0;
% likelihood(any(u < -20)) = 0;

