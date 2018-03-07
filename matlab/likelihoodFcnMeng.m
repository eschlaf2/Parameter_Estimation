function [likelihood, window] = likelihoodFcnMeng(particles, observation, transitionFcn, W, Vth, h, b, delta)
% Calculates probabilities of each particle given observation
% Inputs:
%	particles ...				m x n array where m is the number of
%								variables and n is the number of particles.
%								m should be an even number with the first
%								half of the indices corresponding to
%								magnitudes and second half to phases
%	observation ...				Scalar measurement (number of spikes in
%								bin k with length dt)



%% Set model parameters

[~, NPARTS] = size(particles); % number of particles

%% Get V window (project backward and forward)
if ~mod(W, 2); W = W+1; end
halfWindow = floor(W/2);
k = ceil(W/2);
window = zeros([size(particles) W]);
window(:, :, k) = particles;

% project backward
for i = 1:halfWindow
	temp = window(:, :, k - i + 1);
	window(:, :, k - i) = transitionFcn(temp, true);
% 	window(:,:,k - i) = temp - F;
end
% project forward
for i = 1:halfWindow
	temp = window(:, :, k + i - 1);
	window(:, :, k + i) = transitionFcn(temp, false);
% 	window(:, :, k + i) = temp + F;
end
window = squeeze(window(1, :, :))'; % Keep only the voltage of each particle

lambda = h * ones(1, NPARTS);	% Compute lambda_k of each particle
lambda(window(1,:) >= Vth) = b;
lambda(all(window <= Vth)) = b;

likelihood = ... % (double check about multiplying by delta)
	exp(observation * log(lambda*delta) - lambda*delta);

