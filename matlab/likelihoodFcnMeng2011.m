function [likelihood, window] = ...
	likelihoodFcnMeng2011(particles, observation, window, W, Vth, h, b, delta)
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

window = squeeze(window(1, :, :))'; % Keep only the voltage of each particle

lambda = sum(g(window, Vth) .* f(((1:W)' - k)*delta));
likelihood = exp(observation * log(lambda*delta) - lambda*delta);
likelihood(isnan(likelihood)) = 0;

end

function y = f(x)
% p = 0.9; 
y = 0.9.^ abs(x);

end

function y = g(x, Vth)
eta = 1.622; %
nu = 0.1;

y = eta * exp(nu * (x - Vth))./(1 + exp(nu * (x - Vth)));
end

