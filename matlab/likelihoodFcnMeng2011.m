function [likelihood, window] = ...
	likelihoodFcnMeng2011(window, observation, t, Vth)
% Calculates probabilities of each particle given observation
% Inputs:
%	window ...					m x n x k array where m is the number of
%								variables, n is the number of particles and
%								k is the time span of the window (in
%								samples).
%	observation ...				Indicator measurement (is there a spike in the bin)
%	t ...						k x 1 vector of times
%	Vth ...						scalar voltage threshold [mV]



%% Set model parameters

delta = t(2) - t(1);

window = squeeze(window(1, :, :))'; % Keep only the voltage of each particle

lambda = sum(g(window, Vth) .* f(t(:)));
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

