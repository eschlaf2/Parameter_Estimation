function [posterior, inds] = resamplingMeng(particles, likelihood, obs)
% See Meng, et al., 2014
% At every spike time do a residual sampling scheme.
%	Retain M = nw copies of each particle and then supplement missing
%	particles with iid pulls from the pool of original particles with
%	probability proportional to nw - M residuals
% Otherwise bootstrap

N = length(particles.weights);
fn = fieldnames(particles.params);
wDist = sort(particles.weights, 'descend');
% trigger = sum(wDist(1:N/10)) > .6;
trigger = obs;

% Update weights
weights = particles.weights .* likelihood + 1e-6;
weights = weights/sum(weights);

idx = @(A, ind) A(ind);
cc = structfun(@(x) idx(corrcoef(x, weights), 3), particles.params);
% ll_dist = sort(histcounts(likelihood, 10), 'descend');
% if ll_dist(1) / N > .98
% 	rho = 1.01 * ones(size(cc));
% else
% 	rho = 1.01 - .08 * abs(cc); % discount factor
	rho = ones(size(cc));
	
% end


% Draw new parameters
m = @(theta, rho) rho .* theta + (1 - rho) .* sum(weights .* theta);
h2 = 1 - rho.^2;
sigma = structfun(@std, particles.params);
% noiseStd(noiseStd > 0) = max(h2 .* sigma, noiseStd(noiseStd > 0));
for i = 1:length(fn)
	particles.params.(fn{i}) = m(particles.params.(fn{i}), rho(i));
end
% particles.params = structfun(m, particles.params, 'Uni', 0);
particles.pNoise = max(h2 .* sigma, particles.pNoise);

% Resample if triggered ...
if trigger
	M = floor(N * weights);	 % copies
	p = mod(N * weights, 1);  % residuals
	
	% get M(i) copies of particle i
	inds = zeros(1, sum(M));
	k = 1; % counter for assigned particles
	for particle = 1:N	% for each particle
		for copy = 1:M(particle)	% add M(particle) copies 
			inds(k) = particle;		% ... of that particle
			k = k + 1;
		end
	end
	
	p = p / sum(p); % rescale to probability
	r = rand(1, N - sum(M));
	newParts = floor(interp1(cumsum(p), 1:N, r, 'linear', 0)) + 1;
	inds = [inds, newParts];
	particles.weights = 1/N * ones(1, N);
	
else % bootstrap
	try 
		r = rand(1, N);
		inds = floor(interp1(cumsum(weights), 1:N, r, 'linear', 0)) + 1;
		particles.weights = weights;
	catch ME
		warning('error')
	end
end

posterior = particles;
posterior.params = structfun(@(x) x(inds), particles.params, 'Uni', 0);
for i = 1:length(fn)
	f = fn{i};
	posterior.params.(f) = posterior.params.(f) + ...
		posterior.pNoise(i) * randn(1, N);
end
posterior.weights = posterior.weights(inds) / sum(posterior.weights(inds));
posterior.states = posterior.states(:, inds);
posterior.trigger = trigger;

end

