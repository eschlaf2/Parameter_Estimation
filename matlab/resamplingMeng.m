function [posterior, inds] = resamplingMeng(prior, likelihood, trigger, noiseStd)
% See Meng, et al., 2014
% At every spike time do a residual sampling scheme.
%	Retain M = nw copies of each particle and then supplement missing
%	particles with iid pulls from the pool of original particles with
%	probability proportional to nw - M residuals
% Otherwise bootstrap

N = size(prior, 2); % number of particles
inds = noiseStd > 0;

% Update weights
weights = prior(end, :) .* likelihood + 1e-6;
weights = weights/sum(weights);
prior(end, :) = weights;

cc = corrcoef(prior([inds(1:end-1); true], :)');
ll_dist = sort(histcounts(likelihood, 10), 'descend');
if ll_dist(1) / N > .99
	rho = 1.1;
else
	rho = .999 - .05 * abs(cc(2:end, 1)); % discount factor
end


% Draw new parameters
theta = prior(noiseStd > 0, :);
m = rho .* theta + (1 - rho) .* sum(weights .* theta, 2);
h2 = 1 - rho.^2;
sigma = std(theta, [], 2);
noiseStd(noiseStd > 0) = max(h2 .* sigma, noiseStd(noiseStd > 0));
prior(noiseStd > 0, :) = m;

% Resample if triggered ...
if trigger
	M = floor(N * weights);	% copies
	p = mod(N * weights, 1);	% residuals
	
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
	prior(end,:) = 1/N;
	
else % bootstrap
	try 
		r = rand(1, N);
		inds = floor(interp1(cumsum(weights), 1:N, r, 'linear', 0)) + 1;
	catch ME
		warning('error')
	end
end


posterior = prior(:, inds) + ... % get new particles
	noiseStd .* randn(size(prior)); % add jitter
posterior(end, :) = posterior(end, :) / sum(posterior(end, :));

end

