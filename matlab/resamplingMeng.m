function posterior = resamplingMeng(prior, likelihood, trigger, noiseStd)
% See Meng, et al., 2014
% At every spike time do a residual sampling scheme.
%	Retain M = nw copies of each particle and then supplement missing
%	particles with iid pulls from the pool of original particles with
%	probability proportional to nw - M residuals
% Otherwise bootstrap

n = size(prior, 2); % number of particles

% Update weights
weights = prior(end, :) .* likelihood + 1e-6;
% weights = weights / sum(weights);
weights = weights/sum(weights);
prior(end, :) = weights;

% Resample if triggered ...
if trigger
	M = floor(n * weights);	% copies
	p = mod(n * weights, 1);	% residuals
	
	% get M(i) copies of particle i
	inds = zeros(1, sum(M));
	k = 1; % counter for assigned particles
	for particle = 1:n	% for each particle
		for copy = 1:M(particle)	% add M(particle) copies 
			inds(k) = particle;		% ... of that particle
			k = k + 1;
		end
	end
	
	p = p / sum(p); % rescale to probability
	r = rand(1, n - sum(M));
	newParts = floor(interp1(cumsum(p), 1:n, r, 'linear', 0)) + 1;
	inds = [inds, newParts];
	
else % bootstrap
	try 
		r = rand(1, n);
		inds = floor(interp1(cumsum(weights), 1:n, r, 'linear', 0)) + 1;
	catch ME
		warning('error')
	end
end


posterior = prior(:, inds) + ... % get new particles
	noiseStd .* randn(size(prior)); % add jitter
posterior(end,:) = 1/n;

% posterior = cell2mat(arrayfun(@(i) prior(:, ...
% 	interp1(cumsum(p), prior(i, :), rand(1, n), 'linear', 0), ...
% 	(1:size(prior, 1)), 'uniformoutput', false)');
% posterior(end, :) = posterior(end, :) / sum(posterior(end, :));

end

