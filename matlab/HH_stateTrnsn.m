function particles = HH_stateTrnsn(particles, t, delta, paramNames)
% Euler integration with noise

% noiseParam = strcmp('mNoise', paramNames);
% default_noiseStd = [1 0 0 0];
% default_delta = 0.01;

bounds = [ ...
	-Inf Inf; ... % V
	0 1; ... % n
	0 1; ... % h
	0 1; ... % B
	];

% if ~exist('noiseStd', 'var') || isempty(noiseStd)
% 	noiseStd = default_noiseStd;
% end
% if ~exist('delta', 'var') || isempty(delta)
% 	delta = default_delta;
% end

%% Set parameters

NOISE = zeros(size(particles));
% NOISE(1:length(noiseStd), :) = noiseStd(:) .* particles(noiseParam, :);

dF = HH_dynamics([], particles);

big_slopes = max(abs(dF)) > 5000;
if any(big_slopes)
	options = odeset('InitialStep', delta/5);
	for p = find(big_slopes)
% 		display(['big_slopes at t = ', num2str(t)]);
		[~, temp] = ode45(@HH_dynamics, [0 delta], particles(:, p));
		dF(:, p) = temp(end, :)' - particles(:, p);
	end
end
particles = particles + delta * dF;
% particles = particles + delta * F + NOISE .* randn(size(particles));
% for i = 1:size(bounds, 1)
% 	particles(i, particles(i, :) > bounds(i, 2)) = bounds(i, 2);
% 	particles(i, particles(i, :) < bounds(i, 1)) = bounds(i, 1);
% end
