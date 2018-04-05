function particles = HH_stateTrnsn(particles, measNoise, delta, noiseParam)
% Euler integration with noise

% noiseParam = strcmp('mNoise', paramNames);
noiseParam = 11;
default_measNoise = [1 0 0 0];
% default_delta = 0.01;

% bounds = [ ...
% 	-Inf Inf; ... % V
% 	0 1; ... % n
% 	0 1; ... % h
% 	0 1; ... % B
% 	];

if ~exist('measNoise', 'var') || isempty(measNoise)
	measNoise = default_measNoise;
end
% if ~exist('delta', 'var') || isempty(delta)
% 	delta = default_delta;
% end

%% Set parameters

NOISE = particles(noiseParam, :);

dF = HH_dynamics([], particles);

big_slopes = max(abs(dF)) > 5000;
if any(big_slopes)
	options = odeset('InitialStep', delta/5);
	for p = find(big_slopes)
% 		display(['big_slopes at t = ', num2str(t)]);
		[~, temp] = ode45(@HH_dynamics, [0 delta], particles(:, p), options);
		dF(:, p) = temp(end, :)' - particles(:, p);
	end
end

particles = particles + delta * dF;
particles(1, :) = particles(1, :) + NOISE .* randn(1, size(particles, 2));
% for i = 1:size(bounds, 1)
% 	particles(i, particles(i, :) > bounds(i, 2)) = bounds(i, 2);
% 	particles(i, particles(i, :) < bounds(i, 1)) = bounds(i, 1);
% end
