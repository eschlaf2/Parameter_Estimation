function particles = HH_stateTrnsn(particles, noiseStd, delta, reverse)
% Euler integration with noise

default_noiseStd = [1 0 0 0];
default_delta = 0.01;

bounds = [-Inf Inf; ... % V
	0 1; ... % n
	0 1; ... % h
	0 1; ... % B
	];

if ~exist('noiseStd', 'var') || isempty(noiseStd)
	noiseStd = default_noiseStd;
end
if ~exist('delta', 'var') || isempty(delta)
	delta = default_delta;
end
if ~exist('reverse', 'var') || isempty(reverse)
	reverse = false;
end

%% Set parameters

NOISE = zeros(size(particles, 1), 1);
NOISE(1:length(noiseStd)) = noiseStd;

F = HH_dynamics(particles);
if ~reverse
	particles = particles + delta * F + NOISE .* randn(size(particles));
else
	particles = particles - delta * F + NOISE .* randn(size(particles));
end

for i = 1:size(bounds, 1)
	particles(i, particles(i, :) > bounds(i, 2)) = bounds(i, 2);
	particles(i, particles(i, :) < bounds(i, 1)) = bounds(i, 1);
end
