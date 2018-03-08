%% Set parameters

model = 'HH';	% select which model to use
SPIKETIMES = 'sim'; % simulate ('sim') or 'load' spike times

%% Load firing times

switch SPIKETIMES
	case 'load'
		load firings.mat	% sorted and curated firings from binary_small.dat (preictal); 
							% sorting by MountainSort
		spiketimes = firings(2,firings(3,:) == 95);	% spike times from unit 95
		fs = 3e4;	% sampling frequency [Hz]
	case 'sim'
		HHSim;
		fs = 1e5;
end

%% Select model
switch model
	case 'why'	% easter egg
		why
	case 'HH'	% Hodgkin-Huxley
		
		measNoise = [1 0 0 0]; % std of measurement noise (vector with meas noise for each state)
		procNoise = 0.1; % std of process noise (scalar or vector with proc noise for each param)
		dt = 0.01;	% integration step [ms]
		
		W = 5 / dt + 1;	% window [samples] (ms * fs)
		Vth = 30; % Voltage threshold [mV]
		h = 1/W; % weight
		b = h/10; % allowance
		t = ((1:W) - ceil(W/2)) * dt; % time [ms]
		
		delta = 5e-4; % binwidth [s]

		
		transitionFcn = @(particles, reverse) HH_stateTrnsn(particles, measNoise, dt, reverse);
		likelihoodFcn = @(window, obsn) ...
			likelihoodFcnMeng2011(window, obsn, t, Vth);
% 		likelihoodFcn = @(particles, obsn) normpdf(obsn - particles(1,:), 0, 1) + 1e-6; 
		resamplingFcn = @resamplingMeng;
		
		[s0, stateBounds] = HH_stateBounds(); % get initial conditions and parameter bounds
		N = 1e3; % number of particles; Meng used 1e4, but start at 1e3 for speed
		
		stateNames = {'V', 'n', 'h', 'B'};
		paramNames = {'gB', 'EB', 'VBth', 'SB', 'tauB', 'I'};
end

%% Bin spike times
binwidth = fs * delta; % samples per bin
binedges = 0:binwidth:(max(spiketimes) + binwidth);
tSpan = (1: length(binedges)-1) * delta;	% time [s]
obsn = histcounts(spiketimes, binedges);
% obsn = mean(reshape(sim(1,:)', binwidth, []));

%% Set convenience variables
NUM_STATES = size(s0, 1);
NUM_PARAMS = size(stateBounds, 1);
NUM_ALL = NUM_STATES + NUM_PARAMS + 1;
K = length(obsn);


%% Initialize PF

% Initialize parameter particles
posterior = cell2mat(arrayfun(@(a,b) unifrnd(a, b, 1, N), ...
	stateBounds(:, 1), stateBounds(:, 2), ...
	'uniformoutput', false));

% Initialize estimated states, weights, first prior
estimates = zeros(NUM_ALL, K);
bounds = zeros(NUM_ALL, K, 2); % for holding 95% bounds
estimates(1:NUM_STATES, 1) = s0;
weights = 1/N * ones(1, N);
prior = [s0 .* ones(NUM_STATES, N); posterior; weights];
procNoise = [0 * s0; procNoise .* ones(NUM_PARAMS,1); 0]; 
measNoise = [measNoise(:); zeros(NUM_ALL - numel(measNoise), 1)];

% Initialize window
% if ~mod(W, 2); W = W+1; end
% halfWindow = floor(W/2);
% t0 = ceil(W/2);
% window = zeros([size(prior) W]);


%% Run PF

for k = 1:min(K, 1e3)		% for each observation

	prediction = prior;
	for i = 1:binwidth
		prediction = transitionFcn(prediction, false);	% ... integrate states
	end
	if any(isnan(prediction))
		warning('debug')
	end
	
	% Update window
	window = updateWindow(prediction, t, transitionFcn);
	
	likelihood = likelihoodFcn(window, obsn(k)); % ... calculate likelihood
	[posterior, ~] = resamplingFcn(prediction, likelihood, 0, procNoise); % ... resample particles
		
	% Get next estimates using weighted mean
	estimates(:, k) = [sum(posterior(1:end-1, :) .* posterior(end, :), 2); mean(posterior(end, :))];
	temp = sort(posterior, 2);
	bounds(:, k, :) = temp(:, [floor(N*0.025) ceil(N*0.975)]);
	

	prior = posterior; % ... get the next prior 
	
	if k == 100
		disp('look around')
	end
	
	if 0
		figure(999)
		x = 1; y = 1;
		scatter(posterior(NUM_STATES + x, :), posterior(NUM_STATES + y, :), 16, posterior(end,:), 'filled');
		hold on; plot(sim(x + NUM_STATES, k * binwidth), sim(y + NUM_STATES, k*binwidth), 'r*'); hold off;
		hold on; plot(estimates(x + NUM_STATES, k), estimates(y + NUM_STATES, k), 'b*'); hold off;
		xlim(stateBounds(x,:)); ylim(stateBounds(y,:));
		colormap('cool'); colorbar
		title(sprintf('%d: %d', k, obsn(k)))
		drawnow;
		pause(1e-6)
	end
end

%% Plot results
colors = lines(7);

figure(3); fullwidth()
stem(tSpan(1:k), obsn(1:k)'); hold on;
plot(tSpan(1:k), estimates(1, 1:k)/Vth, 'Color', colors(2,:));
plot(tSpan(1:k), sim(1, 1:binwidth:k*binwidth)/Vth, 'Color', colors(1,:))
plot(tSpan(1:k), ones(1, k), '--', 'Color', .5*[1 1 1]); hold off;
ylim(1.1*(get(gca,'YLim') - mean(get(gca, 'Ylim'))) + mean(get(gca, 'YLim')));
legend('Spikes','Estimated Voltage', 'True Voltage', 'location', 'best');

figure(4); fullwidth()
plot(estimates(2:4, 1:k)', 'linewidth', 2); legend(stateNames(2:4))
hold on; set(gca, 'ColorOrderIndex', 1)
plot(sim(2:4, 1:binwidth:k*binwidth)', 'o'); hold off;

figure(5); fullwidth()
dummy = NUM_STATES + 1;
names = [stateNames(:); paramNames(:)];
shape = [ [1:k, k:-1:1]; [bounds(dummy:end-1, 1:k, 1), bounds(dummy:end-1, k:-1:1, 2)] ];
plot(estimates(dummy:end - 1, 1:k)', 'linewidth', 2, 'displayname', sprintf('%s Est', names{dummy:end})); 
patch(shape(1,:), shape(2,:), 'b', 'FaceColor', colors(1,:), 'facealpha', .5, 'edgecolor', 'none', 'displayname', '95CI');
hold on; set(gca, 'ColorOrderIndex', 1)
plot(sim(dummy:end, 1:binwidth:k*binwidth)', '--', 'displayname', sprintf('%s Truth', names{dummy:end})); 
legend(); hold off;


figure(6); fullwidth()
plot(estimates(end, 1:k)); title('Mean Likelihood');


%% Supplementary functions

function window = updateWindow(prediction, t, transitionFcn)

t0 = find(t == 0);
window = zeros([size(prediction) numel(t)]);
fwd = numel(t) - t0;
bkwd = t0 - 1;

window(:, :, t0) = prediction;
	for i = 1:min(fwd, bkwd)
		window(:, :, t0 - i) = transitionFcn(window(:, :, t0 - i + 1), true); % ... project backward
		window(:, :, t0 + i) = transitionFcn(window(:, :, t0 + i - 1), false); % ... project forward
	end
	for i = min(fwd, bkwd):max(fwd, bkwd)
		sgn = -1^(bkwd > fwd);
		window(:, :, t0 + sgn * i) = transitionFcn(window(:, :, t0 + sgn * (i - 1)), sgn == -1);
	end
		
end
