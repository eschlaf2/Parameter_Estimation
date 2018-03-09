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
		
		noiseStd = [1 0 0 0]; % std of measurement noise [mV]
		dt = 0.01;	% integration step [ms]
		
		W = 5 / dt + 1;	% window [samples] (ms * fs)
		Vth = 30; % Voltage threshold [mV]
		h = 1/W; % weight
		b = h/10; % allowance
		
		delta = 5e-4; % binwidth [s]

		
		transitionFcn = @(particles, reverse) HH_stateTrnsn(particles, noiseStd, dt, reverse);
		likelihoodFcn = @(particles, obsn) ...
			likelihoodFcnMeng2011(particles, obsn, transitionFcn, W, Vth, h, b, dt);
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


%% Run PF

% Initialize parameter particles
posterior = cell2mat(arrayfun(@(a,b) unifrnd(a, b, 1, N), ...
	stateBounds(:, 1), stateBounds(:, 2), ...
	'uniformoutput', false));

% Initialize estimated states, weights, first prior, window
estimates = zeros(NUM_ALL, K);
estimates(1:NUM_STATES, 1) = s0;
weights = 1/N * ones(1, N);
prior = [s0 .* ones(NUM_STATES, N); posterior; weights];
measNoise = [0 * s0; .1 * ones(NUM_PARAMS,1); 0]; 

for k = 2:min(length(obsn), 1e3)		% for each observation

	prediction = prior;
	for i = dt:dt:delta*1e3
		prediction = transitionFcn(prediction, false);	% ... integrate states
	end
	if any(isnan(prediction))
		warning('debug')
	end
	likelihood = likelihoodFcn(prediction, obsn(k)); % ... calculate likelihood
	posterior = resamplingFcn(prediction, likelihood, 0, measNoise); % ... resample particles
	p = [sum(posterior(1:end-1, :) .* posterior(end, :), 2); mean(posterior(end, :))];

	estimates(:, k) = [sum(posterior(1:end-1, :) .* posterior(end, :), 2); mean(posterior(end,:))];

	prior = posterior; % ... get the next prior 
	
	if k == 100
		disp('look around')
	end
	
	if 1
		figure(999)
		x = 1; y = 1;
		scatter(posterior(NUM_STATES + x, :), posterior(NUM_STATES + y, :), 16, posterior(end,:), 'filled');
		hold on; plot(sim(x + NUM_STATES, k * binwidth), sim(y + NUM_STATES, k*binwidth), 'r*'); hold off;
		hold on; plot(estimates(x + NUM_STATES, k), estimates(y + NUM_STATES, k), 'b*'); hold off;
		xlim(stateBounds(x,:)); ylim(stateBounds(y,:));
		colormap('cool'); caxis([0 1/N]); colorbar
		title(sprintf('%d: %d', k, obsn(k)))
		drawnow;
		pause(1e-6)
	end
end

%% Plot results
figure(3); fullwidth()
stem(tSpan(1:k), obsn(1:k)'); hold on;
plot(tSpan(1:k), estimates(1, 1:k)/Vth,'r');
plot(tSpan(1:k), sim(1, 1:binwidth:k*binwidth)/Vth, 'b')
plot(tSpan(1:k), ones(1, k), '--', 'Color', .5*[1 1 1]); hold off;
ylim(1.1*(get(gca,'YLim') - mean(get(gca, 'Ylim'))) + mean(get(gca, 'YLim')));
legend('Spikes','Estimated Voltage', 'location', 'best');

figure(4); fullwidth()
plot(estimates(2:4, 1:k)', 'linewidth', 2); legend(stateNames(2:4))
hold on; set(gca, 'ColorOrderIndex', 1)
plot(sim(2:4, 1:binwidth:k*binwidth)', 'o'); hold off;

figure(5); fullwidth()
plot(estimates(5:end - 1, 1:k)', 'linewidth', 2); legend(paramNames)
hold on; set(gca, 'ColorOrderIndex', 1)
plot(sim(5:end, 1:binwidth:k*binwidth)', 'o'); hold off;

figure(6); fullwidth()
plot(estimates(end, 1:k)); legend('weights');
