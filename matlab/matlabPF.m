%% NOT WORKING... DON'T USE THIS!!!

%% Set parameters

model = 'HH';	% select which model to use

%% Select model
switch model
	case 'why'	% easter egg
		why
	case 'HH'	% Hodgkin-Huxley
		ops = set_options();
		transitionFcn = @(particles) HH_dynamics(particles, ops);
		stateBounds = HH_stateBounds();
		numParticles = 1e3; % Meng used 1e4, but start at 1e3 for speed
		
		lgnd = {...
			'V', 'n', 'h', 'B', 'gB', 'EB', 'VBth', 'SB', 'tauB', 'I'};
end

%% Load firing times
load firings.mat	% sorted and curated firings from binary_small.dat (preictal); 
					% sorting by MountainSort
spiketimes = firings(2,firings(3,:) == 95);	% spike times from unit 95
fs = 3e4;	% sampling frequency [Hz]
delta = str2double(ops.delta);
binwidth = fs * delta; % samples per bin
binedges = 0:binwidth:(max(spiketimes) + binwidth);
tSpan = (1: length(binedges)-1) * delta;	% time [s]
observations = histcounts(spiketimes, binedges);


%% Run filter
myPF = particleFilter(...
	@(particles) transitionFcn(particles), ...
	@(particles, measurement) likelihoodFcnMeng(particles, measurement, ops));
									

initialize(myPF, numParticles, stateBounds);
myPF.ResamplingMethod = 'systematic';
myPF.ResamplingPolicy.MinEffectiveParticleRatio = 0.5;

estimate = zeros(size(observations, 1), myPF.NumStateVariables);
for k=1:1e3
    estimate(k,:) = correct(myPF,observations(k));
    predict(myPF);
	
	if 0
		figure(99);
		plot(k * ones(size(myPF.Weights)), myPF.Weights, 'g.'); hold on;
		xlim([0 1e3]);
	% 	histogram(myPF.Weights(1,:), 30, 'normalization', 'probability')
	% 	ylim([0 1]); title(num2str(k)); % xlim([-80 60]); 
		drawnow;
	end
end

%% Plot results
figure(3); fullwidth()
plot(tSpan(1:k), observations(1:k)','b',...
	tSpan(1:k), estimate(1:k,1)/str2double(ops.Vth),'r'); hold on;
plot(tSpan(1:k), ones(1, k), '--', 'Color', .5*[1 1 1]); hold off;
ylim(1.1*(get(gca,'YLim') - mean(get(gca, 'Ylim'))) + mean(get(gca, 'YLim')));
legend('Spikes','Estimated Voltage', 'location', 'best');

figure(4); fullwidth()
plot(estimate(:, 2:4)); legend(lgnd(2:4))

figure(5); fullwidth()
plot(estimate(:, 5:end)); legend(lgnd(5:end))

