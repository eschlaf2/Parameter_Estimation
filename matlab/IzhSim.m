% sim

% noiseStd = [1 0 0 0]; % ... Measurment noise
delta = 0.01; % integration step [ms]
Vth = 20; % count spikes when voltage goes above Vth
TOTAL_TIME = 1e3 * 1/delta; % time steps to simulate (ms * fs)
% TOTAL_TIME = K_MAX * 1/delta;

v =	-65;
u = -10; 

p = default_Izh_params();
% p.mNoise = 0.1;
% EST = exist('estimates', 'var');
% if ~DEFAULTS
% 	for f = fieldnames(particles.params)'
% 		p.(f{:}) = particles.params.(f{:});
% 	end
% end

simParams = structfun(@(x) x * ones(1, TOTAL_TIME, 'single'), p, 'Uni', 0);
% simParams.EB = linspace(-90, -60, TOTAL_TIME);
% simParams.I = linspace(1, 5, TOTAL_TIME); 
% simParams.gB = linspace(0.5, 6, TOTAL_TIME);
% simParams.mNoise = p.mNoise * randn(1, TOTAL_TIME);

i = 1;
if exist('estimates','var')
	try
		originalSim = sim;
		for f = fieldnames(estimates.params)'
			temp = repmat(estimates.params.(f{:})(1:k)', 1/delta,1);
			simParams.(f{:})(1:k/delta) = temp(:);
		% 	simParams.(f{:})(1:k/delta) = interp1((1:1/delta:k/delta), estimates.params.(f{:})(1:k), (1:k/delta));
			figure(10)
			subplot(length(fieldnames(estimates.params)), 1, i)
			plot(simParams.(f{:}));
			i = i + 1;
		end
	catch ME
	end
end

simParams.I = simParams.I + 4*pinknoise(TOTAL_TIME);


s0 = [v; u];


%% Run Sim
sim = zeros(numel(s0), TOTAL_TIME); % to hold simulation results

sim(:, 1) = s0;
for t = 2:TOTAL_TIME
	p = structfun(@(x) x(t - 1), simParams, 'Uni', 0);
	sim(:, t) = Izh_stateTrnsn(sim(:, t-1), p, delta);
end
t = (1:TOTAL_TIME) * delta;


%% Plot Results
if ~exist('PLOT_RESULTS', 'var') || PLOT_RESULTS
    figure(100-i); fullwidth() % Voltage
    ax = subplot(8,1,2:7);
    plot(t, sim(1,:)); ylabel('Voltage [mV]'); xlabel('Time [ms]'); 
    hold on; plot([0 t(end)], Vth *  [1 1], '--', 'color', .5 * [1 1 1]); hold off
    spikes = [false logical((sim(1,1:end-1) < Vth) .* sim(1,2:end) > Vth)];
    hold on; plot(t(spikes), sim(1,spikes), 'r*'); hold off;

    subplot(8,1,1); % Spike raster
    spiketimes = find(spikes);
    spikeT = t(spikes);
    plot([spikeT(:) spikeT(:)]', [zeros(size(spikeT(:))) ones(size(spikeT(:)))]', 'color', lines(1))
    xlim(get(ax, 'xlim'))
    yticks([]); xticks([]);
	
	
end
