% sim

% noiseStd = [1 0 0 0]; % ... Measurment noise
delta = 0.01; % integration step [ms]
Vth = 30; % count spikes when voltage goes above Vth
TOTAL_TIME = 1e3 * 1/delta; % time steps to simulate (ms * fs)
% TOTAL_TIME = K_MAX * 1/delta;

V =	-71;
n = 0.0147; 
h = 0.7497;
B = 0.0326;

p = default_HH_params();

p.mNoise = 0.5;

simParams = structfun(@(x) x * ones(1, TOTAL_TIME, 'single'), p, 'Uni', 0);


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

else
	% simParams.EB = linspace(-90, -60, TOTAL_TIME);
% 	simParams.I = linspace(1, 5, TOTAL_TIME); 
	% simParams.gB = linspace(0.5, 6, TOTAL_TIME);
end


simParams.I = simParams.I + 8*pinknoise(TOTAL_TIME);
% simParams.mNoise = p.mNoise * randn(1, TOTAL_TIME);

s0 = [V; n; h; B];


%% Run Sim
sim = zeros(numel(s0), TOTAL_TIME); % to hold simulation results

sim(:, 1) = s0;
for t = 2:TOTAL_TIME
	p = structfun(@(x) x(t - 1), simParams, 'Uni', 0);
	sim(:, t) = HH_stateTrnsn(sim(:, t-1), p, delta);
end
t = (1:TOTAL_TIME) * delta;


%% Plot Results
if ~exist('PLOT_RESULTS', 'var') || PLOT_RESULTS
    figure(99); fullwidth() % Voltage
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
