function [sim, spiketimes, simParams] = modelSim(model, Vth, estimates, binwidth, varargin)

if ~exist('model', 'var') || isempty(model)
	model = 'HH';
end
if ~exist('Vth', 'var') || isempty(Vth)
	Vth = 30;
end

dt = 1e-2;
TOTAL_TIME = 1e3 * 1/dt; % time steps to simulate (ms * fs)

for a = 1:2:length(varargin)
	switch lower(varargin{a})
		case 'dt'
			dt = varargin{a+1};
		case 'total_time'
			TOTAL_TIME = varargin{a+1};
	end
end

p = default_params(model);

switch model
	case 'Izh'
		v =	-65;
		u = -10; 
		s0 = [v; u];
		
		noise = 4 * pinknoise(TOTAL_TIME);
		transitionFcn = @(state, p) Izh_stateTrnsn(state, p, dt);

		
	case 'HH'
		V =	-71;
		n = 0.0147; 
		h = 0.7497;
		B = 0.0326;
		s0 = [V; n; h; B];
		
		noise = 8 * pinknoise(TOTAL_TIME);
		noise = 0;
		transitionFcn = @(state, p) HH_stateTrnsn(state, p, dt);
		
end

% Create time series for each paramter
simParams = structfun(@(x) x * ones(1, TOTAL_TIME, 'single'), p, 'Uni', 0);

% Use estimated parameters if given
i = 1;
if exist('estimates','var') && ~isempty(estimates)
	if ~exist('binwidth', 'var') || isempty(binwidth)
		binwidth = 1;
	end
	k = size(estimates.weights, 2);
	for f = fieldnames(estimates.params)'
		temp = repmat(estimates.params.(f{:})(1:k)', binwidth,1);
		simParams.(f{:})(1:k * binwidth) = temp(:);
		simParams.(f{:})(k * binwidth:end) = simParams.(f{:})(k*binwidth);
		figure(10)
		subplot(length(fieldnames(estimates.params)), 1, i)
		plot(simParams.(f{:}));
		i = i + 1;
	end
end

simParams.I = simParams.I + noise;

%% Run Sim
sim = zeros(numel(s0), TOTAL_TIME); % to hold simulation results

sim(:, 1) = s0;
for t = 2:TOTAL_TIME
	p = structfun(@(x) x(t - 1), simParams, 'Uni', 0);
	sim(:, t) = transitionFcn(sim(:, t-1), p);
end
t = (1:TOTAL_TIME) * dt;



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
