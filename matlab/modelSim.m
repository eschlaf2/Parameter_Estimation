function [sim, spiketimes, simParams] = modelSim(model, thresh, estimates, binwidth, spike_method, p, varargin)

if ~exist('model', 'var') || isempty(model)
	model = 'HH';
end
if ~exist('thresh', 'var') || isempty(thresh)
	thresh = 30;
end
if ~exist('spike_method', 'var') || isempty(spike_method)
	spike_method = 'diff';
	thresh = 1;
end
if ~exist('p', 'var') || isempty(p)
	p = default_params(model);
end
dt = 1e-2;
TOTAL_TIME = 1e3 * 1/dt; % time steps to simulate (ms * fs)

for a = 1:2:length(varargin)
	switch lower(varargin{a})
		case 'dt'
			dt = varargin{a+1};
		case 'total_time'
			TOTAL_TIME = varargin{a+1} * 1 / dt;
		case 'total_steps'
			TOTAL_TIME = varargin{a + 1} * binwidth;
	end
	if isfield(p, varargin{a})
		p.(varargin{a}) = varargin{a+1};
	end
end

switch model
	case 'Izh'
		v =	-65;
		u = -10; 
		s0 = [v; u];
		
		noise = @(T) 4 * pinknoise(T);
		transitionFcn = @(state, pStruct) Izh_stateTrnsn(state, pStruct, dt, p);

		
	case 'HH'
		V =	-71;
		n = 0.0147; 
		h = 0.7497;
		B = 0.0326;
		s0 = [V; n; h; B];
		
		noise = @(T) 2 * pinknoise(T);
% 		noise = 0;
		transitionFcn = @(state, pStruct) HH_stateTrnsn(state, pStruct, dt, p);
		
end

% Create time series for each paramter
simParams = structfun(@(x) x * ones(1, TOTAL_TIME, 'single'), p, 'Uni', 0);
simParams.gB = linspace(0.5, 6, TOTAL_TIME);

% Use estimated parameters if given
i = 1;
if exist('estimates','var') && ~isempty(estimates)
	if ~exist('binwidth', 'var') || isempty(binwidth)
		binwidth = 1;
	end
	k = size(estimates.weights, 2);
% 	TOTAL_TIME = k;
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

simParams.I = simParams.I + noise(TOTAL_TIME);

%% Run Sim
sim = zeros(numel(s0), TOTAL_TIME); % to hold simulation results

sim(:, 1) = s0;
for t = 2:TOTAL_TIME
	p = structfun(@(x) x(t - 1), simParams, 'Uni', 0);
	sim(:, t) = transitionFcn(sim(:, t-1), p);
end
t = (1:TOTAL_TIME) * dt;

%% Get spiketimes
spiketimes = get_spiketimes(sim(1, :), thresh, spike_method);

%% Plot Results
if ~exist('PLOT_RESULTS', 'var') || PLOT_RESULTS
	plot_sim(sim, spiketimes, dt, thresh, 100-i);
end
