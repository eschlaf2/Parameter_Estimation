% sim

% noiseStd = [1 0 0 0]; % ... Measurment noise
delta = 0.01; % integration step [ms]
Vth = 30; % count spikes when voltage goes above Vth
TOTAL_TIME = 1e3 * 1/delta; % time steps to simulate (ms * fs)

V =	-71;
n = 0.0147; 
h = 0.7497;
B = 0.0326;

gB = 3.5 * ones(1, TOTAL_TIME);
EB = linspace(-90, -60, TOTAL_TIME);
VBth = -2.2 * ones(1, TOTAL_TIME);
SB = 9.6 * ones(1, TOTAL_TIME);
tauB = 64 * ones(1, TOTAL_TIME);
I = 2 * ones(1, TOTAL_TIME);
mNoise = 0.1 * ones(1, TOTAL_TIME);

s0 = [V; n; h; B];

params = [gB; EB; VBth; SB; tauB; I; mNoise];

%% Run Sim
sim = [zeros(numel(s0), TOTAL_TIME); params]; % to hold simulation results

sim(:, 1) = [s0; params(:, 1)];
for t = 2:TOTAL_TIME
	sim(:, t) = HH_stateTrnsn(sim(:, t-1), [], delta, []);
	sim(numel(s0) + 1:end, t) = params(:, t);
end
t = (1:TOTAL_TIME) * delta;


%% Plot Results
if ~exist('PLOT_RESULTS', 'var') || PLOT_RESULTS
    figure(99); fullwidth() % Voltage
    ax = subplot(8,1,2:7);
    plot(t, sim(1,:)); ylabel('Voltage [mV]'); xlabel('Time [ms]'); 
    hold on; plot([0 t(end)],Vth *  [1 1], '--', 'color', .5 * [1 1 1]); hold off
    spikes = [false logical((sim(1,1:end-1) < Vth) .* sim(1,2:end) > Vth)];
    hold on; plot(t(spikes), sim(1,spikes), 'r*'); hold off;

    subplot(8,1,1); % Spike raster
    spiketimes = find(spikes);
    spikeT = t(spikes);
    plot([spikeT(:) spikeT(:)]', [zeros(size(spikeT(:))) ones(size(spikeT(:)))]', 'color', lines(1))
    xlim(get(ax, 'xlim'))
    yticks([]); xticks([]);
end
