% sim

noiseStd = [1 0 0 0]; % Integration options
delta = 0.01;
Vth = 30;

V =	-70;
n = 0.0147; 
h = 0.7497;
B = 0.0326;

gB = 3.5;
EB = -74.8;
VBth = -2.2;
SB = 9.6;
tauB = 64;
I = -2.5;

s0 = [V; n; h; B];

params = [gB; EB; VBth; SB; tauB; I];


%% Run Sim
TOTAL_TIME = 1e3 * 1/delta; % time steps to simulate (ms * fs)
states = [s0; params]; % initial states
sim = zeros(numel(states), TOTAL_TIME); % to hold simulation results
noise = [noiseStd'; zeros(numel(params), 1)];

for t = 1:TOTAL_TIME
	dS = HH_dynamics(states); % dynamics
	states = states + dS * delta + noise .* randn(size(states)); % updated states
	states([false; states(2:4) > 1]) = 1;
	states([false; states(2:4) < 0]) = 0;
	sim(:, t) = states;
end

%% Plot Results
figure(99); fullwidth()
subplot(8,1,2:7);
tt = (1:TOTAL_TIME) * delta;
plot(tt, sim(1,:)); ylabel('Voltage [mV]'); xlabel('Time [ms]'); 
hold on; plot([0 TOTAL_TIME] * delta,Vth *  [1 1], '--', 'color', .5 * [1 1 1]); hold off
spikes = [false logical((sim(1,1:end-1) < Vth) .* sim(1,2:end) > Vth)];
hold on; plot(tt(spikes), sim(1,spikes), 'r*'); hold off;

subplot(8,1,1);
spiketimes = find(spikes);
plot([spiketimes(:) spiketimes(:)]' * delta, [zeros(size(spiketimes(:))) ones(size(spiketimes(:)))]', 'color', lines(1))
yticks([]); xticks([]);
