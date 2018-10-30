% function [] = pf_settings()

model = 'HH';	% select which model to use
SPIKETIMES = 'sim'; % simulate ('newSim') or 'load' spike times
likelihood = 'spikes';  % 'voltage' or 'spikes'
spike_method = 'Vth';  % 'Vth' or 'diff'
N = 2e3;  % number of particles; Meng used 1e4, but start at 1e3 for speed
Nanneal = 2e3;  % number of particles for annealing
M = 5;  % annealing layers (try using only 200 particles with 10 layers)
PLOT = false;  % Plot particles while algorithm is running
PLOT_RESULTS = true;  % Create summary plots when analysis is complete
K_MAX = Inf;  % Maximum number of time steps
alpha = 0.6;  % particle survival rate during annealing 
filename = 'pf.gif'; newgif = false;
% PARAMS = getfield(load('alternate_params.mat'), 'p');  % Use non-default parameters 
PARAMS = default_params(model);  % Use default parameters for those not being estimated
