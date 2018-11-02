% function [] = pf_settings()

model =          'HH';			% select which model to use
SPIKETIMES =     'load';			% simulate ('newSim') or 'load' spike times
likelihood =     'spikes';		% 'voltage' or 'spikes'
spike_method =   'Vth';			% 'Vth' or 'diff'
N =              10e3;			% number of particles; Meng used 1e4, but start at 1e3 for speed
Nanneal =        2e3;			% number of particles for annealing
M =              0;				% annealing layers (try using only 200 particles with 10 layers)
K_MAX =          1000;			% Maximum number of time steps
alpha =          0.6;			% particle survival rate during annealing 
PARAMS =         default_params(model);  % Use default parameters for those not being estimated
% PARAMS =         getfield(load('alternate_params.mat'), 'p');  % Use non-default parameters 

[s0, boundsStruct] = HH_stateBounds(); % get initial conditions and parameter bounds

%% 
PLOT =           false;			% Plot particles while algorithm is running
% filename =       'pf.gif'; newgif = false;
