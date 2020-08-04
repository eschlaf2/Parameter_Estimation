% function [] = pf_settings()

model =          'HH';			% select which model to use
SPIKETIMES =     'sim1';			% simulate ('newSim') or 'load' spike times
likelihood =     'voltage';		% 'voltage' or 'spikes'
spike_method =   'Vth';			% 'Vth' or 'diff'
N =              4e3;			% number of particles; Meng used 1e4, but start at 1e3 for speed
Nanneal =        2e3;			% number of particles for annealing
M =              5;				% annealing layers (try using only 200 particles with 10 layers)
K_MAX =          Inf;			% Maximum number of time steps
alpha =          0.6;			% particle survival rate during annealing 
% PARAMS =         default_params(model);  % Use default parameters for those not being estimated
PARAMS =         getfield(load('alternate_params.mat'), 'p');  % Use non-default parameters 

[s0, boundsStruct] = HH_stateBounds4(); % get initial conditions and parameter bounds

%% 
PLOT =           false;			% Plot particles while algorithm is running
% filename =       'pf.gif'; newgif = false;
