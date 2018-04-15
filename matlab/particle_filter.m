%% Set parameters

model = 'HH';	% select which model to use
SPIKETIMES = 'sim'; % simulate ('sim') or 'load' spike times
N = 1e3; % number of particles; Meng used 1e4, but start at 1e3 for speed
PLOT = false;
PLOT_RESULTS = true;
K_MAX = 1e2;

%% Load firing times

switch SPIKETIMES
	case 'load'
		load firings.mat	% sorted and curated firings from binary_small.dat (preictal); 
							% sorting by MountainSort
		spiketimes = firings(2,firings(3,:) == 95);	% spike times from unit 95
		fs = 3e4;	% sampling frequency [Hz]
	case 'sim'
		load sim.mat
		fs = 1e5;
	case 'new_sim'
		HHSim;
		fs = 1e5;
end

%% Select model
switch model
	case 'why'	% easter egg
		why
	case 'HH'	% Hodgkin-Huxley
		[s0, boundsStruct] = HH_stateBounds(); % get initial conditions and parameter bounds
		
		stateNames = {'V', 'n', 'h', 'B'};
		
		procNoise = 0.01; % std of process noise as proportion of range
		dt = 0.01;	% integration step [ms]
		
		W = 3;	% window (ms)
		Vth = 30; % Voltage threshold [mV]
		
		delta = 1; % binwidth [ms]

		
		transitionFcn = @(states, particles) HH_stateTrnsn(states, particles, dt);
		likelihoodFcn = @(window, obsn) ...
			likelihoodFcnMeng(window, obsn, W, Vth, delta); 
% 			likelihoodFcnMeng2011(window, obsn, ((dt:dt:W) - W/2), Vth);
			
		resamplingFcn = @resamplingMeng;
		
		
end

%% Bin spike times
binwidth = fs * delta * 1e-3; % samples per bin
binedges = 0:binwidth:(max(spiketimes) + binwidth);
tSpan = (1: length(binedges)-1) * delta * 1e-3;	% time [s]
obsn = histcounts(spiketimes, binedges);
% obsn = mean(reshape(sim(1,:)', binwidth, []));

%% Set convenience variables
NUM_STATES = size(s0, 1);
NUM_PARAMS = size(boundsStruct, 1);
NUM_ALL = NUM_STATES + NUM_PARAMS + 1;
K = length(obsn) - W;
fn = fieldnames(boundsStruct)';


%% Initialize PF

% Initialize parameter particles
% posterior = cell2mat(arrayfun(@(a,b) unifrnd(a, b, 1, N), ...
% 	paramBounds(:, 1), paramBounds(:, 2), ...
% 	'uniformoutput', false));

particles.params = structfun(@(x) unifrnd(x(1), x(2), 1, N), boundsStruct, ...
	'uniformoutput', false);

% Initialize estimated states, weights, first prior
estimates.states = zeros(NUM_STATES, K); % for holding estimates
estimates.params = structfun(@(x) zeros(K, 1), boundsStruct, 'Uni', 0);
estimates.ci = structfun(@(x) zeros(K, 2), boundsStruct, 'Uni', 0); % for holding 95% CI
estimates.states(:, 1) = s0; % set initial conditions
estimates.weights = zeros(1, K);
particles.weights = 1/N * ones(1, N); % initialize all weights to 1/N

particles.states = s0 .* ones(NUM_STATES, N); % states associated with each particle
particles.pNoise = structfun(@(x) procNoise .* range(x(1:2)), boundsStruct); % noise to be injected into the filter

% paramDistX = cell2mat(arrayfun(@(i) linspace(boundsStruct(i, 1), boundsStruct(i, 2), 1e3), pEst, ...
% 	'UniformOutput', false));
paramDist = structfun(@(x) zeros(1e3, K, 'single'), boundsStruct, 'uni', 0); % for holding (interpolated) posteriors

paramDistX = structfun(@(x) linspace(x(1), x(2), 1e3), boundsStruct, ...
	'UniformOutput', false);


%% Run PF
for k = 1:min(K, K_MAX)		% for each observation
	
	prediction = particles.states;
	
	for i = 1:binwidth % for each integration step within a bin
		prediction = transitionFcn(prediction, particles.params);	% ... integrate states
	end
	
	% Update window of V surrounding t_k
	window = updateWindow(prediction, W*binwidth, @(s) transitionFcn(s, particles.params));
	
	probability = likelihoodFcn(window, sum(obsn(k:k+W-1))); % ... calculate likelihood
% 	[posterior, inds] = resamplingFcn(particles, probability, obsn(k)); % ... resample particles
	[particles, inds] = resamplingFcn(particles, probability, 0); % ... resample particles
	
	
	particles.params = keep_in_bounds(particles.params, boundsStruct);
	weighted_mean = @(x) sum(x .* particles.weights, 2);
	estimates.states(:, k) = weighted_mean(particles.states);
	estimates.weights(k) = mean(particles.weights); 
	temp = structfun(@sort, particles.params, 'uni', 0);
	
	for f = fn
		f = f{:};
		
		paramDist.(f)(:, k) = interp1(double(particles.params.(f)), particles.weights, ...
			paramDistX.(f), 'linear', 0);  % ... save the distribution
		
		% Get next estimates using weighted mean
		estimates.params.(f)(k) = weighted_mean(particles.params.(f));
		
		% Store the CIs
		
		estimates.ci.(f)(k, :) = temp.(f)(:, [floor(N*0.025) ceil(N*0.975)]);
		
	end
% 	for p = 1:numel(pEst)
% 		paramDist(p, :, k) = ...
% 			interp1(posterior(NUM_STATES+pEst(p), :), posterior(end,:), paramDistX(p, :), 'linear', 0); % ... save distribution
% 	end	
	
	if k == 10
		disp('look around')
	end
	
	if PLOT && ~mod(k, 5)
		figure(999)
		x = pEst(1); y = pEst(2);
		inds = randsample(N, 100);
		scatter(particles(NUM_STATES + x, inds), particles(NUM_STATES + y, inds), 16, particles(end, inds), 'filled');
		hold on; plot(sim(x + NUM_STATES, k * binwidth), sim(y + NUM_STATES, k*binwidth), 'r*'); hold off;
		hold on; plot(estimates(x + NUM_STATES, k), estimates(y + NUM_STATES, k), 'b*'); hold off;
		xlim(boundsStruct(x,1:2)); ylim(boundsStruct(y,1:2));
		colormap('cool'); caxis([0 2/N]); colorbar
		title(sprintf('%1.2f: %d', tSpan(k) * 1e3, obsn(k)))
		xlabel(paramNames{pEst(1)}); ylabel(paramNames{pEst(2)});
		drawnow;
		pause(1e-6)
	end
end

%% Plot results
if PLOT_RESULTS
    colors = lines(7);
    k = k-1;

    figure(3); fullwidth()
    stem(tSpan(obsn(1:k) > 0), obsn(obsn(1:k) > 0)', 'k', 'linewidth', 2); hold on;
    plot(tSpan(1:k), estimates.states(1, 1:k)/Vth, 'Color', colors(2,:));
    plot(tSpan(1:k), sim(1, 1:binwidth:k*binwidth)/Vth, 'Color', colors(1,:))
    plot(tSpan(1:k), ones(1, k), '--', 'Color', .5*[1 1 1]); hold off;
    ylim(1.1*(get(gca,'YLim') - mean(get(gca, 'Ylim'))) + mean(get(gca, 'YLim')));
    legend('Spikes','Estimated Voltage', 'True Voltage', 'location', 'best');
    xlabel('Time [s]'); ylabel('Voltage [mV]');
    title('Voltage')

    figure(4); fullwidth()
    plot((tSpan(1:k) .* ones(3, k))', estimates.states(2:4, 1:k)', 'linewidth', 2); 
    legend(stateNames(2:4))
    hold on; set(gca, 'ColorOrderIndex', 1)
    plot((tSpan(1:k) .* ones(3, k))', sim(2:4, 1:binwidth:k*binwidth)', ':', 'linewidth', 2); hold off;
    xlabel('Time [s]');
    title('Hidden States')

    figure(5); clf; fullwidth(1)
    for i = 1:numel(fn)
		f = fn{i};
        subplot(numel(fn), 1, i)
        cla;
        shapeX = tSpan([1:k, k:-1:1]);
        shapeY = [estimates.ci.(f)(1:k, 1)', estimates.ci.(f)(k:-1:1, 2)'];
        patch(shapeX, shapeY, 'b', 'FaceColor', colors(i, :), 'facealpha', .3, ...
            'edgecolor', 'none', 'displayname', '95CI');
        hold on;
        plot(tSpan(1:k), estimates.params.(f)(1:k), 'linewidth', 2, ...
            'displayname', sprintf('%s Est', f), ...
            'Color', colors(i, :)); 
        plot(tSpan(1:k), simParams.(f)(1:binwidth:k*binwidth)', '--', ...
            'displayname', sprintf('%s True', f), ...
            'Color', colors(i, :)); 
        hold off;
        legend();
        title(['Estimated ', f])
        axis('tight')
    end
    xlabel('Time [s]');


    figure(6); fullwidth()
%     plot(estimates.weights(1:k)); title('Mean Likelihood');

    figure(7); fullwidth(numel(fn) > 1)
    for i = 1:numel(fn)
		f = fn{i};
        subplot(numel(fn), 1, i)
        contourf(tSpan(1:k-1), paramDistX.(f), squeeze(paramDist.(f)(:, 1:k-1)), 'linestyle', 'none')
        caxis([0 2/N]); % colormap('gray')
        colorbar()
        hold on; 
        plot(tSpan(1:k)', simParams.(f)(1:binwidth:k*binwidth)', 'r--', 'linewidth', 3); 
        hold off; 
        title([f ' Estimate'])
    end
    xlabel('Time [s]');
end

if exist('outfile', 'var')
   save(outfile, 'sim', 'estimates')
   print(5, [outfile '_5'], '-dpng')
   print(7, [outfile '_7'], '-dpng')
   print(99, [outfile '_sim'], '-dpng')
end
%% Supplementary functions

function window = updateWindow(prediction, W, transitionFcn)

	window = zeros([size(prediction) W]);

	window(:, :, 1) = prediction;
	for i = 2:W
		window(:, :, i) = transitionFcn(window(:, :, i - 1));
	end
		
end

function bounded = keep_in_bounds(data, bounds)
% hit the wall

bounded = data;

fn = fieldnames(bounds)';
for i = 1:numel(fn)
	f = fn{i};
	b = bounds.(f);
	if ~b(3), continue; end  % if parameter is not bounded continue
	b(1:2) = (b(1:2) - mean(b(1:2))) * b(3) + mean(b(1:2));  % adjust range
	d = data.(f);

	indsHi = d > b(2);
	indsLo = d < b(1);
	% High values move inside range
	d(indsHi) = b(2) - abs(randn(1, sum(indsHi)) * .02 * diff(b(1:2)));  
	% ... and Low values
	d(indsLo) = b(1) + abs(randn(1, sum(indsLo)) * .02 * diff(b(1:2)));
	bounded.(f) = d;
end

end
