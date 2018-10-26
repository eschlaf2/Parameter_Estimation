%% Set parameters
if ~exist('outfile', 'var')
	clear
else
	rng(sum(double(outfile)));
end

model = 'HH';	% select which model to use
SPIKETIMES = 'sim'; % simulate ('newSim') or 'load' spike times
likelihood = 'spikes';  % voltage or spikes
N = 2e3;  % number of particles; Meng used 1e4, but start at 1e3 for speed
Nanneal = 2e3;  % number of particles for annealing
M = 5;  % annealing layers (try using only 200 particles with 10 layers)
PLOT = false;  % Plot particles while algorithm is running
PLOT_RESULTS = false;  % Create summary plots when analysis is complete
% STATE_BOUNDS = [-200, 200, -20, 20];
K_MAX = Inf;  % Maximum number of time steps
alpha = 0.6;  % particle survival rate during annealing 
filename = 'pf.gif'; newgif = false;
% injectedCurrent = pinknoise(K_MAX);

%% Load firing times
 
switch SPIKETIMES
	case 'load'
		load firings.mat	% sorted and curated firings from binary_small.dat (preictal); 
							% sorting by MountainSort
		spiketimes = firings(2,firings(3,:) == 95);	% spike times from unit 95
		fs = 3e4;	% sampling frequency [Hz]
	case 'sim'
		load sim_2e3_noise2_gbVariable.mat
% 		load sim0
		fs = 1e5; 
	case 'newSim'
		switch model
			case 'HH'
				Vth = 30;
		% 		HHSim;
				fs = 1e5;
			case 'Izh'
				Vth = 20;
				fs = 1e5;
		end
		[simV, spiketimes, simParams] = modelSim(model, Vth);
		meanSpike = mean(cell2mat(arrayfun(@(i) simV(1, i - 50:i+100), spiketimes, 'uni', 0)'));
end

%% Select model
switch model
	case 'why'	% place holder
		why
	case 'Izh'
		[s0, boundsStruct] = Izh_stateBounds();
		
		stateNames = {'v', 'u'};
		
		dt = 0.01;	% integration step [ms]
		
		W = 3;	% window (ms)
		
		Vth = 20; % Voltage threshold [mV]
		
		delta = 1; % binwidth [ms]
		
		transitionFcn = @(states, particles) Izh_stateTrnsn(states, particles, dt);
		switch likelihood
			case 'spikes'
				likelihoodFcn = @(window, obsn) ...
					likelihoodFcnMeng(window, obsn, W, Vth, delta); 
			case 'voltage'
				likelihoodFcn = @(window, obsn) ...
					likelihood_voltage(window, obsn, sigma);
			case '2011'
				likelihoodFcn = @(window, obsn) ...
					likelihoodFcnMeng2011(window, obsn, ((dt:dt:W) - W/2), Vth);% 			likelihood_voltage(window, obsn, sigma);
		end
			
		resamplingFcn = @resamplingMeng;
		
		
	case 'HH'	% Hodgkin-Huxley
		[s0, boundsStruct] = HH_stateBounds(); % get initial conditions and parameter bounds
		
		stateNames = {'V', 'n', 'h', 'B'};
		
% 		procNoise = 0.02; % covariance of process noise as proportion of range
		dt = 0.01;	% integration step [ms]
		
		W = 3;	% window (ms)
		Vth = 30; % Voltage threshold [mV]
		
		delta = .3; % binwidth [ms]
		sigma = 1;
		
		transitionFcn = @(states, particles) HH_stateTrnsn(states, particles, dt);
		switch likelihood
			case 'spikes'
				likelihoodFcn = @(window, obsn) ...
					likelihoodFcnMeng(window, obsn, W, Vth, delta); 
			case 'voltage'
				likelihoodFcn = @(window, obsn) ...
					likelihood_voltage(window, obsn, sigma);
			case '2011'
				likelihoodFcn = @(window, obsn) ...
					likelihoodFcnMeng2011(window, obsn, ((dt:dt:W) - W/2), Vth);% 			likelihood_voltage(window, obsn, sigma);
		end

		sigma = 1;
			
		resamplingFcn = @resamplingMeng;
		
		
end

%% Bin spike times
binwidth = fs * delta * 1e-3; % samples per bin
binedges = 0:binwidth:(max(spiketimes) + binwidth);
% tSpan = (1: length(binedges)-1) * delta * 1e-3;	% time [s]
obsn = histcounts(spiketimes, binedges);

% obsn = simV(1, 1:binwidth:end);
tSpan = (1:length(obsn) - W) * delta * 1e-3;
% obsnV = zeros(1, ceil(spiketimes(end)/binwidth) * binwidth);
% obsnV(spiketimes) = 1;
% obsnV = conv(obsnV, meanSpike - min(meanSpike), 'same') + min(meanSpike);
% obsnV = max(reshape(obsnV, binwidth, []), [], 1);

%% Set convenience variables
NUM_STATES = size(s0, 1);
NUM_PARAMS = size(fieldnames(boundsStruct), 1);
NUM_ALL = NUM_STATES + NUM_PARAMS + 1;
K = length(obsn) - W;
fn = fieldnames(boundsStruct)';


%% Initialize PF

% Initialize parameter particles
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
particles.pNoise = structfun(@(x) x(4) .* range(x(1:2)), boundsStruct); % noise to be injected into the filter

% paramDistX = cell2mat(arrayfun(@(i) linspace(boundsStruct(i, 1), boundsStruct(i, 2), 1e3), pEst, ...
% 	'UniformOutput', false));
paramDist = structfun(@(x) zeros(1e3, K, 'single'), boundsStruct, 'uni', 0); % for holding (interpolated) posteriors

paramDistX = structfun(@(x) linspace(x(1), x(2), 1e3), boundsStruct, ...
	'UniformOutput', false);
triggered = zeros(1, K);
beta = ones(M, 1);  % Initial exponents for annealing
P0 = diag(particles.pNoise);  % Covariance matrix - see APF paper for ideas on how to improve this
P = @(m, i) particles.pNoise(i) * alpha ^ m;
ESS = 1;

%% Run PF
for k = 1:min(K, K_MAX)		% for each observation
	
	prediction = particles.states;
	for i = 1:binwidth  % for each integration step within a bin (advance to t + 1)
		prediction = transitionFcn(prediction, particles.params);	% ... integrate states
% 		try
% 			prediction(prediction < STATE_BOUNDS(1)) = STATE_BOUNDS(1);
% 			prediction(prediction > STATE_BOUNDS(2)) = STATE_BOUNDS(2);
% 		catch ME
% 		end
	end
	particles.states = prediction;
	
	params = particles.params;
	window = updateWindow(prediction, max(W*binwidth, 1), @(s) transitionFcn(s, params));
	probability = likelihoodFcn(window, sum(obsn(k:k+W-1))); % ... calculate likelihood
% 	probability = likelihoodFcn(window(:, :, 1:binwidth:end), obsn(k:k+W-1)); % ... calculate likelihood
	
% 	if M < 1  % if no annealing, use the resampling function
		trig = obsn(k);
		[particles, inds] = resamplingFcn(particles, probability, trig); % ... resample particles	
		triggered(k) = trig;
		wts = particles.weights;
% 	else  % else assign weights based on single probability calculation
% 		wts = probability + 1e-6;
% 		wts = wts / sum(wts);
% 	end
% 	
	particles.params = keep_in_bounds(particles.params, boundsStruct);
	for m = M:-1:1  % anneal
		
		if m > 1
			Nm = Nanneal;
		else
			Nm = N;
		end
		
		[inds, beta(m)] = anneal(wts, alpha, beta(m), Nm);  % get betas and sample indices

		for ii = 1:NUM_PARAMS  % add noise
			params.(fn{ii}) = params.(fn{ii})(inds) + ...
				P(m, ii) * randn(1, Nm);
		end
		params = keep_in_bounds(params, boundsStruct);
		
		prediction = prediction(:, inds); % resample
		% Update window of V surrounding t_k
		window = updateWindow(prediction, max(W*binwidth, 1), @(s) transitionFcn(s, params));
		probability = likelihoodFcn(window, sum(obsn(k:k+W-1))); % ...
% 		calculate likelihood base
% 		probability = likelihoodFcn(window(:, :, 1:binwidth:end), obsn(k:k+W-1)); % ... calculate likelihood based on voltage

		
		wts = particles.weights(inds) .* probability + 1e-6;
		wts = wts / sum(wts);
		particles.params = params;
		particles.states = prediction;
% 		particles.weights = particles.weights(inds) / sum(particles.weights(inds));
		particles.weights = wts;
	end
	
	weighted_mean = @(x) sum(x .* particles.weights, 2);
	estimates.states(:, k) = weighted_mean(particles.states);
	estimates.weights(k) = median(particles.weights); 
	ESS = 1 / sum(particles.weights.^2) / N; 
	estimates.ESS(k) = ESS;
	
% 	if M > 0  % if annealing, reset particle weights
% 		particles.weights = 1/N * ones(1, N);
% 	end
	
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

	if ~mod(k, 10)
		disp(['k = ' num2str(k)])
	end
	
	if PLOT && ~mod(k, 1)
		h = figure(999);
		x = fn{1}; y = fn{2};
		inds = randsample(N, 100);
		scatter(particles.params.(x)(inds), particles.params.(y)(inds), 16, particles.weights(inds), 'filled');
		hold on; plot(simParams.(x)(k * binwidth), simParams.(y)(k*binwidth), 'r*'); hold off;
		hold on; plot(estimates.params.(x)(k), estimates.params.(y)(k), 'b*'); hold off;
		xlim((boundsStruct.(x)(1:2) - mean(boundsStruct.(x)(1:2))) * 2*boundsStruct.(x)(3) + mean(boundsStruct.(x)(1:2))); 
		ylim((boundsStruct.(y)(1:2) - mean(boundsStruct.(y)(1:2))) * 2*boundsStruct.(y)(3) + mean(boundsStruct.(y)(1:2))); 
		colormap('cool'); caxis([0 2/N]); colorbar
		title(sprintf('%1.2f: %d', tSpan(k) * 1e3, obsn(k)))
		xlabel(x); ylabel(y);
		drawnow;
% 		pause(1e-6)
		
		if 1
			frame = getframe(h); 
			im = frame2im(frame); 
			[imind,cm] = rgb2ind(im,256); 
			% Write to the GIF File 
			if newgif
			  imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'delaytime', 0); 
			  newgif = false;
			else 
			  imwrite(imind,cm,filename,'gif','WriteMode','append', 'delaytime', 0); 
			end 
		end
	end
end

for f = fn
	f = f{:};
	estimates.params.(f)(k:end) = estimates.params.(f)(k);
end

beep

%% Plot results
if PLOT_RESULTS
    colors = lines(max(NUM_PARAMS, 7));
%     k = k-1;

    figure(3); fullwidth(0)
    stem(tSpan(obsn(1:k) > 0), obsn(obsn(1:k) > 0)', 'k', 'linewidth', 2); hold on;
    plot(tSpan(1:k), estimates.states(1, 1:k), 'Color', colors(2,:));
	if ~strcmp(SPIKETIMES, 'load')
		plot(tSpan(1:k), obsn(1:k), 'Color', colors(1,:))
	end
    plot(tSpan(1:k), ones(1, k), '--', 'Color', .5*[1 1 1]); hold off;
    ylim(1.1*(get(gca,'YLim') - mean(get(gca, 'Ylim'))) + mean(get(gca, 'YLim')));
    legend('Spikes','Estimated Voltage', 'True Voltage', 'location', 'best');
    xlabel('Time [s]'); ylabel('Voltage [mV]');
    title('Voltage')

    figure(4); fullwidth()
    plot((tSpan(1:k) .* ones(3, k))', estimates.states(2:end, 1:k)', 'linewidth', 2); 
    legend(stateNames(2:end))
    hold on; set(gca, 'ColorOrderIndex', 1)
	try
		plot((tSpan(1:k) .* ones(3, k))', simV(2:end, 1:binwidth:k*binwidth)', ':', 'linewidth', 2); hold off;
	catch ME
	end
    xlabel('Time [s]');
    title('Hidden States')

    figure(5); clf; fullwidth(numel(fn) > 3)
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
		if ~strcmp(SPIKETIMES, 'load')
			plot(tSpan(1:k), simParams.(f)(1:binwidth:k*binwidth)', '--', ...
				'displayname', sprintf('%s True', f), ...
				'Color', colors(i, :)); 
		end
        hold off;
        legend();
        title(['Estimated ', f])
        axis('tight')
    end
    xlabel('Time [s]');


    figure(6); fullwidth()
	wts = estimates.weights(1:k);
	subplot(2,1,1)
    plot(wts); title('Mean Weight');
	hold on;
	xx = 1:k;
	plot(xx(triggered(1:k) > 0), wts(triggered(1:k) > 0), 'r*'); hold off;
	subplot(2,1,2)
	plot(estimates.ESS(1:k) * 100); title('ESS%')

	PLOT7 = false;
	if PLOT7
		figure(7); clf; fullwidth(numel(fn) > 3)
		for i = 1:numel(fn)
			f = fn{i};
			subplot(numel(fn), 1, i)
			contourf(tSpan(1:k-1), paramDistX.(f), squeeze(paramDist.(f)(:, 1:k-1)), 'linestyle', 'none')
			caxis([0 2/N]); % colormap('gray')
			colorbar()
			if ~strcmp(SPIKETIMES, 'load')
				hold on; 
				plot(tSpan(1:k)', simParams.(f)(1:binwidth:k*binwidth)', 'r--', 'linewidth', 3); 
				hold off; 
			end
			title([f ' Estimate'])
		end
	end
    xlabel('Time [s]');
end
[simEst, stEst, ~] = modelSim(model, Vth, estimates, binwidth, 'total_steps', size(simV, 2));  % simulate based on estimated parameters

if exist('outfile', 'var')
%    save(outfile, 'sim', 'estimates')
   save(outfile);
   if PLOT_RESULTS
	   print(5, [outfile '_5'], '-dpng')
	   if PLOT7
		   print(7, [outfile '_7'], '-dpng')
	   end
	   if strcmp(SPIKETIMES, 'newSim')
		   print(99, [outfile '_sim'], '-dpng')
	   end
   end
end

disp('success')