% function [] = forward_density()

model = 'HH';
Vth = 30;
delta = 1e-2;  % integration step (ms)

[sim, spiketimes, simParams] = modelSim(model, Vth);

s = ceil(spiketimes(1) * delta);  % next spike time (ms)
hiddenStates = sim(2:end, 1:(dt/delta):(s/delta));
dt = .1;  % time step of probability space (ms)
[s0, ~] = HH_stateBounds(); % get initial conditions and parameter bounds
paramStruct = default_params(model);
mNoise = 1;  % measurement noise (for now assume that voltage lands in a gaussian around deterministic value)
dV = 1;
Vrange = (-100:dV:Vth);  % Set the space

pF = zeros(length(Vrange), s/dt);
gaussian = @(x, mu, sigma) 1 ./ (sqrt(2 * pi * sigma^2)) .* exp(-(x - mu).^2 / (2 * sigma^2));

vProj = zeros(size(pF));  % save for use with backward density
tempState = s0;
for dummy = 1:(dt / delta)
	tempState = HH_stateTrnsn(tempState, paramStruct, delta);
end
% vProj = repmat(s0, 1, length(Vrange));
vProj(:, 1) = tempState(1);

pF(:, 1) = arrayfun(@(v) gaussian(v, vProj(1, 1), mNoise), Vrange);
pF(:, 1) = pF(:, 1) / sum(pF(:, 1));

% vProj(1, :) = Vrange;
path = zeros(1, s/dt);
path(1) = Vrange(find(rand(1) < cumsum(pF(:, 1)), 1) - 1);

%% Compute forward density
for t = 2:(s/dt)
	
	% Compute projections from all V1s
	tempState = [Vrange; repmat(hiddenStates(:, t), 1, length(Vrange))];
	for Vi = 1:length(Vrange)
		for dummy = 1:(dt / delta)
			tempState(:, Vi) = HH_stateTrnsn(tempState(:, Vi), paramStruct, delta);
		end
	end
	vProj(:, t) = tempState(1, :);
% 	vProj(1, vProj(1, :) > Vth) = Vth;
	
	for Vi = 1:length(Vrange)  
		
		% Does not use simulated hidden states
		if 0  % wrong... needs to be updated
			for dummy = 1:(dtT / delta)  % advance one full step (dtT)
				vProj(:, i) = HH_stateTrnsn(vProj(:, i), paramStruct, delta);
			end

			% Compute forward density
			p(i, t) = sum(arrayfun(@(v) ...
				gaussian(vProj(1, i), v, mNoise), ...
				Vrange)' .* p(:, t - 1));

			% Reset states using interpolation for next update
			% ... method is 'pchip' to avoid NaN's in states
			for state = 2:4
				vProj(state, :) = interp1(vProj(1,:), vProj(state,:), Vrange, 'pchip');
			end
			vProj(1, :) = Vrange;

		end
		
		% Uses simulated hidden states
		if 1
			% Compute forward density
			pF(Vi, t) = sum(arrayfun(@(v) ...
				gaussian(Vrange(Vi), v, mNoise), vProj(:, t)) .* pF(:, t - 1));
		end
		
		if path(t-1) == Vrange(Vi)
			tempdist = arrayfun(@(v) gaussian(v, vProj(Vi, t), mNoise), Vrange) + 1e-6;
			tempdist(isnan(tempdist)) = 1e-6;
			tempdist = tempdist / sum(tempdist);
			path(t) = Vrange(find(rand(1) <= cumsum(tempdist), 1));
		end
		
	end
	pF(:, t) = pF(:, t) / sum(pF(:, t));  % normalize p
	
	% Show progress
	if mod(t, 10^(log10(s / dt) >= 1)) == 0
		disp(['t = ' num2str(t) '/' num2str(s/dt)])
	end
end

%% Plot results
figure(84411201); imagesc((1:t), Vrange, (pF(:, 1:t)))
axis('xy'); colorbar
hold on; plot(sim(1, 1:dt/delta:t/(dt/delta)), 'r', 'linewidth', 2); hold off
disp('success')
beep

%% Sample path
pp = cumsum(pF);
figure(1); plot(Vrange, pp(:, 1))
pathRnd = rand(1, t);

figure(2); plot(path)

%% backward density

pB = zeros(size(pF));
pB(end, end) = 1;

for k = s/dt-1:-1:1
	for vi = 1:length(Vrange)
		tempInt = arrayfun(@(i) gaussian(Vrange(i), vProj(vi, k), mNoise) * ...
			pF(vi, k) / pF(i, k+1) * pB(i, k+1), 1:length(Vrange));
		pB(vi, k) = sum(tempInt);
	end
	pB(:, k) = pB(:, k) + 1e-6;
	pB(:, k) = pB(:, k) / sum(pB(:, k));
	
end
		
%% Plot results
figure(97411201); imagesc((1:s/dt), Vrange, (pB(:, 1:s/dt)))
% figure(97411202); imagesc((1:s/dt), Vrange, pB(:, 1:s/dt) ./ pF(:, 1:s/dt))
axis('xy'); colorbar
% hold on; plot(sim(1, 1:dt/delta:t/(dt/delta)), 'r', 'linewidth', 2); hold off


%% Path simulation

p = zeros(1, s/dt);
state = s0;
p(1) = state(1);

for t = 2:s/dt
	for dummy = 1:(dt / delta)
		state = HH_stateTrnsn(state, paramStruct, delta);
	end
	tempDist = gaussian(Vrange, state(1), mNoise)' + 1e-6;
	tempDist(isnan(tempDist)) = 1e-6;
	tempDist = tempDist / sum(tempDist); 
	tempDist = tempDist ./ pF(:, t) .* pB(:, t) + 1e-6;
	figure(2); hold off; plot(Vrange, tempDist / sum(tempDist)); pause(.1)
	tempDist = cumsum(tempDist / sum(tempDist));
% 	[~, uniqInd] = unique(tempDist);
	tempRnd = rand(1);
	p(t) = Vrange(find(tempDist >= tempRnd, 1));
% 	p(t) = interp1(tempDist(uniqInd), Vrange(uniqInd), tempRnd);
% 	p(t) = interp1(tempDist, Vrange, tempRnd);
	state(1) = p(t);
end
figure(1); plot(p)
beep