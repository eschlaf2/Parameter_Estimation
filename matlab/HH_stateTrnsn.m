function stateMat = HH_stateTrnsn(stateMat, paramStruct, delta, p)
% Euler integration with noise

%% Set parameters

dF = model_dynamics([], stateMat, paramStruct, 'HH', p);
% if ~isfield(paramStruct, 'mNoise')
% 	noise = default_params('HH');
% 	noise = noise.mNoise;
% else
% 	noise = paramStruct.mNoise;
% end

big_slopes = max(abs(dF)) > 5000;
if any(big_slopes)  % Use ode45 to integrate more slowly
	options = odeset('InitialStep', delta/5);
	for s = find(big_slopes)
		subset = structfun(@(x) x(s), paramStruct, 'uni', 0);
		[~, temp] = ode45(@(t, y) model_dynamics(t, y, subset, 'HH', p), ...
			single([0 delta]), stateMat(:, s), options);
		dF(:, s) = temp(end, :)' - stateMat(:, s);
	end
end

stateMat = stateMat + delta * dF;  % integrate
if isfield(paramStruct, 'mNoise')
	stateMat(1, :) = stateMat(1, :) + ...
		paramStruct.mNoise * randn(1, size(stateMat, 2));  % add noise to V
end

