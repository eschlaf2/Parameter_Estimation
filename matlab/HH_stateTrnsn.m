function stateMat = HH_stateTrnsn(stateMat, paramStruct, delta)
% Euler integration with noise

%% Set parameters

dF = HH_dynamics([], stateMat, paramStruct);
if ~isfield(paramStruct, 'mNoise')
	noise = 0.1;
else
	noise = paramStruct.mNoise;
end

big_slopes = max(abs(dF)) > 5000;
if any(big_slopes)
	options = odeset('InitialStep', delta/5);
	for p = find(big_slopes)
		subset = structfun(@(x) x(p), paramStruct, 'uni', 0);
		[~, temp] = ode45(@(t, y) HH_dynamics(t, y, subset), ...
			[0 delta], stateMat(:, p), options);
		dF(:, p) = temp(end, :)' - stateMat(:, p);
	end
end

stateMat = stateMat + delta * dF;  % integrate
stateMat(1, :) = stateMat(1, :) + ...  % add noise to V
	noise .* randn(1, size(stateMat, 2));

