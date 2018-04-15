function stateMat = HH_stateTrnsn(stateMat, paramStruct)
% Euler integration with noise

%% Set parameters

dF = HH_dynamics([], stateMat, paramStruct);

big_slopes = max(abs(dF)) > 5000;
if any(big_slopes)
	options = odeset('InitialStep', delta/5);
	for p = find(big_slopes)
		[~, temp] = ode45(@(t, y) HH_dynamics(t, y, paramStruct), ...
			[0 delta], stateMat(:, p), options);
		dF(:, p) = temp(end, :)' - stateMat(:, p);
	end
end

stateMat = stateMat + delta * dF;  % integrate
stateMat(1, :) = stateMat(1, :) + ...  % add noise to V
	paramStruct.mNoise .* randn(1, size(stateMat, 2));

