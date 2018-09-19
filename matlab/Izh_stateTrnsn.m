function stateMat = Izh_stateTrnsn(stateMat, paramStruct, delta)
% Euler integration with noise

%% Set parameters

dF = model_dynamics([], stateMat, paramStruct, 'Izh');

stateMat = stateMat + delta * dF;  % integrate

spikeMask = stateMat(1, :) > 30;
stateMat(1, spikeMask) = paramStruct.c(spikeMask);
stateMat(2, spikeMask) = stateMat(2, spikeMask) + paramStruct.d(spikeMask);

if isfield(paramStruct, 'mNoise')
	stateMat(1, :) = stateMat(1, :) + ...
		paramStruct.mNoise * randn(1, size(stateMat, 2));  % add noise to V
end

