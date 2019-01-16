function [spiketimes] = get_spiketimes(V, thresh, method)

%%

if ~exist('method', 'var') || isempty(method)
	method = 'diff';  % 'Vth' or 'diff'
end

switch method
	case 'Vth'		
% 		spikes = [false logical((V(1:end-1) < thresh) .* ...
% 			V(2:end) > thresh)];
		spikes = diff((V > thresh)) > 0;
	case 'diff'
		dV = diff(V);
		spikes = diff(dV > thresh) > 0;
end

spiketimes = find(spikes); % spike times in samples
