function [inds, beta] = anneal(weights, alpha, beta0)

if ~exist('alpha', 'var') || isempty(alpha)
	alpha = 0.5;  % particle survival rate
end
if ~exist('beta0', 'var') || isempty(beta0)
	beta0 = 1;
end
N = numel(weights);  % number of particles

% Pm = alpha^m * P0;  % covariance at layer m

D = @(beta) 1 / sum(weights .^ (beta * 2)) - alpha * N;
% D = @(beta) weights .^ (beta * 2) - 1 / (alpha * N);
% options = optimoptions('lsqnonlin', 'Display', 'none');
options = optimset('Display', 'off');

try
	[beta, ~, exitflag, ~] = fzero(D, beta0, options);
	if exitflag ~= 1
		beta = fzero(D, 1);
	end
catch ME
	try 
		beta = fzero(D, 1); 
	catch MEinner
		disp('debug')
	end
end

weights = (weights .^ beta) + 1e-6 * randn(size(weights));
weights = weights / sum(weights);

try
	inds = floor(interp1(cumsum(weights), 1:N, rand(1, N), 'linear', 0)) + 1;
catch ME
	disp('debug')
end
	
end





