function F = HH_dynamics(state, varargin)


V = state(1, :);
n = state(2, :);
h = state(3, :);
B = state(4, :);
gB = state(5, :);
EB = state(6, :);
VBth = state(7, :);
SB = state(8, :);
tauB = state(9, :);
I = state(10, :);

F = zeros(size(state));

%% Set parameters

mInf = 1./(1 + exp((-V-34.5)/10));	% HH options
nInf = 1./(1 + exp((-V-29.5)/10));
hInf = 1./(1 + exp((V+59.4)/10.7));
tauN = 0.25 + 4.35*exp(-abs(V+10)/10);
tauH = 0.15 + 1.15./(1 + exp((V+33.5)/15));
BInf = 1./(1 + exp(-(V-VBth)./SB));
tauB = tauB;
C = 0.9;
EK = -95;
ENa = 50;
EL = -70;
gNa = 100;
gK = 7;
gL = 0.25;

for i = 1:2:nargin-1
	if ~ischar(varargin{i}) || ~ischar(varargin{i + 1})
		error('All arguments should be given as strings.')
	end
	eval([varargin{i} '=' varargin{i+1} ';']);
end


%% Calculate changes
Vdot = (... % F1(V, n, h, B)
	I - gK * n.^4 .* (V - EK) ... % drive current minus Potassium current
	- gNa * mInf.^3 .* h .* (V - ENa) ... % Sodium current
	- gB .* B .* (V - EB) ... % mystery current
	- gL .* (V - EL)) / C; % Leak

ndot = (nInf - n) ./ tauN; % F2(n)

hdot = (hInf - h) ./ tauH; % F3(h)

Bdot = (BInf - B) ./ tauB; % F4(B)

F(1:4, :) = [Vdot; ndot; hdot; Bdot]; % state change

end
