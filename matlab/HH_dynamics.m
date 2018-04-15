function F = HH_dynamics(~, stateMat, paramStruct)
% state should be a matrix

%% Default parameter values
C = 0.9;
ENa = 50;
EL = -70;
gNa = 100;
gL = 0.25;
gB = 3.5;
EB = -74.8;
VBth = -2.2;
SB = 9.6;
tauB = 64;
I = 2;
EK = -95;
gK = 7;

%% Parse input
V = stateMat(1, :);
n = stateMat(2, :);
h = stateMat(3, :);
B = stateMat(4, :);

fn = fieldnames(paramStruct);
eval(sprintf('%s = paramStruct.%s;\n', fn{:}, fn{:}))

%% Set parameters

mInf = 1./(1 + exp((-V-34.5)/10));	% HH options
nInf = 1./(1 + exp((-V-29.5)/10));
hInf = 1./(1 + exp((V+59.4)/10.7));
tauN = 0.25 + 4.35*exp(-abs(V+10)/10);
tauH = 0.15 + 1.15./(1 + exp((V+33.5)/15));
BInf = 1./(1 + exp(-(V-VBth)./SB));
% tauB = tauB;

%% Calculate changes
Vdot = (... % F1(V, n, h, B)
	I - gK * n.^4 .* (V - EK) ... % drive current minus Potassium current
	- gNa * mInf.^3 .* h .* (V - ENa) ... % Sodium current
	- gB .* B .* (V - EB) ... % mystery current
	- gL .* (V - EL)) / C; ... % Leak

ndot = (nInf - n) ./ tauN; % F2(n)

hdot = (hInf - h) ./ tauH; % F3(h)

Bdot = (BInf - B) ./ tauB; % F4(B)

F = [Vdot; ndot; hdot; Bdot]; % state change

end
