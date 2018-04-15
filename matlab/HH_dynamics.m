function F = HH_dynamics(~, stateMat, paramStruct)
% state should be a matrix

%% Default parameter values
p.C = 0.9;
p.ENa = 50;
p.EL = -70;
p.gNa = 100;
p.gL = 0.25;
p.gB = 3.5;
p.EB = -74.8;
p.VBth = -2.2;
p.SB = 9.6;
p.tauB = 64;
p.I = 2;
p.EK = -95;
p.gK = 7;
p.mNoise = 0.1;

%% Parse input
V = stateMat(1, :);
n = stateMat(2, :);
h = stateMat(3, :);
B = stateMat(4, :);

for f = fieldnames(paramStruct)'
	p.(f{:}) = paramStruct.(f{:});
end

%% Set parameters

mInf = 1./(1 + exp((-V-34.5)/10));	% HH options
nInf = 1./(1 + exp((-V-29.5)/10));
hInf = 1./(1 + exp((V+59.4)/10.7));
tauN = 0.25 + 4.35*exp(-abs(V+10)/10);
tauH = 0.15 + 1.15./(1 + exp((V+33.5)/15));
BInf = 1./(1 + exp(-(V-p.VBth)./p.SB));
tauB = p.tauB;

%% Calculate changes
Vdot = (... % F1(V, n, h, B)
	p.I - p.gK * n.^4 .* (V - p.EK) ... % drive current minus Potassium current
	- p.gNa * mInf.^3 .* h .* (V - p.ENa) ... % Sodium current
	- p.gB .* B .* (V - p.EB) ... % mystery current
	- p.gL .* (V - p.EL)) / p.C; ... % Leak

ndot = (nInf - n) ./ tauN; % F2(n)

hdot = (hInf - h) ./ tauH; % F3(h)

Bdot = (BInf - B) ./ tauB; % F4(B)

F = [Vdot; ndot; hdot; Bdot]; % state change

end
