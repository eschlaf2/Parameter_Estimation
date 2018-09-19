function p = default_params(model)
%% Default parameter values

switch model
	case 'Izh'
		p.a = 0.02;
		p.b = 0.2;
		p.c = -65;
		p.d = 8;
		p.I = 10;
		p.mNoise = 0;
	case 'HH'
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
		p.mNoise = 0;
end