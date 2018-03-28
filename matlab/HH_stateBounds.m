function [s0, stateBounds] = HH_stateBounds()
% STATES = (V, n, h, B, gB, EB, VBth, SB, tauB, I) 

s0 = [	-71;	... % V
		0.0147; ...	% n
		0.7497; ...	% h
		0.0326; ...	% B
		];
	
% gB = 3.5;
% EB = -74.8;
% VBth = -2.2;
% SB = 9.6;
% tauB = 64;
% I = -2.5;

stateBounds = [	3.5,		3.5,		1;	...	% gB
				-74.8,	-74.8,	1;	...	% EB
				-2.2,	-2.2,		1;	...	% VBth
				9.6,	9.6,		1;	... % SB
				64,		64,		1;	... % tauB
				-5,		5,		1;	... % I
				0.,		.5,		1;	... % mNoise
				];

end

% stateBounds = [	0,		10;		...	% gB
% 				-110,	110;	...	% EB
% 				-95,	5;		...	% VBth
% 				-10,	10;		... % SB
% 				0,		80;		... % tauB
% 				-10,	0;		... % I
% 				];
