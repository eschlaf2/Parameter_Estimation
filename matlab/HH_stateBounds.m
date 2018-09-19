function [s0, boundsStruct] = HH_stateBounds()
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
% I = 2;
% EK = -95;
% gK = 7;


paramBounds = {...  % 'parameter', [lowerBound, upperBound, rangeExpand, procNoise]
				'gB',	[0,	10,	1, .02];	... % (0 10)
				'EB',	[-90, 110, 2, .02];	...	% (-110 110)
				'VBth', [-95, 5, 2, .02];	... % (-95 5)
				'SB',	[-10, 15, 1, .02];	... % (-10 10)
% 				'tauB', [0, 80,	1, .02];	... % (0 80)
				'I',	[-5, 5,	1, .05];	... % (-5 5)
% 				'mNoise', [-2,	2,		2, .25];	... 
				'EK',	[-115, -50, 1, .02]; ...
				'gK',	[0, 15, 1, .02]; ...
				'EL',	[-90, -20, 1, .02]; ...
				};

boundsStruct = cell2struct(paramBounds(:, 2), paramBounds(:, 1));

end
