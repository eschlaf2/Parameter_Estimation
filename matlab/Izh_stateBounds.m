function [s0, boundsStruct] = Izh_stateBounds()
% STATES = (V, n, h, B, gB, EB, VBth, SB, tauB, I) 

s0 = [	-50;	... % v
		-10; ...	% u
		];


paramBounds = {...  % 'parameter', [lowerBound, upperBound, rangeExpand, procNoise]
				'a',	[0, .2, 4, .01];	... % ()
				'b',	[0, 1, 4, .01];	...	% ()
				'c',	[-80, -40, 4, .01];	... % ()
				'd',	[-40, 12, 4, .01];	... % ()
% 				'I',	[-5, 80, 4, .05];	... % ()
				};

boundsStruct = cell2struct(paramBounds(:, 2), paramBounds(:, 1));

end
