function bounded = keep_in_bounds(data, bounds)
% hit the wall

bounded = data;

fn = fieldnames(bounds)';
for i = 1:numel(fn)
	f = fn{i};
	b = bounds.(f);
	if ~b(3), continue; end  % if parameter is not bounded continue
	b(1:2) = (b(1:2) - mean(b(1:2))) * b(3) + mean(b(1:2));  % adjust range
	d = data.(f);

	indsHi = d > b(2);
	indsLo = d < b(1);
	% High values move inside range
	d(indsHi) = b(2) - abs(randn(1, sum(indsHi)) * .02 * diff(b(1:2)));  
	% ... and Low values
	d(indsLo) = b(1) + abs(randn(1, sum(indsLo)) * .02 * diff(b(1:2)));
	bounded.(f) = d;
end

end