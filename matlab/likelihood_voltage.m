function likelihood = likelihood_voltage(window, obsn, sigma)

%%
likelihood = arrayfun(@(part) ...
	median((1 ./ (sqrt(2 * pi * sigma^2)) .* ...
	exp(-(squeeze(window(1, part, :)) - obsn(:)).^2 / (2 * sigma^2)))), ...
	1:size(window, 2));