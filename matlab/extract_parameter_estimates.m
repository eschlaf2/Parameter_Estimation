% function [] = extract_parameter_estimates()

result = struct();  % store all results in a structure
load('1', 'spiketimes');

% load each trial and store the results in res
fprintf(' 0%%\n');
N = size(dir('*.mat'), 1);
for ii = 1:N
	fname = num2str(ii);
	try
		load(fname, 'estimates', 'res');	
		for f = fieldnames(estimates.params)'
			param = f{:};
			result.(param){ii} = estimates.params.(param);
		end
		result.compare{ii} = res;
		clear estimates 
		fprintf('\b\b\b\b%02d%%\n', round(ii / N * 100))
	catch ME
% 		display(['Skipping file: ' fname])
		fprintf('\b\b\b\bSkipping file: %s\n%02d%%\n', fname, round(ii / N * 100))
	end
	
end

save('res')

% Clean up the following to visualize results in 3D

% figure(990);
% for i = 1:900
% scatter3(T.Var1(i, :), T.Var2(i, :), T.Var3(i, :), 1000 * [res.compare{1, 1}.p res.compare{1, 2}.p res.compare{1, 3}.p], lines(3), 'o')
% xlim([min(T.Var1(:)), max(T.Var1(:))]);
% ylim([min(T.Var2(:)), max(T.Var2(:))]);
% zlim([min(T.Var3(:)), max(T.Var3(:))]);
% drawnow; grid on;
% pause(.01)
% end