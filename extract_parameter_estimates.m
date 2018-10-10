% function [] = extract_parameter_estimates()

res = struct();  % store all results in a structure
load('1', 'spiketimes');

% load each trial and store the results in res
for ii = 1:3  
	fname = num2str(ii);
	load(fname, 'estimates', 'stEst', 'dt');	
	for f = fieldnames(estimates.params)'
		param = f{:};
		res.(param){ii} = estimates.params.(param);
	end
	res.compare{ii} = compare_spiketimes(spiketimes, stEst, dt);
	clear estimates stEst dt
end


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