% function [] = visualize_pf_results()

% load res
% load('1.8783638.mat', 'binwidth', 'simParams')

N = length(res.compare);
p = zeros(N, 1);
acfDiff = p;
phatEst = zeros(N, 2);

for i = 1:N
	p(i) = res.compare{i}.p;
	acfDiff(i) = res.compare{i}.isiAcfDiff;
	phatEst(i, :) = res.compare{i}.phat2;
end
pReject = p < .05;
acfReject = acfDiff > 0;
reject = pReject - 2 * acfReject;

pScaled = p / sum(p);

fn = fieldnames(res)';

figure(); fullwidth(true);
i = 1;
for f = fn
	if strcmp(f{:}, 'compare')
		continue
	end
	resMat.(f{:}) = cell2mat(res.(f{:}));
	T = size(resMat.(f{:}), 1);
	subplot(size(fn, 2) - 1, 1, i); i = i + 1;
% 	plot(resMat.(f{:})(:, pReject), 'r'); hold on;
% 	plot(resMat.(f{:})(:, acfReject), 'b'); hold on;
	plot(resMat.(f{:})(:, reject == 0), 'color', .5 * [1 1 1]); hold on;
	ax3 = plot(simParams.(f{:})(1:binwidth:T*binwidth), 'k--', 'linewidth', 3, ...
		'Displayname', 'Actual');
% 	ax1 = plot(sum(resMat.(f{:}) .* repmat(pScaled', T, 1), 2), 'linewidth', 3, ...
% 		'Displayname', 'Weighted mean');
	ax2 = plot(median(resMat.(f{:}), 2), 'linewidth', 3, ...
		'Displayname', 'Median Estimate'); 
	hold off
	title(f{:})
end
legend([ax3, ax2])

figure(); stem(reject); fullwidth()