function result = compareSpiketimes(st1, st2, dt, outfile)

%%
ISI1 = sort(diff(st1 * dt));
ISI2 = sort(diff(st2 * dt));

[~, result.p, result.k] = kstest2(ISI1, ISI2);

if ~exist('isiFig', 'var')
	isiFig = figure();
	acfFig = figure();
end

figure(isiFig); clf
subplot(2, 2, 1);
histogram(ISI1, ceil(max(ISI1))); % hold on;
title('Sample 1')
subplot(222);
histogram(ISI2, ceil(max(ISI2))); % hold on;
title('Sample 2')
subplot(2, 2, 3:4)
[f1, x1, flo, fhi] = ecdf(ISI1); hold on;
e1 = plot(x1, f1);
plot(x1, flo, 'b:', x1, fhi, 'b:')
[f2, x2] = ecdf(ISI2); e2 = plot(x2, f2); hold off;
legend([e1, e2], {'Sample 1', 'Sample 2'}, 'location', 'southeast')

% Save figure
if exist('outfile', 'var')
	savefig(isiFig, [outfile '_cdf']);
end
%%
acf = cell(2, 1);
numlags = min(numel(ISI1), numel(ISI2));
for i = [1 2]
	isi = eval(['ISI' num2str(i)]);
	[acf{i}, lags] = xcorr(isi - mean(isi), numlags, 'coeff');
	acf{i} = acf{i}(lags >= 0); lags = lags(lags >= 0);
	
	figure(acfFig);
	subplot(2, 4, 2*i-1:2*i)
	stem(lags, acf{i}); hold on
	line([lags(1) lags(end)], 2 / sqrt(numel(isi)) * [1 1]); 
	line([lags(1) lags(end)], -2 / sqrt(numel(isi)) * [1 1]); hold off;
	xlabel('Lags (spikes)')
	if i == 1; ylabel('ACF'); end
	title(['ACF of ISIs: sample ' num2str(i)]);
end

subplot(2,4,6:7)
acfDiff = acf{1} - acf{2};
stem(lags, acfDiff); hold on
line([lags(1) lags(end)], 2 * sqrt(1/numel(ISI1) + 1/numel(ISI2)) * [1 1]); 
line([lags(1) lags(end)], -2 * sqrt(1/numel(ISI1) + 1/numel(ISI2)) * [1 1]); hold off;
xlabel('Lags (spikes)')
ylabel('Difference')
title('ACF1 - ACF2')

result.acf1 = acf{1};
result.acf2 = acf{2};
result.isiAcfDiff = sum(abs(acfDiff) > 2 * sqrt(1/numel(ISI1) + 1/numel(ISI2)));

% Save figure
if exist('outfile', 'var')
	savefig(acfFig, [outfile '_acf']);
end

%%
bins = 0:max([ISI1(:); ISI2(:)]);
edges = [-.5, (bins + .5)];

[phat1, pci1] = mle(ISI1, 'distribution', 'inverse gaussian');
[phat2, pci2] = mle(ISI2, 'distribution', 'inverse gaussian');
figure(990); clf
histogram(ISI1, edges, 'Normalization', 'pdf'); hold on;
histogram(ISI2, edges, 'Normalization', 'pdf');
set(gca,'ColorOrderIndex',1)
plot(pdf('inverse gaussian', bins, phat1(1), phat1(2)), 'b', ...
	'linewidth', 2);
plot(pdf('inverse gaussian', bins, phat2(1), phat2(2)), 'r:', ...
	'linewidth', 2);
hold off;
legend('ISI 1', 'ISI 2', 'Fit 1', 'Fit 2')

result.phat1 = phat1;
result.phat2 = phat2;
result.pci1 = pci1;
result.pci2 = pci2;

% Save figure
if exist('outfile', 'var')
	savefig(isiFig, [outfile '_cdf']);
end

