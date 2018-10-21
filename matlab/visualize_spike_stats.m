% visualize_spike_stats


load firings.mat	% sorted and curated firings from binary_small.dat (preictal); 
							% sorting by MountainSort
spiketimes = firings(2,firings(3,:) == 95);	% spike times from unit 95
fs = 3e4;	% sampling frequency [Hz]
st = spiketimes * 1e3 / fs;  % convert firing times to ms
st = st(st < 60e3);


%% BSIs
figure(1); fullwidth(true)
bsi = zeros(ceil(st(end)), 1);
t = 1:length(bsi);  % time (ms)
for i = 1:length(st)
	bsi(floor(st(i))) = bsi(floor(st(i))) + 1;
end
subplot(7, 4, 1:4)
stem(t(bsi > 0), bsi(bsi > 0), '.')
title('BSI (1 ms)')
xlabel('Time (ms)')
ylabel('Spike count'); yticks(1);
axis('tight')


%% histogram of ISIs
isi = diff(st);
subplot(7, 4, [5, 6, 9, 10] + 4)
histogram(isi, 0:10:max(isi));
title('ISIs (10 ms bins)')
ylabel('Counts')
xlabel('ISI length (ms)')

subplot(7, 4, [7 8 11 12] + 4)
histogram(isi, 0:1:200);
title('ISIs (1 ms bins)')
ylabel('Counts')
xlabel('ISI length (ms)')

%% MLE

if 0  % ignore this section
	phat = mle(isi, 'distribution', 'inverse gaussian');
	pd = fitdist(isi', 'inverse gaussian');
	x_values = 1:10:max(isi);
	y = pdf(pd,x_values);
	subplot(7, 4, [5 6 9 10] + 4); hold on;
	plot(x_values,y,'LineWidth',2); hold off;

	[f, x] = ecdf(isi);
end

%% ACF

nlags = 100;
[acf, lags] = xcorr(bsi - mean(bsi), nlags, 'coeff');
subplot(7, 4, [21 22 25 26])
stem(lags(lags > 0), acf(lags > 0))
hold on
plot([1 nlags; 1 nlags]', 1.96 / sqrt(length(bsi)) * [1 1; -1 -1]', 'r--')
hold off
title('ACF of increments')
xlabel('Lags (ms)')
ylabel('CC')

[acfI, lagsI] = xcorr(isi - mean(isi), nlags, 'coeff');
subplot(7, 4, [23 24 27 28])
stem(lagsI(lagsI > 0), acfI(lagsI >0));
hold on
plot([1 nlags; 1 nlags]', 1.96 / sqrt(length(isi)) * [1 1; -1 -1]', 'r--')
hold off
title('ACF of ISIs')
xlabel('Lags')
ylabel('CC')
