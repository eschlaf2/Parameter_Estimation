function [] = plot_sim(sim, spiketimes, dt, Vth, figNum)

%%
if ~exist('dt', 'var') || isempty(dt)
	dt = 1e-2;
end
if ~exist('figNum', 'var') || isempty(figNum)
	figNum = 99;
end
if ~exist('Vth', 'var') || isempty(Vth)
	Vth = 30;
end
t = (1:max(length(sim), spiketimes(end))) * dt;

figure(figNum); fullwidth() % Voltage
ax = subplot(8,1,2:7);
if ~isempty(sim)
    
    plot(t, sim(1,:)); ylabel('Voltage [mV]'); xlabel('Time [ms]'); 
    hold on; plot([0 t(end)], Vth * [1 1], '--', 'color', .5 * [1 1 1]); hold off
    
    hold on; plot(t(spiketimes), sim(1,spiketimes), 'r*'); hold off;
end
    subplot(8,1,1); % Spike raster
	c = lines(2);
    spikeT = t(spiketimes);
    plot([spikeT(:) spikeT(:)]', [zeros(size(spikeT(:))) ones(size(spikeT(:)))]', 'color', c(2,:))
    xlim(get(ax, 'xlim'))
    yticks([]); xticks([]);