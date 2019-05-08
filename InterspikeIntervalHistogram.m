function InterspikeIntervalHistogram(P)
%InterspikeIntervalHistogram This function creates a histogram colormap of
%interspike interval of all neurons in the network over a course of time.
%   

%%
for i=1:1000
    SpikeTimePoints=find(P.Spikes(:,i));
    InstantaneousFiringPeriod=(SpikeTimePoints(2:end)-SpikeTimePoints(1:end-1));
    HistogramMatrix(i,1:100)=hist(InstantaneousFiringPeriod,100);
end

h.Figure=figure('Name', 'Interspike Interval Histogram');
set(gcf,'Color',[0.95,0.98,0.98])
imagesc(HistogramMatrix(:,2:100));
title('Insterspike Interval Histogram','FontName','Candara','FontSize',25,...
    'FontWeight','Bold','Color',[0 .5 .5])

xlabel('Histogram Bins','FontName','Candara','FontSize',20,...
    'FontWeight','Light','Color',[0 .5 .5])
ylabel('Neuron','FontName','Candara','FontSize',20,...
    'FontWeight','Light','Color',[0 .5 .5])
colorbar('EastOutside')
end

