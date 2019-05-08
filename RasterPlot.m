function RasterPlot(Spikes)
%RasterPlot plots the spikes of neurons with neurons on the Y axis and time
%steps on the X axis.
%
%Ex. RasterPlot(P.Spikes)

%%
[row,column]=find(Spikes);
h=plot(row,column,'.','Markersize',0.08);
set(gcf,'Color',[0.95,0.98,0.98])
xlabel('Time Steps','FontSize', 20,'FontName','Candara','Color',[0 0.4 0.4]);
ylabel('Neuron   number','FontSize', 20,'FontName','Candara','Color',[0 0.4 0.4]);
title('RASTER PLOT ','FontSize', 30,'fontWeight','bold','FontName','Candara','Color',[0 0.4 0.4]);
set(h,'color',[30/255,144/255,235/255]);
set(gca,'Ydir','reverse')
    
end

