function  [BurstMovie]=SimulateNetwork2D(P,varargin)
%SIMULATENETWORK2 This function animates neurons firing in theta phi space.
%
%First load your desired P struct.
%If desired, specify a NumberOfConnections,a NodeSizeParameter, a NodeBrightnessParameter, an EdgeParameter, a number >= 1 for the StartTime,
%and a number <= 600,000 for the EndTime.
%
% Ex. SimulateNetwork2D(P,'NumberOfConnections',1000,'NodeSizeParameter',P.NetworkProperties{1}.PositionSphericalCoordinates(:,1),'EdgeParameter',P.Connectivity,'StartTime',3000,'EndTime',3800)

%% Parsing Variables

for i=1:2:length(varargin);
    switch varargin{i}
        case 'NodeSizeParameter'
            NodeSizeParameter=varargin{i+1};
        case 'NodeBrightnessParameter'
            NodeBrightnessParameter=varargin{i+1};
        case 'EdgeParameter'
            EdgeParameter=varargin{i+1};
        case 'StartTime'
            StartTime=varargin{i+1};
        case 'EndTime'
            EndTime=varargin{i+1};
        case 'NumberOfConnections'
            NumberOfConnections=varargin{i+1};
        case 'IsolateNeuron'
            IsolateNeuron=varargin{i+1};
        case 'HighlightConnections'
            HighlightConnections=varargin{i+1};
        case 'BurstRaster'
            BurstRaster=varargin{i+1};
            [InputOutputMatrix, Time]=ConnectBurst(P, BurstRaster);
    end
end

%% Default Variables
Theta=P.NetworkProperties{1,1}.PositionSphericalCoordinates(:,1);
Phi=P.NetworkProperties{1,1}.PositionSphericalCoordinates(:,2);

NumberOfNeurons=P.NetworkProperties{1,1}.NumberOfNeurons;
NumberOfExcitatory=NumberOfNeurons*(1-P.NetworkProperties{1,1}.PercentInhibitoryNeurons);

MinSize=1;
MaxSize=200;

MinBrightness=0.1;
MaxBrightness=1.3;

MinWidth=0.2;
MaxWidth=1;

NodeSizeConstant=1;
NodeBrightnessConstant=1;
EdgeConstant=1;

if exist('StartTime','var')==0
    StartTime=1;
end
if exist('EndTime','var')==0
    EndTime=length(P.Spikes);
end
if exist('NumberOfConnections','var')==0
    NumberOfConnections=1000;
end
if exist('EdgeParameter','var')==0
    EdgeParameter=ones(size(P.Connectivity));
    EdgeConstant=1/MaxWidth;
end

Connectivity=P.NetworkProperties{1,1}.Connectivity;
Position=P.NetworkProperties{1,1}.PositionXYZCoordinates;

%% Creating the Figure
h.Figure=figure('Units','Normalized','Position', [.3, .1, .8, .6],...
    'Name', 'Neuron Firing Simulation');
movegui(h.Figure,'center')
set(gcf,'Color',[0.95,0.98,0.98])
set(h.Figure,'toolbar','figure')
axis([-pi,pi,-pi/2,pi/2])
axis off
h.Title= uicontrol('Style','Text','Units','Normalized','Position',[.1 .87 .8 .1],...
    'BackgroundColor',[0.95,0.98,0.98],'FontName','Candara','FontSize',40,...
    'FontWeight','Light','ForegroundColor',[0 .5 .5],'String','SIMULATED NEURONS');
hold on
%% Plotting the Neurons

if exist('NodeSizeParameter','var')==0 && exist('NodeBrightnessParameter','var')==0
    Excitatory=scatter(Theta(1:NumberOfExcitatory),...
        Phi(1:NumberOfExcitatory),100,'filled','MarkerFaceColor',[0 .3 .75]);
    scatter(Theta(1:NumberOfExcitatory)-2*pi,...
        Phi(1:NumberOfExcitatory),100,'filled','MarkerFaceColor',[0 .3 .75]);
    Inhibitory=scatter(Theta(NumberOfExcitatory+1:NumberOfNeurons),...
        Phi(NumberOfExcitatory+1:NumberOfNeurons),100,'filled','MarkerFaceColor',[145/255;44/255;23/255]);
    scatter(Theta(NumberOfExcitatory+1:NumberOfNeurons)-2*pi,...
        Phi(NumberOfExcitatory+1:NumberOfNeurons),100,'filled','MarkerFaceColor',[145/255;44/255;23/255]);
else
    if exist ('NodeSizeParameter','var')==0
    NodeSizeParameter=ones(1,NumberOfNeurons);
    NodeSizeConstant=100/MaxSize;
    end
    if exist ('NodeBrightnessParameter','var')==0
    NodeBrightnessParameter=ones(1,NumberOfNeurons);
    NodeBrightnessConstant=1/MaxBrightness;
    end
    NodeSizeParameter= abs(NodeSizeParameter);
    Slope=(MaxSize-MinSize)/(max(NodeSizeParameter));
    
    NodeBrightnessParameter=abs(NodeBrightnessParameter);
    Slope2=(MaxBrightness-MinBrightness)/(max(NodeBrightnessParameter));
    
    for Neuron=1:NumberOfExcitatory
        Excitatory= scatter(P.NetworkProperties{1,1}.PositionSphericalCoordinates(Neuron,1),...
            P.NetworkProperties{1,1}.PositionSphericalCoordinates(Neuron,2),(MinSize+Slope*(NodeSizeParameter(Neuron)))...
            *NodeSizeConstant,'Filled','MarkerFaceColor',...
            [0 .3 .75].*(MinBrightness+Slope2*NodeBrightnessParameter(Neuron))*NodeBrightnessConstant);
        scatter(Theta(Neuron)-2*pi,...
            Phi(Neuron),(MinSize+Slope*(NodeSizeParameter(Neuron)))...
            *NodeSizeConstant,'Filled','MarkerFaceColor',...
            [0 .3 .75].*(MinBrightness+Slope2*NodeBrightnessParameter(Neuron))*NodeBrightnessConstant);
    end

    for Neuron=NumberOfExcitatory+1:NumberOfNeurons
        Inhibitory= scatter(P.NetworkProperties{1,1}.PositionSphericalCoordinates(Neuron,1),...
            P.NetworkProperties{1,1}.PositionSphericalCoordinates(Neuron,2),(MinSize+Slope*(NodeSizeParameter(Neuron)))...
            *NodeSizeConstant,'Filled','MarkerFaceColor',...
            [145/255;44/255;23/255].*(MinBrightness+Slope2*NodeBrightnessParameter(Neuron))*NodeBrightnessConstant);
        scatter(Theta(Neuron)-2*pi,...
            Phi(Neuron),(MinSize+Slope*(NodeSizeParameter(Neuron)))...
            *NodeSizeConstant,'Filled','MarkerFaceColor',...
            [145/255;44/255;23/255].*(MinBrightness+Slope2*NodeBrightnessParameter(Neuron))*NodeBrightnessConstant);
    end
    
end
legend([Excitatory,Inhibitory],'excitatory','inhibitory')
set(legend,'FontSize',14,'FontWeight','Bold','Units',...
    'Normalized','Position',[0.8,0.75,0.1,0.08])
%% Plotting the Connections between Neurons

Connectivity= full(Connectivity);
[Postsynaptic,Presynaptic,Strength]=find(Connectivity);

if exist('IsolateNeuron','var')==0
    Index=randperm(length(Postsynaptic),NumberOfConnections);
    Postsynaptic=Postsynaptic(Index);
    Presynaptic=Presynaptic(Index);
    Strength=Strength(Index);
else
    Index=vertcat(find(Postsynaptic==IsolateNeuron), find(Presynaptic==IsolateNeuron));
    Postsynaptic=Postsynaptic(Index);
    Presynaptic=Presynaptic(Index);
    Strength=Strength(Index);
end

[~,~,EdgeParameter]=find(EdgeParameter);
EdgeParameter=EdgeParameter(Index);
EdgeParameter= abs(EdgeParameter);
Slope3=(MaxWidth-MinWidth)/(max(EdgeParameter));


for Connection=1:length(Index)
    LinePoints=[];
    Count=0;
    for t=linspace(0,1,200);
        Count=Count+1;
        LinePoints(Count,:)= Position(Presynaptic(Connection),:)...
            +t*(Position(Postsynaptic(Connection),:)-Position(Presynaptic(Connection),:));
        LinePoints(Count,:)=LinePoints(Count,:)./(sqrt(sum(LinePoints(Count,:).^2)));
    end
    [Th,Ph,~]=cart2sph(LinePoints(:,1),LinePoints(:,2),LinePoints(:,3));
    
    if Strength(Connection)>0
    plot(Th(Th<pi&Th>0&Ph<0),Ph(Th<pi&Th>0&Ph<0),'Color',[0 .2 .75],...
            'LineWidth',(MinWidth+Slope3*(EdgeParameter(Connection)))*EdgeConstant)
    plot(Th(Th<0&Ph<0),Ph(Th<0&Ph<0),'Color',[0 .2 .75],...
            'LineWidth',(MinWidth+Slope3*(EdgeParameter(Connection)))*EdgeConstant)
    plot(Th(Th<pi&Th>0&Ph>0),Ph(Th<pi&Th>0&Ph>0),'Color',[0 .2 .75],...
            'LineWidth',(MinWidth+Slope3*(EdgeParameter(Connection)))*EdgeConstant)
    plot(Th(Th<0&Ph>0),Ph(Th<0&Ph>0),'Color',[0 .2 .75],...
            'LineWidth',(MinWidth+Slope3*(EdgeParameter(Connection)))*EdgeConstant)
    else
    plot(Th(Th<pi&Th>0&Ph<0),Ph(Th<pi&Th>0&Ph<0),'Color',[230/255;44/255;23/255],...
            'LineWidth',(MinWidth+Slope3*(EdgeParameter(Connection)))*EdgeConstant)
    plot(Th(Th<0&Ph<0),Ph(Th<0&Ph<0),'Color',[230/255;44/255;23/255],...
            'LineWidth',(MinWidth+Slope3*(EdgeParameter(Connection)))*EdgeConstant)
    plot(Th(Th<pi&Th>0&Ph>0),Ph(Th<pi&Th>0&Ph>0),'Color',[230/255;44/255;23/255],...
            'LineWidth',(MinWidth+Slope3*(EdgeParameter(Connection)))*EdgeConstant)
    plot(Th(Th<0&Ph>0),Ph(Th<0&Ph>0),'Color',[230/255;44/255;23/255],...
            'LineWidth',(MinWidth+Slope3*(EdgeParameter(Connection)))*EdgeConstant)
    end
end

%% Highlight Connections
if exist('HighlightConnections','var')~=0
for Connection=1:size(HighlightConnections,1)
    LinePoints=[];
    Count=0;
    
    for t=linspace(0,1,100);
        Count=Count+1;
        LinePoints(Count,:)= Position(HighlightConnections(Connection,1),:)...
            +t*(Position(HighlightConnections(Connection,2),:)-Position(HighlightConnections(Connection,1),:));
        LinePoints(Count,:)=LinePoints(Count,:)./(sqrt(sum(LinePoints(Count,:).^2)));
    end
    [Th,Ph,~]=cart2sph(LinePoints(:,1),LinePoints(:,2),LinePoints(:,3));

    plot(Th(Th<pi&Th>0&Ph<0),Ph(Th<pi&Th>0&Ph<0),'Color','g',...
            'LineWidth',2*(MinWidth+Slope3*(EdgeParameter(Connection)))*EdgeConstant)
    plot(Th(Th<0&Ph<0),Ph(Th<0&Ph<0),'Color','g',...
            'LineWidth',2*(MinWidth+Slope3*(EdgeParameter(Connection)))*EdgeConstant)
    plot(Th(Th<pi&Th>0&Ph>0),Ph(Th<pi&Th>0&Ph>0),'Color','g',...
            'LineWidth',2*(MinWidth+Slope3*(EdgeParameter(Connection)))*EdgeConstant)
    plot(Th(Th<0&Ph>0),Ph(Th<0&Ph>0),'Color','g',...
            'LineWidth',2*(MinWidth+Slope3*(EdgeParameter(Connection)))*EdgeConstant)
end
end

%% Firing Animation
if isfield(P, 'Spikes')~=0
    Spikes=P.Spikes;
    Fire=scatter(0,0,90,'filled','Marker','d','MarkerFaceColor','white');
    
    for t=StartTime:EndTime
        %ActiveNeurons=find(Spikes(t,:));
       
        set(Fire,'XData',vertcat(Theta((Spikes(t,:)==1)'&Theta<pi),Theta((Spikes(t,:)==1)'&Theta>pi)-2*pi),'YData',...
            vertcat(Phi((Spikes(t,:)==1)'&Theta<pi),Phi((Spikes(t,:)==1)'&Theta>pi)));
        pause(0.02)
    end
end
end






