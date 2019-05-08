function ActivityRateHeatMap(NeuronX,NeuronY,NeuronZ,NodeParameter,Resolution )
%ACTIVITYRATEHEATMAP This Function creates a heat map based on any given
%node parameter at the given Node Points.
%   Imput the X,Y,Z coordinates of your points, a NodeParameter variable that is a column of values
%   corresponding to each point, and a desired Resolution.
%
%   Ex.
%   ActivityRateHeatMap(P,P.NetworkProperties{1,1}.PositionSphericalCoordinates(:,2),50)

%% 
% These are the X,Y,Z coordinates you
% NeuronX=P.NetworkProperties{1,1}.PositionXYZCoordinates(:,1);
% NeuronY=P.NetworkProperties{1,1}.PositionXYZCoordinates(:,2);
% NeuronZ=P.NetworkProperties{1,1}.PositionXYZCoordinates(:,3);

[A,B] = meshgrid(linspace(0,2*pi,2*Resolution)',linspace(-pi/2,pi/2,Resolution)');
HeatPoints=reshape(cat(2,A',B'),[],2);
[HeatPointX,HeatPointY,HeatPointZ]=sph2cart(HeatPoints(:,1),HeatPoints(:,2),1);

%%
for i=1:length(NeuronX)
    C(i,:)=1./(2*asin(sqrt((NeuronX(i)-HeatPointX).^2+(NeuronY(i)-HeatPointY).^2+(NeuronZ(i)-HeatPointZ).^2)/2));
    [k]=find(C(i,:)-mean(C(i,:))>=4*std(C(i,:)));
    C(i,k)=mean(C(i,:))+4*std(C(i,:));
    D(i,:)=NodeParameter(i).*C(i,:);
end

HeatValues=mean(D,1)';
HeatMapMatrix=flipud(reshape(HeatValues,[2*Resolution,Resolution])');

h.Figure=figure('Name', 'Neuron Heat Map');
set(gcf,'Color',[0.95,0.98,0.98])
imagesc(HeatMapMatrix);
title('Activity Rate Heat Map','FontName','Candara','FontSize',35,...
    'FontWeight','Light','Color',[0 .5 .5])
axis equal
axis off
xlabel('Theta','FontName','Candara','FontSize',20,...
    'FontWeight','Light','Color',[0 .5 .5])
ylabel('Phi','FontName','Candara','FontSize',20,...
    'FontWeight','Light','Color',[0 .5 .5])
end

