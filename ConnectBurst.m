function [InputOutputMatrix, Time]=ConnectBurst(P, Burst)
%ConnectBurst This function determines which neurons are connected to each
%other within 15 time steps

Connectivity=full(P.Connectivity);
NumberOfNeurons=P.NetworkProperties{1,1}.NumberOfNeurons;
Presynaptic=[];
Postsynaptic=[];
Time=[];
for t=1:size(Burst,1)-14
    [~,Neuron1]=find(Burst(t,1:NumberOfNeurons));
    [~,Neuron2]=find(Burst(t+1:t+14,1:NumberOfNeurons));
    [Index2,Index1]=find(Connectivity(Neuron2,Neuron1));
    Index1=reshape(Index1,1,numel(Index1));
    Presynaptic=vertcat(Presynaptic,Neuron1(Index1)');
    Postsynaptic=vertcat(Postsynaptic,Neuron2(Index2));
    Time=vertcat(Time,t*ones(1,numel(Index1))');
end

InputOutputMatrix=[Presynaptic,Postsynaptic];

end