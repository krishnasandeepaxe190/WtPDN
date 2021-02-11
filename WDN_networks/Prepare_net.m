%% This code extract the relevant paramaters for water side
clear all
clc
close all
disp('Availabel water_nets 3, 8')
Net_num = input('Enter the water_net network num:')
if Net_num ==3
    wd_3 = epanet('Wdn_3node.inp');
else
    wd_8 = epanet('tutorial8node_noeff_PRV.inp');
    filename = ['WDN_',num2str(Net_num),'.mat'];
end
%% Extract the relevant parameters from EPANET file
if Net_num == 8
wd = wd_8;
%% Get Reservoir Index
NodeReservoirIndex = (wd.NodeReservoirIndex); 
%% Get Valave Index and setting
LinkValveIndex = wd.LinkValveIndex; %Link having valves
LinkInitialSetting = wd.LinkInitialSetting; 
%% Get Link Name ID and other parameters 
LinkNameID = wd.LinkNameID;
LinkDiameter = wd.LinkDiameter; % in inches
LinkLength = wd.LinkLength; % Link length in ft 
LinkRoughnessCoeff = wd.LinkRoughnessCoeff; % Unit less
LinkMinorLossCoeff = wd.LinkMinorLossCoeff; % Unit less
%% Add link minor loss coeffcient 
LinkResistance  = zeros(1,length(LinkRoughnessCoeff));
for i = 1:length(LinkRoughnessCoeff)
    LinkResistance(i) = 10.4622*LinkRoughnessCoeff(i)^(-1.852)*(LinkDiameter(i))^(-4.871)*LinkLength(i);
end
LinkMinorResistance  = zeros(1,length(LinkMinorLossCoeff));
for i = 1:length(LinkMinorLossCoeff)
    LinkMinorResistance(i) = (0.002596.*LinkMinorLossCoeff(i))/(LinkDiameter(i).^4);
end
%% Get Node Index and Junction Index and elevations
NodeIndex = wd.NodeIndex;
JunctionIndex = wd.NodeJunctionIndex; 
NodeElevations = wd.NodeElevations; % in feet 
%NodeElevations = [700 700 710 700 650 700 700 830]; % Need to replace following EPANET Manual
%% Junction Demands
JunctionDemand = cell2mat(wd.NodeBaseDemands); 
%JunctionDemand = [0 0 150 150 200 150 150 200]; % in gpm
%DemandProfile  = [0.8 0.8 0.7 0.6 0.6 0.5 0.6 0.7 0.7 0.9 0.9 1.0]; 
DemandProfile  = [1 0.8 0.7 0.6 0.6 0.5 0.6 0.7 0.7 0.9 0.9 1.0 ]; 
%% Creating Junction profile for 8:00 AM to 8:00 PM 
JunctionProfile = zeros(length(JunctionDemand),length(DemandProfile));
for i = 1:length(JunctionDemand)
    JunctionProfile(i,:) = JunctionDemand(i)*DemandProfile;
end
%% Get Tank related parameters
TankIndex = wd.NodeTankIndex;
NodeTankDiameter = wd.NodeTankDiameter; % in ft
NodeTankMinimumWaterLevel = wd.NodeTankMinimumWaterLevel+NodeElevations(TankIndex);
NodeTankMaximumWaterLevel = wd.NodeTankMaximumWaterLevel+NodeElevations(TankIndex);
NodeTankInitialWaterLevel = wd.NodeTankInitialLevel+NodeElevations(TankIndex); 
A_tk = zeros(1,length(TankIndex));
for i = 1:length(TankIndex)
    A_tk(i) = 3.1415*(NodeTankDiameter(i)/2)^2;
end
%% Get Pump related parameters 
LinkPumpCount = wd.LinkPumpCount; 
LinkPumpIndex = wd.LinkPumpIndex; 
% Pump parameters 
%h_0 = 200 ; % Pump shut off head in ft
h_0 = 266.67;
r_m = 4.629587e-05;
v_m = 2; % pump curve coefficient from EPANET
c_m = 0.7457/(3960*0.8); % pump power coefficient in Kw (Q*H/3960*eff)*0.745

%% Get Node connecting Link Ids
NodesConnectingLinksID = double(wd.NodesConnectingLinksIndex); 
%NodeNameID = wd.NodeNameID;
%LinkNameIndex = zeros(1,length(LinkNameID));
%for i = 1:length(LinkNameID)
%    LinkNameIndex(i) = str2double(LinkNameID{i});
%end
%FromNodeIndex = zeros(length(LinkNameID),1);
%ToNodeIndex = zeros(length(LinkNameID),1);
%for i = 1:length(LinkNameID)
% FromNodeIndex(i) = str2double(NodesConnectingLinksID{i});
% ToNodeIndex(i) = str2double(NodesConnectingLinksID{i,2});
%end
FromNodeIndex = NodesConnectingLinksID(:,1);
ToNodeIndex = NodesConnectingLinksID(:,2); 
FromNodeIndexNew = zeros(length(FromNodeIndex),1);
ToNodeIndexNew = zeros(length(ToNodeIndex),1);
if NodeReservoirIndex~=1
for i = 1:length(FromNodeIndex)
    if (FromNodeIndex(i)~=NodeReservoirIndex && FromNodeIndex(i)>ToNodeIndex(i))
    FromNodeIndexNew(i) = ToNodeIndex(i);
    ToNodeIndexNew(i) = FromNodeIndex(i);
   else
       FromNodeIndexNew(i) = FromNodeIndex(i);
       ToNodeIndexNew(i) = ToNodeIndex(i);
    end
end
else 
   for i = 1:length(FromNodeIndex)
    if FromNodeIndex(i)>ToNodeIndex(i)
    FromNodeIndexNew(i) = ToNodeIndex(i);
    ToNodeIndexNew(i) = FromNodeIndex(i);
   else
       FromNodeIndexNew(i) = FromNodeIndex(i);
       ToNodeIndexNew(i) = ToNodeIndex(i);
    end
   end
end
else
wd = wd_3; 
LinkNameID = wd.LinkNameID;
LinkDiameter = wd.LinkDiameter; % in inches
LinkLength = wd.LinkLength; % Link length in ft
LinkRoughnessCoeff = wd.LinkRoughnessCoeff; % Unit less
LinkResistance  = zeros(1,length(LinkRoughnessCoeff));
for i = 1:length(LinkRoughnessCoeff)
    LinkResistance(i) = 10.5*LinkRoughnessCoeff(i)^(-1.852)*(LinkDiameter(i))^(-4.874)*LinkLength(i);
end
NodeIndex = wd.NodeIndex;
JunctionIndex = wd.NodeJunctionIndex; 
NodeElevations = wd.NodeElevations; % in feet 
TankIndex = wd.NodeTankIndex; 
LinkPumpCount = wd.LinkPumpCount; 
LinkPumpIndex = wd.LinkPumpIndex;
NodesConnectingLinksID = wd.NodesConnectingLinksID; 
NodeTankDiameter = wd.NodeTankDiameter; % in ft
NodeTankMinimumWaterLevel = wd.NodeTankMinimumWaterLevel+NodeElevations(TankIndex);
NodeTankMaximumWaterLevel = wd.NodeTankMaximumWaterLevel+NodeElevations(TankIndex);
NodeTankInitialWaterLevel = wd.NodeTankInitialLevel+NodeElevations(TankIndex); 
A_tk = zeros(1,length(TankIndex)); % Area of Tank in Ft.^2
for i = 1:length(TankIndex)
    A_tk(i) = 3.14*(NodeTankDiameter(i)/2)^2;
end
%% Pump parameters 
h_0 = 333.34 ; % shut off head in ft
r_m = 3.754e-05; % pump curve coefficient from EPANET
v_m = 2; % pump curve coefficient from EPANET
c_m = 0.7457/(3960*0.8); % pump power coefficient in Kw (Q*H/3960*eff)*0.745
%% Junction Demands
JunctionDemand = 150;% in gpm
%DemandProfile  = [1 0.8 0.7 0.6 0.6 0.5 0.6 0.7 0.7 0.9 0.9 1.0 ]; 
%% Creating Junction profile for 8:00 AM to 8:00 PM 
JunctionProfile = zeros(length(JunctionDemand),length(DemandProfile));
for i = 1:length(JunctionDemand)
    JunctionProfile(i,:) = JunctionDemand(i)*DemandProfile;
end
end
cd('..');
if exist('PumpMatFiles')~=7
    mkdir 'PumpMatFiles'
end
cd('PumpMatFiles'); 
save(filename);
cd('..');