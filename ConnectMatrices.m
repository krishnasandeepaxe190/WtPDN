% This code setups the power distribution system parameters to solve WDN
% problem 
clear
clc
close all
%% Loading grid side parameters
disp('Available cases, 13, 33')
PDN_num = input('Enter the PDN case file number:')
cd('SinglePhaseMatfiles')
if PDN_num == 13
PDN_13 = load('IEEE13SinglePhaseall.mat'); % Single phase version of IEEE 13 bus feeder
else
PDN_33 = load('Case33BWSinglePhaseall.mat'); % Case 33 BAran & Wu
end
cd('..')
%% Loading water side parameters
disp('Available WDN, 3,8')
WDN_num = input('Enter the WDN case file number:')
cd('PumpMatFiles')
if WDN_num == 3
    WDN_3 = load('WDN_3.mat');
else
    WDN_8 = load('WDN_8.mat');
end
cd('..')
%% Form the following connection matrices
if PDN_num==13
    PDN = PDN_13;
    clear PDN_13
%% Grid Side
%1. Psi : Bus to load connection matrix NBXNc
%2. Gamma: Bus to PV connection matrix  NBXNpv
%3. A : Branch to Bus incidence matrix (We get from A_telda. See Kekatos
%voltage Regulation) NbXNB
NB = PDN.N;
Psi = eye(NB);
Noload_ids = find(PDN.p==0); % Bus ids without loads
Psi(:,Noload_ids)=[]; 
Gamma = eye(NB); 
Nopv_ids = find(PDN.pv_ids==0); %Bus without PVs
Gamma(:,Nopv_ids)=[];
A = PDN.A; 

end
%% Water Side
%3. Pi:  node to Pipe/Pump mapping matrix
%4  Lambda : Pump to link mapping matrix
%5  Tau :node to tank mapping matrix TK X N
%6  Theta: Reservoir to node mapping matrix RXN
%7  Omega : Diagonal matrix contains the resistance value for each pipe (P
%X P)
%8 Delta : Surface area of each tank in ft.^2 (Tk X Tk)
%9 Zeta : Valve to link mapping matrix  
if WDN_num == 8
    WDN = WDN_8;
    clear WDN_8
 N = length(WDN.NodeIndex); % Number of Nodes (includes junctions, Tanks and Resrvoirs)
 %% Tau
 NprimeTank = zeros(N,1);
 Tau = eye(N);
 TankIndex = WDN.TankIndex;
 NprimeTank([TankIndex])=1;
 NoTankIndex = find(NprimeTank==0);
 Tau(:,NoTankIndex)=[]; 
 
 %% Theta
 NprimeReservoir = zeros(1,N);
 Theta = eye(N);
 ReservoirIndex = WDN.NodeReservoirIndex; % This has to be provided 
 NprimeReservoir([ReservoirIndex])=1;
 NoReservoirIndex = find(NprimeReservoir==0);
 Theta(:,NoReservoirIndex)=[];
 Theta = Theta';
 %% Lambda 
 NprimePump = zeros(1,length((WDN.NodesConnectingLinksID)));
 Lambda = eye(length(NprimePump));
 LinkPumpIndex = WDN.LinkPumpIndex;
 NprimePump([LinkPumpIndex])=1;
 NoPumpIndex = find(NprimePump==0);
 Lambda(:,NoPumpIndex)=[];
 Lambda = Lambda';
 %% Omega is arc X arc diagonal matrix 
 LinkResistance = WDN.LinkResistance; 
 Omega = diag(LinkResistance); 
 Omega(isnan(Omega))=0; % setting the pump & valve pipe resistance to zero. may need to change  
 %% Omega for PRV
 LinkMinorResistance = WDN.LinkMinorResistance;
 Omega_prv = diag(LinkMinorResistance);
 Omega_prv(isnan(Omega_prv))=0;
 %% Pi Arc to Node mapping matrix P X N
 NodeIndex = WDN.NodeIndex;
 FromNodeIndexNew = WDN.FromNodeIndexNew;
 ToNodeIndexNew = WDN.ToNodeIndexNew; 
 NodesConnectingLinksIndex= [FromNodeIndexNew ToNodeIndexNew];
 Pi  = zeros(length(FromNodeIndexNew),length(NodeIndex));
 for i = 1:length(FromNodeIndexNew)
 Pi(i,NodesConnectingLinksIndex(i,1)) = 1;
 Pi(i,NodesConnectingLinksIndex(i,2))= -1;
 end
 Pi = Pi';
 %Pi_telda = Pi(:,2:end); %Taking reservoir node out 
 %% Delta 
 Delta = diag(WDN.A_tk);
 %% Zeta  
 NprimeValve = zeros(1,length((WDN.NodesConnectingLinksID)));
 Zeta = eye(length(NprimeValve));
 LinkValveIndex = WDN.LinkValveIndex;
 NprimeValve([LinkValveIndex])=1;
 NoValveIndex = find(NprimeValve==0);
 Zeta(:,NoValveIndex)=[];
 %% Connection Matrix for Grid and Water 
 % Need to define a connection matrix coupling grid and WDN Pump 
 % Find the max power of the pump and then couple it to the node <= nominal
 % The Pumps are connected to nodes which hosts loads (Assumption) 
BusAvailabe = find(PDN.p~=0)
disp('Choose the Bus to connect pump based on the maximum power of the Pump')
PuBusIndex = [0 0 0 0 0 0 1 0 0 0 0 0]'; 
Xi = eye(NB);
NoPu_ids = find(PuBusIndex==0); % Bus ids without loads
Xi(:,NoPu_ids)=[]; 
end
filename= ['ConnectMatPDN_',num2str(PDN_num),'WDN_',num2str(WDN_num),'.mat'];
if exist('ConnectMatFiles')~=7
    mkdir 'ConnectMatFiles'
end
cd('ConnectMatFiles'); 
save(filename);
cd('..')


 