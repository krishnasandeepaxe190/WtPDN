% This code solves the WDN problem using GP approach for WDN and verifies
% the Pump power and heads of the nodes with nonlinear solver EPANET 
clear
clc
close all
%% Loading grid side parameters
disp('Available cases, 13, 33')
disp('Available WDN, 3,8')
PDN_num = input('Enter the PDN case file number:')
WDN_num =  input('Enter the WDN case file number:')
if PDN_num==13&&WDN_num==8
cd('ConnectMatfiles')
Connect_Mat = load('ConnectMatPDN_13WDN_8.mat'); 
filename = ['WaterDecoupled_',num2str(WDN_num),'.mat'];
else
end
cd('..')
%% WDN side Quantities
WDN = Connect_Mat.WDN;
Pi = Connect_Mat.Pi; % Node to arc
Lambda = Connect_Mat.Lambda; % Pump to arc
Tau = Connect_Mat.Tau; %Node to tank
Theta = Connect_Mat.Theta; % reservoir to node 
Omega = Connect_Mat.Omega; % arc to arc (diagonal with resistance coefficients for the pipes) 
Delta = Connect_Mat.Delta; % diag(tank_Atk) 
Zeta = Connect_Mat.Zeta;  % Arc-to-valve matrix 
Pi_telda = Pi';
PumpLinkIndex = find(Lambda==1);
ValveLinkIndex = find(Zeta==1); 
%Pi_telda(Lambda==1,:) =[]; % Link to Node Ignoring Pump Link 
%Pi_telda(Zeta==1,:)=[]; % Link to Node Ignoring valve 
Pi_telda([PumpLinkIndex,ValveLinkIndex],:)=[];
Pipes = size(Pi_telda,1); % without valves and pumps whose flows can be changed as needed 
%% Pi_reduced giving juctions only 
Pi_reduced = Pi; 
Pi_reduced([find(Tau==1),find(Theta==1)],:)=[];
%Pi_reduced(Tau==1,:)=[]; % ignoring tank and reservoir 
%Pi_reduced(Theta==1,:) = []; 
%% Pump Parameter
h_0 = WDN.h_0;
r_m = WDN.r_m;
v_m = WDN.v_m ; 
%c_m = WDN.c_m; %commented to include the pump efficiency as the variable
c_m = 0.7457;
c_m2 = 3960*0.9; 
%% efficiency coefficients
Q_bep = 1200; 
a1 = (0.9-0.00001)./Q_bep.^2;
a2 = (2*(0.9-0.00001))./Q_bep;
a3 = 0.00001; 
NodalHeads = (WDN.NodeElevations)'; 
JunctionDemandProfile = WDN.JunctionProfile; 
JunctionDemandProfile = JunctionDemandProfile(:,1:12); 
JunctionDemandProfile([find(Tau==1),find(Theta==1)],:)=[];
Time = size(JunctionDemandProfile,2);
Nodes = length(NodalHeads); % No of Nodes
Links = length(Lambda);    % No of arc (Pump + Pipes+valves)
Tanks = length(find(Tau==1));  % No of Tanks
Pumps = length(find(Lambda==1));  % No of pumps
Reservoirs = length(find(Theta==1)); % No of reservoirs  
Valves = length(find(Zeta==1)); % No of Valves ; 
Pi_prime = eye(Links); 
Pi_prime([find(Lambda==1),find(Zeta==1)],:) = []; 
%% New matrix Kappa 
Kappa = eye(Nodes);
Kappa([find(Tau==1),find(Theta==1)],:)=[];
%Kappa(Tau==1,:) = [];
%Kappa(Theta==1,:) = [];
%disp('Nodes',num2str(Nodes),'Links',num2str(Links),'Tanks',num2str(Tanks),'Pumps',num2str(Pumps),'Reservoirs',num2str(Reservoirs))
%% Nodal Head constraints Repmating all the nodal heads to T instances
MinNodalHeads = repmat(NodalHeads,[1,Time]); 
MaxNodalHeads = repmat(1300,[Nodes,Time]);  

%% Tank State Space Matrices for T instances 
%% This has to written for General Tank State Space 
Atk = ones(Time,1);
del_tk = 1*8.0208; %(ft3/hr*hr) time step converting gpm-ft3/hr because flow is in gpm 
Btk = (del_tk/Delta)*tril(ones(Time,Time));
NodeTankMaximumHead = Connect_Mat.WDN.NodeTankMaximumWaterLevel;
NodeTankMinimumHead = Connect_Mat.WDN.NodeTankMinimumWaterLevel;
NodeTankInitialHead  = Connect_Mat.WDN.NodeTankInitialWaterLevel; 
MaxNodalHeads([find(Tau==1)],:) = repmat(NodeTankMaximumHead,[Tanks,Time]); 
MaxNodalHeads([find(Theta==1)],:) = MinNodalHeads([find(Theta==1)],:);
%% Pump Flow Constraints 
MinFlow = 0; % gpm
MaxFlow = sqrt(h_0/r_m);  % If it has two pumps then we need to put it into vector notation
MinFlowPump = zeros(Links,Time);
MaxFlowPump = repmat(MaxFlow,[Links,Time]); 
%% Initial flows for the links 
%IntFlows = 150*ones(Pipes+Pumps,Time); 
IntFlowPipes = 200*ones(Pipes,Time);
IntFlowPumps = 250*ones(Pumps,Time); 
IntFlowValves = 200*ones(Valves,Time);
IntFlows = [IntFlowPipes;IntFlowPumps;IntFlowValves];
IntHeads = MinNodalHeads; 
IntHeads([find(Tau==1)],:) = repmat(NodeTankInitialHead,[Tanks,Time]);  
%% Pump Inital Binary condition for Big-M 
IntOnOff  = ones(Pumps,Time);
%% Valve Inital Status condition for Big-M  
IntValveStatusOnOff = ones(Valves,Time); %x in the proof ON/OFF
IntValveSubStatey = zeros(Valves,Time); % Open, z=1 in the proof OPEN/ACTIVE
IntValveSubStatez = ones(Valves,Time); % Open,y=0 in the proof OPEN/ACTIVE 
%IntValveHeadLoss = zeros(Valves,Time); % valve Initial Headloss, make sure you satisfy the condition in proof
ValveSetting =repmat(MinNodalHeads([find(Zeta'*Pi'==-1)]),[Valves,Time])+ repmat(WDN.LinkInitialSetting([ValveLinkIndex])*2.30724939,[Valves,Time]); %changing from Psi to ft of head
ValveDownStreamNode = zeros(Nodes,Nodes);
ValveDownStreamNode(find(Zeta'*Pi'==-1),find(Zeta'*Pi'==-1))=1; 


%% Estimates for Pipes based on GP approximation 
% Need to use Omega P X P diagonal matrix here 
Omega([find(Zeta==1)],:) = []; %ignoring row
Omega(:,[find(Zeta==1)]) = [] ; % ignoring coloumn Link to Link ignoring pump Link
Omega([find(Lambda==1)],:) = []; %ignoring row
Omega(:,[find(Lambda==1)]) = [] ; % ignoring coloumn Link to Link ignoring pump Link
IntCp = IntFlowPipes.*(Omega*abs(IntFlowPipes).^(0.852)-ones(Pipes,Time));
%% Estimates for Pump Head Gain   
IntC1M = -h_0.*ones(Pumps,Time); 
IntC2M = r_m.*IntFlowPumps.^(v_m-1);
%% First Order Taylor Approximation for Pump Power
IntPumpPower = c_m*(h_0-r_m.*(Lambda*IntFlows).^2).*(Lambda*IntFlows) ;
IntPumpPrimePower = c_m*(h_0-3*r_m.*(Lambda*IntFlows).^2);
IntAPrimePump = IntPumpPrimePower ;
IntBPrimePump = IntPumpPower-IntPumpPrimePower.*(Lambda*IntFlows); 
%% Initial efficiency and fraction values 
eff_int = c_m2.*(-a1.*(Lambda*IntFlows).^2+a2.*(Lambda*IntFlows)+a3);
y_int  = (IntAPrimePump.*(Lambda*IntFlows)+IntBPrimePump.*IntOnOff)./(eff_int);
%% Inital Variables and saving from Iterates  
IntEps = [IntHeads;IntFlows;IntOnOff;IntValveStatusOnOff;IntValveSubStatey;IntValveSubStatez;y_int];
M_pump = 3e03;
M_valve = 4e03; 
Cp = [] ;
C1M = [];
C2M = []; 
Flows = [];
Heads = [];
OnOff = [] ;
ValveStatusOnOff = [];
ValveSubStatey = [];
ValveSubStatez = []; 
ValveHeadLoss = [] ;
%HeadGain = [] ;
Eps = [] ;
Error = 1; 
n_iter = 1;
Max_iter =400; 
%% Saving Iteration Values 
Eps_iter = [];
Flows_iter = [];
Cp_iter = [];
C1M_iter = [] ;
C2M_iter = [];  
APrimePump_iter = [];
BPrimePump_iter = [];
Ppump_iter = []; 
Ppower_iter = []; % newly added
y_iter = [];   % newly added
eff_iter = [] ;  % newly added 
OnOff_iter = [] ;
ValveStatusOnOff_iter = [];
ValveSubStatey_iter = [];
ValveSubStatez_iter = [] ;
ValveHeadLoss_iter = [];
PumpPower_iter = [];
Obj_iter = [] ;
Error_iter = [] ;
lambda1 = 1;
lambda2 = 0.000001; 
%lambda2 = 0;
while Error >= 0.0001 && n_iter <= Max_iter
if n_iter == 1
    Cp = IntCp;
    C1M = IntC1M;
    C2M = IntC2M;
    Eps_save = IntEps; 
    APrimePump= IntAPrimePump;
    BPrimePump = IntBPrimePump ; 
    y = y_int; 
end
if n_iter<=40
% cvx loop inside while loop 
cvx_begin 
cvx_solver Gurobi
cvx_precision best
variable Flows(Links,Time) 
variable Heads(Nodes,Time)
variable OnOff(Pumps,Time) binary
variable ValveStatusOnOff(Valves,Time) binary
variable ValveSubStatey(Valves,Time) binary
variable ValveSubStatez(Valves,Time) binary
variable ValveHeadLoss(Valves,Time)
variable Ppump(Pumps,Time)
%variable eff(Pumps,Time)  % newley added 
variable Hdummy(Tanks,Time)
%minimize (norm((APrimePump.*(Lambda*Flows)+BPrimePump)',1))
%minimize (2*norm((Ppump)',1)+0.1*(norm(Heads([find(Zeta'*Pi'==-1)],:)-ValveSetting,2)))
%minimize (lambda1*norm((Ppump)',1))
%%minimize (lambda1*norm((Ppump)',1)+lambda2*(vec(Heads([find(Zeta'*Pi'==-1)],:)-ValveSetting).'*vec(Heads([find(Zeta'*Pi'==-1)],:)-ValveSetting)))
minimize (lambda1*sum((Ppump-y.*(c_m2.*(-a1.*(Lambda*Flows).^2+a2.*(Lambda*Flows)+a3)))')+lambda2*(vec(Heads([find(Zeta'*Pi'==-1)],:)-ValveSetting).'*vec(Heads([find(Zeta'*Pi'==-1)],:)-ValveSetting)))
%minimize (0)
subject to
%% Reservoir Operational constraint
Theta*Heads == Theta*MinNodalHeads; 
%% Pump Head Gain Negative
C1M+C2M.*(Lambda*Flows)<=0; 
% Big M for Pump On/Off 
M_pump*(OnOff-1)<=(Lambda*Pi')*Heads-(C1M+C2M.*(Lambda*Flows))<=M_pump*(1-OnOff); 
%Pump Flow Constraints
Lambda*MinFlowPump.*OnOff<=Lambda*Flows<=(Lambda*MaxFlowPump).*OnOff;
% Junction Heads Constraints 
%Kappa*MinNodalHeads<=Kappa*Heads<=Kappa*MaxNodalHeads 
Ppump == (APrimePump.*(Lambda*Flows)+BPrimePump.*OnOff);
%eff <= c_m2.*(-a1.*(Lambda*Flows).^2+a2.*(Lambda*Flows)+a3); % concave right side 
%% Tank State Space 
Hdummy' == Atk*NodeTankInitialHead+Btk*(Flows'*(-Pi)'*Tau); 
% Tank Head Constraint 
(MinNodalHeads'*Tau)<=Heads'*Tau<=(MaxNodalHeads'*Tau); 
%(MinNodalHeads'*Tau)<=Heads'*Tau;
Tau'*Heads == [NodeTankInitialHead, Hdummy(Tanks,1:Time-1)];
%% Pipe Head Loss/Energy Equation  
Pi_telda*Heads == Cp+ Pi_prime*Flows ;
%% Mass Balance only for Junctions ignoring Tanks 
Pi_reduced*Flows  == -JunctionDemandProfile; 
%% Include Valve Constraints 
M_valve*(ValveStatusOnOff-1)<=(Zeta'*Pi'*Heads)-ValveHeadLoss<=M_valve*(1-ValveStatusOnOff); 
(Zeta'*MaxFlowPump).*(ValveStatusOnOff-1)<=Zeta'*Flows<=(Zeta'*MaxFlowPump).*ValveStatusOnOff;
%%0<=ValveSubStatey<=ValveStatusOnOff
%%0<=ValveSubStatez<=ValveStatusOnOff
%%ValveSubStatez+ValveSubStatey<=ones(Valves,Time);
%%M_valve*(ValveSubStatey-1)+ValveSetting<=Heads([find(Zeta'*Pi'==1)],:)<=M_valve*(1-ValveSubStatez)+ValveSetting
%%0<=ValveHeadLoss<=M_valve.*ValveSubStatey;
0<=ValveHeadLoss<=M_valve*ValveStatusOnOff;
%%M_valve*(ValveSubStatey-1)+ValveSetting<=Heads([find(Zeta'*Pi'==-1)],:)<=M_valve*(1-ValveSubStatey)+ValveSetting
cvx_end

%% calculating new values for next iteration. 
Cp = (Pi_prime*Flows).*(Omega*abs(Pi_prime*Flows).^(0.852)-ones(Pipes,Time));
C1M = -h_0.*ones(Pumps,Time);
C2M = r_m.*(Lambda*Flows).^(v_m-1);
PumpPower = c_m*(h_0-r_m.*(Lambda*Flows).^2).*(Lambda*Flows) ;
PumpPrimePower = c_m*(h_0-3*r_m.*(Lambda*Flows).^2);
APrimePump = PumpPrimePower ;
BPrimePump = PumpPower-PumpPrimePower.*(Lambda*Flows); 
eff = c_m2.*(-a1.*(Lambda*Flows).^2+a2.*(Lambda*Flows)+a3);
y = (APrimePump.*(Lambda*Flows)+BPrimePump.*OnOff)./eff;
y(y<0)=0;
Eps = [Heads;Flows;OnOff;ValveStatusOnOff;ValveSubStatey;ValveSubStatez;y];
Obj = cvx_optval;
Error = norm(Eps-Eps_save)
Eps_save = Eps;
n_iter = n_iter+1;
%% saving variables from  iterates 
Flows_iter = [Flows_iter;Flows];
Eps_iter = [Eps_iter;Eps_save];
Error_iter = [Error_iter;Error];
Cp_iter = [Cp_iter;Cp];
C1M_iter = [C1M_iter;C1M];
C2M_iter = [C2M_iter;C2M];
OnOff_iter = [OnOff_iter;OnOff];
Obj_iter = [Obj_iter;Obj];
APrimePump_iter = [APrimePump_iter;APrimePump];
BPrimePump_iter = [BPrimePump_iter;BPrimePump];
PumpPower_iter = [PumpPower_iter;PumpPower];
Ppump_iter = [IntPumpPower;Ppump];
y_iter = [y_iter;y];
eff_iter = [eff_iter;eff];
ValveStatusOnOff_iter = [ValveStatusOnOff_iter;ValveStatusOnOff];
ValveSubStatey_iter = [ValveSubStatey_iter;ValveSubStatey];
ValveSubStatez_iter = [ValveSubStatez_iter;ValveSubStatez] ;
ValveHeadLoss_iter = [ValveHeadLoss_iter;ValveHeadLoss];
else 
cvx_begin 
cvx_solver Gurobi
cvx_precision best
variable Flows(Links,Time) 
variable Heads(Nodes,Time)
variable OnOff(Pumps,Time) binary
variable ValveStatusOnOff(Valves,Time) binary
variable ValveSubStatey(Valves,Time) binary
variable ValveSubStatez(Valves,Time) binary
variable ValveHeadLoss(Valves,Time)
variable Ppump(Pumps,Time)
variable Hdummy(Tanks,Time)
%minimize (norm((APrimePump.*(Lambda*Flows)+BPrimePump)',1))
%minimize (lambda1*norm((Ppump)',1))
%minimize (lambda1*norm((Ppump)',1)+lambda2*(vec(Heads([find(Zeta'*Pi'==-1)],:)-ValveSetting).'*vec(Heads([find(Zeta'*Pi'==-1)],:)-ValveSetting)))
%minimize (0)
minimize (lambda1*sum((Ppump-y.*(c_m2.*(-a1.*(Lambda*Flows).^2+a2.*(Lambda*Flows)+a3)))')+lambda2*(vec(Heads([find(Zeta'*Pi'==-1)],:)-ValveSetting).'*vec(Heads([find(Zeta'*Pi'==-1)],:)-ValveSetting)))
subject to
%% Reservoir Operational constraint
Theta*Heads == Theta*MinNodalHeads; 
%% Pump Head Gain Negative
C1M+C2M.*(Lambda*Flows)<=0; 
% Big M for Pump On/Off 
M_pump*(OnOff-1)<=(Lambda*Pi')*Heads-(C1M+C2M.*(Lambda*Flows))<=M_pump*(1-OnOff); 
%Pump Flow Constraints
Lambda*MinFlowPump.*OnOff<=Lambda*Flows<=(Lambda*MaxFlowPump).*OnOff;
% Junction Heads Constraints 
Kappa*MinNodalHeads<=Kappa*Heads<=Kappa*MaxNodalHeads 
Ppump == (APrimePump.*(Lambda*Flows)+BPrimePump.*OnOff);
%% Tank State Space 
Hdummy' == Atk*NodeTankInitialHead+Btk*(Flows'*(-Pi)'*Tau); 
% Tank Head Constraint 
(MinNodalHeads'*Tau)<=Heads'*Tau<=(MaxNodalHeads'*Tau); 
%(MinNodalHeads'*Tau)<=Heads'*Tau;
Tau'*Heads == [NodeTankInitialHead, Hdummy(Tanks,1:Time-1)];
%% Pipe Head Loss/Energy Equation  
Pi_telda*Heads == Cp+ Pi_prime*Flows ;
%% Mass Balance only for Junctions ignoring Tanks 
Pi_reduced*Flows  == -JunctionDemandProfile; 
%% Include Valve Constraints 
M_valve*(ValveStatusOnOff-1)<=(Zeta'*Pi'*Heads)-ValveHeadLoss<=M_valve*(1-ValveStatusOnOff); 
(Zeta'*MaxFlowPump).*(ValveStatusOnOff-1)<=Zeta'*Flows<=(Zeta'*MaxFlowPump).*ValveStatusOnOff;
%%0<=ValveSubStatey<=ValveStatusOnOff
%%0<=ValveSubStatez<=ValveStatusOnOff
%%ValveSubStatez+ValveSubStatey<=ones(Valves,Time);
%%M_valve*(ValveSubStatey-1)+ValveSetting<=Heads([find(Zeta'*Pi'==1)],:)<=M_valve*(1-ValveSubStatez)+ValveSetting
%0<=ValveHeadLoss<=M_valve.*ValveSubStatey;
0<=ValveHeadLoss<=M_valve*ValveStatusOnOff;
%%M_valve*(ValveSubStatey-1)+ValveSetting<=Heads([find(Zeta'*Pi'==-1)],:)<=M_valve*(1-ValveSubStatey)+ValveSetting
cvx_end

%% calculating new values for next iteration. 
Cp = (Pi_prime*Flows).*(Omega*abs(Pi_prime*Flows).^(0.852)-ones(Pipes,Time));
C1M = -h_0.*ones(Pumps,Time);
C2M = r_m.*(Lambda*Flows).^(v_m-1);
PumpPower = c_m*(h_0-r_m.*(Lambda*Flows).^2).*(Lambda*Flows) ;
PumpPrimePower = c_m*(h_0-3*r_m.*(Lambda*Flows).^2);
APrimePump = PumpPrimePower ;
BPrimePump = PumpPower-PumpPrimePower.*(Lambda*Flows); 
eff = c_m2.*(-a1.*(Lambda*Flows).^2+a2.*(Lambda*Flows)+a3);
y = (APrimePump.*(Lambda*Flows)+BPrimePump.*OnOff)./eff;
y(y<0)=0;
Eps = [Heads;Flows;OnOff;ValveStatusOnOff;ValveSubStatey;ValveSubStatez;y];
Obj = cvx_optval;
Error = norm(Eps-Eps_save)
Eps_save = Eps;
n_iter = n_iter+1;
%% saving variables from  iterates 
Flows_iter = [Flows_iter;Flows];
Eps_iter = [Eps_iter;Eps_save];
Error_iter = [Error_iter;Error];
Cp_iter = [Cp_iter;Cp];
C1M_iter = [C1M_iter;C1M];
C2M_iter = [C2M_iter;C2M];
OnOff_iter = [OnOff_iter;OnOff];
Obj_iter = [Obj_iter;Obj];
APrimePump_iter = [APrimePump_iter;APrimePump];
BPrimePump_iter = [BPrimePump_iter;BPrimePump];
PumpPower_iter = [PumpPower_iter;PumpPower];
Ppump_iter = [IntPumpPower;Ppump];
ValveStatusOnOff_iter = [ValveStatusOnOff_iter;ValveStatusOnOff];
ValveSubStatey_iter = [ValveSubStatey_iter;ValveSubStatey];
ValveSubStatez_iter = [ValveSubStatez_iter;ValveSubStatez] ;
ValveHeadLoss_iter = [ValveHeadLoss_iter;ValveHeadLoss];
end
end
cd('Decoupled_MatFiles')
save('water_decoupled_eff','Flows','Heads','OnOff','PumpPower','Ppump');
save(filename)
cd('..')