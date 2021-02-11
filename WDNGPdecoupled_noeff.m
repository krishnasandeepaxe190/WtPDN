
% This code solves the WDN problem using GP approach for WDN considers PRV with variable opening and verifies
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
filename = ['WaterDecoupled_noeff_PRV_s',num2str(WDN_num),'.mat'];
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
Omega_prv = Connect_Mat.Omega_prv;
Delta = Connect_Mat.Delta; % diag(tank_Atk)
Zeta = Connect_Mat.Zeta;  % Arc-to-valve matrix
Pi_telda = Pi';
PumpLinkIndex = find(Lambda==1);
ValveLinkIndex = find(Zeta==1);
Pi_telda([PumpLinkIndex,ValveLinkIndex],:)=[];
Pipes = size(Pi_telda,1); % without valves and pumps whose flows can be changed as needed
%% Pi_reduced giving juctions only
Pi_reduced = Pi;
Pi_reduced([find(Tau==1),find(Theta==1)],:)=[];
%% Pump Parameter
h_0 = WDN.h_0;
r_m = WDN.r_m;
v_m = WDN.v_m ;
%c_m = WDN.c_m; %commented to include the pump efficiency as the variable
c_m = 0.7457/(3960*0.81);
%% efficiency coefficients
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
%MaxFlow = 1380;
MinFlowPump = zeros(Links,Time);
MaxFlowPump = repmat(MaxFlow,[Links,Time]);
%% Initial flows for the links
IntFlowPipes = 200*ones(Pipes,Time);
IntFlowPumps = 250*ones(Pumps,Time);
IntFlowValves = 300*ones(Valves,Time);
IntFlows = [IntFlowPipes;IntFlowPumps;IntFlowValves];
IntHeads = MinNodalHeads;
IntHeads([find(Tau==1)],:) = repmat(NodeTankInitialHead,[Tanks,Time]);
%% Pump Inital Binary condition for Big-M
IntOnOff  = ones(Pumps,Time);
%% Valve Inital Status condition for Big-M
IntValveStatusOnOff = ones(Valves,Time); %x in the proof ON/OFF
ValveSetting =repmat(MinNodalHeads([find(Zeta'*Pi'==-1)]),[Valves,Time])+ repmat(WDN.LinkInitialSetting([ValveLinkIndex])*2.30724939,[Valves,Time]); %changing from Psi to ft of head
ValveDownStreamNode = zeros(Nodes,Nodes);
ValveDownStreamNode(find(Zeta'*Pi'==-1),find(Zeta'*Pi'==-1))=1;
ValveUpStreamNode = zeros(Nodes,Nodes);
ValveUpStreamNode(find(Zeta'*Pi'==1),find(Zeta'*Pi'==1))=1;

%% Estimates for Pipes based on GP approximation
% Need to use Omega P X P diagonal matrix here
Omega([find(Zeta==1)],:) = []; %ignoring row
Omega(:,[find(Zeta==1)]) = [] ; % ignoring coloumn Link to Link ignoring pump Link
Omega([find(Lambda==1)],:) = []; %ignoring row
Omega(:,[find(Lambda==1)]) = [] ; % ignoring coloumn Link to Link ignoring pump Link
IntCp = IntFlowPipes.*(Omega*abs(IntFlowPipes).^(0.852)-ones(Pipes,Time));
%% Estimate for PRV based on monomial approximation
Omega_prv = Omega_prv([ValveLinkIndex:ValveLinkIndex], [ValveLinkIndex:ValveLinkIndex]);
IntValveOpen_iter = 1*ones(Valves,Time);
IntCp_prv = zeros(Valves, Time);
for i=1:length(Valves)
IntCp_prv(i,:) = IntFlowValves(i,:).*((Omega_prv(i,i)./(IntValveOpen_iter(i,:).^2)).*(IntFlowValves(i,:))-(ones(1,Time)./(IntValveOpen_iter(i,:))));
end
%% Estimates for Pump Head Gain
IntC1M = -h_0.*ones(Pumps,Time);
IntC2M = r_m.*IntFlowPumps.^(v_m-1);
%% First Order Taylor Approximation for Pump Power
IntPumpPower = c_m*(h_0-r_m.*(Lambda*IntFlows).^2).*(Lambda*IntFlows) ;
IntPumpPrimePower = c_m*(h_0-3*r_m.*(Lambda*IntFlows).^2);
IntAPrimePump = IntPumpPrimePower ;
IntBPrimePump = IntPumpPower-IntPumpPrimePower.*(Lambda*IntFlows);
%% Inital Variables and saving from Iterates
IntEps = [IntHeads;IntFlows;IntOnOff;IntValveStatusOnOff;IntValveOpen_iter];
M_pump = 4e05;
M_valve = 3e05;
M_valve_1 = 3e05;
Cp = [] ;
C1M = [];
C2M = [];
Cp_prv = [];
Flows = [];
Heads = [];
OnOff = [] ;
ValveStatusOnOff = [];
ValveOpen = [] ;
ValveSubStatusy =[];
ValveSubstatusz = [];
Eps = [] ;
Error = 1;
n_iter = 1;
Max_iter =4000;
%% Saving Iteration Values
Eps_iter = [];
Flows_iter = [];
Cp_iter = [];
Cp_prv_iter = [];
C1M_iter = [] ;
C2M_iter = [];
APrimePump_iter = [];
BPrimePump_iter = [];
Ppump_iter = [];
OnOff_iter = [] ;
ValveStatusOnOff_iter = [];
ValveSubStatusy_iter =[];
ValveSubStatusz_iter = [];
ValveOpen_iter = [];
PumpPower_iter = [];
Obj_iter = [] ;
Error_iter = [] ;
PRVHLoss_iter =[];
lambda1 = 10000;
lambda2 =0;
%lambda2 = 0.00001;
%lambda2 = 0.000001;
%lambda2 = 0.00001;
while Error >= 0.0001 && n_iter <= Max_iter
if n_iter == 1
    Cp = IntCp;
    Cp_prv = IntCp_prv;
    C1M = IntC1M;
    C2M = IntC2M;
    Eps_save = IntEps;
    APrimePump= IntAPrimePump;
    BPrimePump = IntBPrimePump ;
end
if n_iter<=30
% cvx loop inside while loop
cvx_begin
params.NumericFocus =3;
cvx_solver Gurobi
cvx_precision high
variable Flows(Links,Time)
variable Heads(Nodes,Time)
variable OnOff(Pumps,Time) binary
variable ValveStatusOnOff(Valves,Time) binary
variable ValveSubStatusy(Valves, Time) binary
variable ValveSubStatusz(Valves, Time) binary
variable ValveOpen(Valves,Time)
variable Ppump(Pumps,Time)
variable Hdummy(Tanks,Time)
variable PRVHLoss(Valves, Time)
%variable PRVOnHLoss(Valves,Time)
minimize (norm(Ppump,1))
subject to
%% Reservoir Operational constraint
Theta*Heads == Theta*MinNodalHeads;
%% Pump Head Gain Negative
C1M+C2M.*(Lambda*Flows)<=0;
% Big M for Pump On/Off
M_pump*(OnOff-1)<=(Lambda*Pi')*Heads-(C1M+C2M.*(Lambda*Flows))<=M_pump*(1-OnOff);
%Pump Flow Constraints
(Lambda*MinFlowPump).*OnOff<=Lambda*Flows<=(Lambda*MaxFlowPump).*OnOff;
% Junction Heads Constraints
Ppump == (APrimePump.*(Lambda*Flows)+BPrimePump.*OnOff);
%% Tank State Space
Hdummy' == Atk*NodeTankInitialHead+Btk*(Flows'*(-Pi)'*Tau);
% Tank Head Constraint
(MinNodalHeads'*Tau)<=Heads'*Tau<=(MaxNodalHeads'*Tau);
Tau'*Heads == [NodeTankInitialHead, Hdummy(Tanks,1:Time-1)];
%% Pipe Head Loss/Energy Equation
Pi_telda*Heads == Cp+ Pi_prime*Flows ;
%% Mass Balance only for Junctions ignoring Tanks
Pi_reduced*Flows  == -JunctionDemandProfile;
%% Junction Head Constraints
%Kappa*MinNodalHeads<=Kappa*Heads<=Kappa*MaxNodalHeads;
%% Include Valve Constraints
%(Zeta'*Pi')*Heads==(Cp_prv.*ValveOpen+(Zeta'*Flows))
%M_valve*(ValveSubStatusz-1)<=(Zeta'*Pi')*Heads-PRVOnHLoss<=M_valve*(1-ValveSubStatusz);

M_valve*(ValveStatusOnOff-1)<=(Zeta'*Pi')*Heads-PRVHLoss<=M_valve*(1-ValveStatusOnOff);
PRVHLoss==(Cp_prv.*ValveOpen+(Zeta'*Flows));
%0<=ValveSubStatusy<=ValveStatusOnOff;
%0<=ValveSubStatusz<=ValveStatusOnOff;
ValveSubStatusy+ValveSubStatusz==ValveStatusOnOff;
M_valve_1*(ValveSubStatusy-1)+ValveSetting<=Heads(5,:)<=M_valve_1*(1-ValveSubStatusz)+ValveSetting;
M_valve_1*(ValveSubStatusy-1)+ValveSetting<=Heads(7,:)<=M_valve_1*(1-ValveSubStatusy)+ValveSetting;
0*(ValveStatusOnOff)<=Zeta'*Flows<=(Zeta'*MaxFlowPump).*ValveStatusOnOff;
ValveSubStatusz+0.00001.*ValveSubStatusy+0.00000001.*(1-ValveStatusOnOff)<=ValveOpen<=ValveSubStatusz+ValveSubStatusy+0.00000001.*(1-ValveStatusOnOff);
%0.00001<=ValveStatusOnOff<=ValveOpen<=1.*ValveStatusOnOff;
%0.0001.*ValveSubStatusy<=ValveOpen<=1.*ValveStatusOnOff
%0<=PRVHLoss<=M_valve.*(1-ValveSubStatusz)
%0<=PRVOnHLoss<=M_valve.*(1-ValveSubStatusz)
cvx_end

%% calculating new values for next iteration.
% FlowValves = Zeta'*Flows;
% Cp_prv = zeros(Valves,Time);
% for i=1:length(Valves)
% Cp_prv(i,:) = FlowValves(i,:).*(((Omega_prv(i,i)./(ValveOpen(i,:).^2)).*(FlowValves(i,:)))-(ones(1,Time)./(ValveOpen(i,:))));
% end
% Cp = (Pi_prime*Flows).*(Omega*abs(Pi_prime*Flows).^(0.852)-ones(Pipes,Time));
% C1M = -h_0.*ones(Pumps,Time);
% C2M = r_m.*(Lambda*Flows).^(v_m-1);
% PumpPower = c_m*(h_0-r_m.*(Lambda*Flows).^2).*(Lambda*Flows) ;
% PumpPrimePower = c_m*(h_0-3*r_m.*(Lambda*Flows).^2);
% APrimePump = PumpPrimePower ;
% BPrimePump = PumpPower-PumpPrimePower.*(Lambda*Flows);
% Eps = [Heads;Flows;OnOff;ValveStatusOnOff;ValveOpen];
% Obj = cvx_optval;
% Error = norm(Eps-Eps_save)
% Eps_save = Eps;
% n_iter = n_iter+1;
% %% saving variables from  iterates
% Flows_iter = [Flows_iter;Flows];
% Eps_iter = [Eps_iter;Eps_save];
% Error_iter = [Error_iter;Error];
% Cp_iter = [Cp_iter;Cp];
% Cp_prv_iter = [Cp_prv_iter;Cp_prv];
% C1M_iter = [C1M_iter;C1M];
% C2M_iter = [C2M_iter;C2M];
% OnOff_iter = [OnOff_iter;OnOff];
% Obj_iter = [Obj_iter;Obj];
% APrimePump_iter = [APrimePump_iter;APrimePump];
% BPrimePump_iter = [BPrimePump_iter;BPrimePump];
% PumpPower_iter = [PumpPower_iter;PumpPower];
% Ppump_iter = [IntPumpPower;Ppump];
% ValveStatusOnOff_iter = [ValveStatusOnOff_iter;ValveStatusOnOff];
% ValveSubStatusy_iter = [ValveSubStatusy_iter;ValveSubStatusy];
% ValveSubStatusz_iter = [ValveSubStatusz_iter;ValveSubStatusz];
% ValveOpen_iter = [ValveOpen_iter;ValveOpen];
% PRVHLoss_iter = [PRVHLoss_iter;PRVHLoss];
else
cvx_begin
params.NumericFocus =3;
cvx_solver Gurobi
cvx_precision high
variable Flows(Links,Time)
variable Heads(Nodes,Time)
variable OnOff(Pumps,Time) binary
variable ValveStatusOnOff(Valves,Time) binary
variable ValveSubStatusy(Valves, Time) binary
variable ValveSubStatusz(Valves, Time) binary
variable ValveOpen(Valves,Time)
variable Ppump(Pumps,Time)
variable Hdummy(Tanks,Time)
variable PRVHLoss(Valves,Time)
minimize (norm(Ppump,1))
subject to
%% Reservoir Operational constraint
Theta*Heads == Theta*MinNodalHeads;
%% Pump Head Gain Negative
C1M+C2M.*(Lambda*Flows)<=0;
% Big M for Pump On/Off
M_pump*(OnOff-1)<=(Lambda*Pi')*Heads-(C1M+C2M.*(Lambda*Flows))<=M_pump*(1-OnOff);
%Pump Flow Constraints
(Lambda*MinFlowPump).*OnOff<=Lambda*Flows<=(Lambda*MaxFlowPump).*OnOff;
% Junction Heads Constraints
Ppump == (APrimePump.*(Lambda*Flows)+BPrimePump.*OnOff);
%% Tank State Space
Hdummy' == Atk*NodeTankInitialHead+Btk*(Flows'*(-Pi)'*Tau);
% Tank Head Constraint
(MinNodalHeads'*Tau)<=Heads'*Tau<=(MaxNodalHeads'*Tau);
Tau'*Heads == [NodeTankInitialHead, Hdummy(Tanks,1:Time-1)];
%% Pipe Head Loss/Energy Equation
Pi_telda*Heads == Cp+ Pi_prime*Flows ;
%% Mass Balance only for Junctions ignoring Tanks
Pi_reduced*Flows  == -JunctionDemandProfile;
%% Junction Head Constraints
Kappa*MinNodalHeads<=Kappa*Heads<=Kappa*MaxNodalHeads;
%% Include Valve Constraints
%(Zeta'*Pi')*Heads==(Cp_prv.*ValveOpen+(Zeta'*Flows));
%M_valve*(ValveStatusOnOff-1)<=(Zeta'*Pi'*Heads)-ValveOpen<=M_valve*(1-ValveStatusOnOff);
%M_valve.*(ValveStatusOnOff-1)<=(Zeta'*Pi')*Heads-(Cp_prv.*ValveOpen+(Zeta'*Flows))<=M_valve.*(1-ValveStatusOnOff);
%0<=ValveSubStatusy<=ValveStatusOnOff;
%0<=ValveSubStatusz<=ValveStatusOnOff;
M_valve*(ValveStatusOnOff-1)<=(Zeta'*Pi')*Heads-PRVHLoss<=M_valve*(1-ValveStatusOnOff);
PRVHLoss==(Cp_prv.*ValveOpen+(Zeta'*Flows));
ValveSubStatusy+ValveSubStatusz==ValveStatusOnOff;
M_valve_1*(ValveSubStatusy-1)+ValveSetting<=Heads(5,:)<=M_valve_1*(1-ValveSubStatusz)+ValveSetting;
M_valve_1*(ValveSubStatusy-1)+ValveSetting<=Heads(7,:)<=M_valve_1*(1-ValveSubStatusy)+ValveSetting;
0*(ValveStatusOnOff)<=Zeta'*Flows<=(Zeta'*MaxFlowPump).*ValveStatusOnOff;
%ValveSubStatusz+0.00001.*ValveSubStatusy<=ValveOpen<=ValveStatusOnOff;
ValveSubStatusz+0.00001.*ValveSubStatusy+0.00000001.*(1-ValveStatusOnOff)<=ValveOpen<=ValveSubStatusz+ValveSubStatusy+0.00000001.*(1-ValveStatusOnOff);
%0.0001.*ValveSubStatusy<=ValveOpen<=1.*ValveStatusOnOff;
%0<=PRVHLoss<=M_valve.*(1-ValveSubStatusz)
%0.001.*ValveStatusOnOff<=ValveOpen<=1.*ValveStatusOnOff
cvx_end
end
%% calculating new values for next iteration.
FlowValves = Zeta'*Flows;
Cp_prv = zeros(Valves,Time);
for i=1:length(Valves)
Cp_prv(i,:) = FlowValves(i,:).*(((Omega_prv(i,i)./(ValveOpen(i,:).^2)).*(FlowValves(i,:)))-(ones(1,Time)./(ValveOpen(i,:))));
end
Cp = (Pi_prime*Flows).*(Omega*abs(Pi_prime*Flows).^(0.852)-ones(Pipes,Time));
C1M = -h_0.*ones(Pumps,Time);
C2M = r_m.*(Lambda*Flows).^(v_m-1);
PumpPower = c_m*(h_0-r_m.*(Lambda*Flows).^2).*(Lambda*Flows) ;
PumpPrimePower = c_m*(h_0-3*r_m.*(Lambda*Flows).^2);
APrimePump = PumpPrimePower ;
BPrimePump = PumpPower-PumpPrimePower.*(Lambda*Flows);
Eps = [Heads;Flows;OnOff;ValveStatusOnOff;ValveOpen];
Obj = cvx_optval;
Error = norm(Eps-Eps_save)
Eps_save = Eps;
n_iter = n_iter+1;
%% saving variables from  iterates
Flows_iter = [Flows_iter;Flows];
Eps_iter = [Eps_iter;Eps_save];
Error_iter = [Error_iter;Error];
Cp_iter = [Cp_iter;Cp];
Cp_prv_iter = [Cp_prv_iter;Cp_prv];
C1M_iter = [C1M_iter;C1M];
C2M_iter = [C2M_iter;C2M];
OnOff_iter = [OnOff_iter;OnOff];
Obj_iter = [Obj_iter;Obj];
APrimePump_iter = [APrimePump_iter;APrimePump];
BPrimePump_iter = [BPrimePump_iter;BPrimePump];
PumpPower_iter = [PumpPower_iter;PumpPower];
Ppump_iter = [IntPumpPower;Ppump];
ValveStatusOnOff_iter = [ValveStatusOnOff_iter;ValveStatusOnOff];
ValveSubStatusy_iter = [ValveSubStatusy_iter;ValveSubStatusy];
ValveSubStatusz_iter = [ValveSubStatusz_iter;ValveSubStatusz];
ValveOpen_iter = [ValveOpen_iter;ValveOpen];
PRVHLoss_iter = [PRVHLoss_iter;PRVHLoss];
end
 cd('Decoupled_MatFiles')
%save('water_decoupled_noeff_PRV','Flows','Heads','OnOff','PumpPower','Ppump');
save(filename)
cd('..')
