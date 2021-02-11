%% This code takes the input from optimization of decoupled and verifies it with EPANET solver.
clc
clear
close all
cd('Decoupled_MatFiles')   
d_cvx = load('WaterDecoupled_Cubic_lineff_PRV_8.mat'); 
cd('..')
cd('WDN_networks')
cd('WDN_files')
d = epanet('tutorial8node_PRV_Cubic_EPANET.inp');
cd('..')
cd('..')
% d.openHydraulicAnalysis;
% d.initializeHydraulicAnalysis;
% tstep=1;
% TimeStep = d.getTimeHydraulicStep;
Results = d.getComputedHydraulicTimeSeries;
Flow_epanet = Results.Flow'; 
Head_epanet = Results.Head'; 
%Heads_epanet_dummy = Head_epanet;
%Heads_epanet_dummy(1,:) = Head_epanet(7,:);
%Heads_epanet_dummy(2:7,:) = Head_epanet(1:6,:);
%Head_epanet = Heads_epanet_dummy; 
Flows_cvx = d_cvx.Flows;
Flows_cvx(1,:) = -Flows_cvx(1,:);
Heads_cvx = d_cvx.Heads; 
OnOff_cvx = d_cvx.OnOff;
PumpPower_epanet = Results.Energy(:,10); 
PumpPower_epanet = PumpPower_epanet';
PumpPower_cvx = d_cvx.y;
eff_cvx = 0.9.*(d_cvx.eff);
eff_epanet = Results.Efficiency(:,10);
eff_epanet = eff_epanet(eff_epanet>0)
cd('Decoupled_MatFiles')
save('Decoupled_Cubic_eff_PRV_8_epanet')
cd('..')
%save('Pumppowerdcoupled','PumpPower_epanet');
%save('PumppowerDecoupled','PumpPower_epanet');