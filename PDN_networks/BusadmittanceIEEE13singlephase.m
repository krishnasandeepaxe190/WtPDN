%% This script creates Bus admittance matrix for IEEE13 single phase feeder
clc;
close all
clear all; 
%cd('PDN_networks')
Network = load('IEEE13Singlephase_data');
SBase = Network.SBase; % base KVA,  5 MVA
VBase = Network.VBase; % base voltage 2.4018 KV
IBase = SBase./VBase;  %  base current  2.08 KA
ZBase = VBase.^2/SBase;  % Zbase   1.157 ohms
NodeID = Network.NodeID;
NodeID_n = NodeID(NodeID~=650); % taking the slack out
Branchdata = Network.Branchdata; 

Fromnode = Branchdata(:,2); % From node id
Tonode = Branchdata(:,1);  % To node id
Capdata = Network.Capdata; % Cap data
Capid = Capdata(:,1); % Capid
Cap_value = Capdata(:,2); % cap value
%Cap_value = [0;0];
Cap_value = Cap_value./SBase; % cap value in pu. 
N = length(NodeID_n);
p = Branchdata(:,5);
q = Branchdata(:,6);
p = p./SBase; %  p in pu
q = q./SBase; % q in pu
s_abs = sqrt(p.^2+q.^2); % s in pu
r = Branchdata(:,3);
x = Branchdata(:,4);
r = r./ZBase; % pu resistance
x = x./ZBase; % pu reactance
Tonodenew =  [1;2;3;4;5;6;7;8;9;10;11;12]; % labels as per DFS
Fromnodenew =[0;1;2;1;4;1;6;6;6;7;7;8]; % labels as per DFS
labels = Branchdata(:,1);
Capid_new = [10;12]; 
pv_ids = Network.pv_ids;
bshpu = zeros(N,1);
bshpu(Capid_new) = Cap_value;
p_new = p(p~=0);
q_new = q(q~=0); 
s_new = sqrt(p_new.^2+q_new.^2);
pf = p_new./s_new;
z = r+1j*x; % series line impedence
y = 1./z; % shunt admittance
v02=(2.4018*10^3)^2/(VBase.^2);
%% Step ii, Construct Full branch bus and reduced branch bus incidence
%matrix based on network topology and DFS algorithm
pie_n = [0 1 2 1 4 1 6 6 6 7 7 8]; 
[A_telda, a_0, A, F,V_k, R, X, Qsx,Qsxtelda, V_nsh, V_n1] = func_branchbus(pie_n,... 
N,v02,r,x,p,q,bshpu);
%% Zbus sanity check before connection matrices
%[Vmag,V_new1,iter,I_PQ] = func_zbussan(A_telda, y, p,q,N,qshpu);
ybus = A_telda.'*diag(y)*A_telda; % N+1 X N+1; 
% qshpu = zeros(N,1); 
[Vmag,V_new1,iter,I_PQ] = func_zbussan(A_telda, y, p,q,N,bshpu);
results = runpf('case13test');
bus = results.bus;
Vmag_matpower = bus(:,8);
Vmag_matpower = Vmag_matpower(2:end,:); 
cd('..');
if exist('SinglePhaseMatFiles')~=7
    mkdir 'SinglePhaseMatFiles'
end
cd('SinglePhaseMatFiles'); 
save('IEEE13SinglePhaseall');
cd('..');

