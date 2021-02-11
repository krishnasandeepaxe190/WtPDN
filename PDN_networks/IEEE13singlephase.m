%% Parameters for IEEE-13 node single phase 
clc;
clear all;
close all
%% IEEE-13 Node Network, refer to MS thesis  
NodeID = [650;632; 633; 634; 645; 646; 671; 684; 652; 680; 611;692; 675];
Frombus = [650 632 633 632 645 632 671 671 671 684 684 692].';
Tobus =   [632 633 634 645 646 671 684 692 680 611 652 675].';
% Tonodenew =  [1;2;3;4;5;6;7;8;9;10;11;12]; % labels as per DFS
% Fromnodenew =[0;1;2;1;4;1;6;6;6;7;7;8]; % labels as per DFS
r = [0.0704 0.0561 0.1269 0.1061 0.0636 0.0704 0.0636 0.0000 0.0352 0.0755...
     0.2034 0.0462].';
x = [0.2260 0.0720 0.2307 0.0846 0.0507 0.2260 0.0507 0.0001 0.1130 0.0776...
     0.0776 0.0393].';
% r   = [0.0705 0.0561 0.0017 0.1061 0.0636 0.0705 0.0636 0.2034 0.0352 0.0755...
%        0.0017 0.0462].';
% x   = [0.2261 0.072 0.0017 0.0798 0.0479 0.2261 0.0507 0.0776 0.113 0.0766...
%        0.0017 0.0393].';
% p   = [100 0 400 170 230 1255 0 170 0 170 128 843].'*10^3;
% q   = [58  0 290 125 132 718  0 151 0 80  86 462].'*10^3;
%% Accounting distributed loads and single-phase load = sum(3-phase load)/3. 
p = [33.33 0 133.33 56.67 76.67 418.33 0 56.67 0 56.67 42.67 281].'*10^3; 
q = [19.33 0 96.67  41.67 44    239.3  0 50.33 0 26.67 28.67 154].'*10^3; 
%% Selecting PV nodes (The nodes which hosts loads should have PV)
pv_ids = [0;0;1;1;1;0;0;1;0;0;1;0];
Branchdata = [Tobus Frombus r x p q];
SBase = 5000000; % 5000 KVA
VBase = 4.16*10^3/sqrt(3); 
Capnode_id = [611 675].';
Capsize = [100 200].'*10^3; %200 KVAr
Capdata = [Capnode_id Capsize]; %100 KVar.
save('IEEE13singlephase_data')
