clc
close all
clear
%% Plotting cubic efficiency curve
Q_bep = 1200; 
Q_cutoff = 2*Q_bep;
A_prime = [Q_bep.^3 Q_bep.^2 Q_bep;
           3*Q_bep.^2 2*Q_bep 1;
           Q_cutoff.^2 Q_cutoff 1];
b_coeff = [0.9;0;0];
%b_coeff = [0.9;0.0000;0];

coeff = inv(A_prime)*b_coeff ; 
Flows = 0:1:2400;
Flows = Flows(:,2:2400);
ww_cubic = 0.9.*(coeff(1).*(Flows).^3+coeff(2).*(Flows).^2+coeff(3).*(Flows)+0.00001);
curve_points = [1 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400];
% ww_cubic_points = zeros(1,length(curve_points));
% for i = 1:length(curve_points)
%     ww_cubic_points(i) = ww_cubic(curve_points(i))*100;
% end
%% Including speeds
Speeds = 0:1;
Flows_eff = 0:2400;
%Speeds = Speeds(:,2:2400);
[X, Y] = meshgrid(Flows_eff,Speeds);
ww_cubic_speedflows = 0.9.*(coeff(1).*(X./Y).^3+coeff(2).*(X./Y).^2+coeff(3).*(X./Y)+0.00001);
surf(X,Y,ww_cubic_speedflows)
%% plotting quadratic efficiency curve
a1 = (0.9-0.00001)./Q_bep.^2;
a2 = (2*(0.9-0.00001))./Q_bep;
a3 = 0.00001; 
ww_Quadratic = 0.9.*(-a1.*(Flows).^2+a2.*(Flows)+a3);

hold on
plot(ww_cubic)
plot(ww_Quadratic)
legend('cubic','Quadratic')
