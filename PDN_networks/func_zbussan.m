function [Vmag,V_new1, iter,I_PQ] = func_zbussan(A_telda, y, p,q,N,bshpu)
 ybus = A_telda.'*(diag(y))*A_telda;%ybus computation
 nbus = length(ybus);
 VV = ones(N,1);
 P_new = -p; % making them injections
 Q_new = -q; % making them injections
 Y = ybus(2:nbus,2:end)+1j*diag(bshpu);
 YNS = ybus(2:nbus,1);
 Vs =1.00;
 ww = -Y\(YNS*Vs);
 I_PQ = zeros(N,1);
 V_new1 = VV; % flat start
 toler1= 1;                  % initial Tolerence.
 iter = 1;  
 maxiter = 1000;
 % iteration starting
 while (toler1 > 0.00001)&&(iter<maxiter)    % Start of while loop, tol=1.0*e-5 p.u.
 I_PQ  = conj(P_new+1j*Q_new)./conj(V_new1);
  VV = (Y\I_PQ)+ww;
 iter = iter + 1;    % Increment iteration count.
 toler1 = max(abs(abs(VV) - abs(V_new1)));     % Calculate tolerance.
 V_new1 = (VV);
 end
iter ;      % Total iterations.
V_new1 ;              % Bus Voltages in Complex form.
Vmag = abs(V_new1);