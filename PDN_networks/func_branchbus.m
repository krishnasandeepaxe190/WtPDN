%% Function to generate branch bus incidence matrix based on topology
function[A_telda, a_0, A, F,V_k, R, X, Qsx,Qsxtelda, V_nsh, V_n1] = func_branchbus(pie_n, N, v02, r,x,p, q, bshpu)
A_telda = zeros(N,N+1); % Full branch bus incidence matrix, see proof
for i = 1:N
    A_telda(i,pie_n(i)+1) =1;
    A_telda(i,i+1) = -1;
end
a_0 = A_telda(:,1);
A = A_telda(:,2:end); % reduced branch bus incidence matrix N x N
F = -inv(A); % negative inverse of the reduced branch incidence matrix, kekatos formulation.
%% Voltage equation parameters R & X, R = 2*F*diag(r)*F', X = 2*F*diag(x)*F', V = R*p+X*q+V_k
V_k = F*a_0*v02; % constant part in the voltage equation, also F*a_0 = 1_N
R =  2*F*diag(r)*F.';% matrix related to active power injection
X =  2*F*diag(x)*F.';% matrix related to reactive power injection.
Qsx = X*diag(bshpu);% reactive power flow using shunt cap for voltage term with nominal voltage 1.0 Pu2
V_nsh = sqrt(((eye(N)-Qsx))\(R*(-p)+X*(-q)+V_k));% kekatos voltage with shunt
q = q-bshpu;
V_n1 = sqrt(R*(-p)+X*(-q)+V_k); % checking the nodal nomial voltages using net injections.
Qsxtelda = inv(eye(N)-Qsx);