function [phi_u, phi_v] = FSFOV_2021_third(phi_u,phi_v,u_plus,u_minus,v_plus,v_minus,u_cap,v_cap,ep,tempmat,thetaOld,lambdaFrac,lev)
%THIS FUNCTION EVALUATES NEW ITERATES FOR PHI_U_NEW AND PHI_V_NEW, USING
%THE PREVIOUSLY ESTIMATED ITERATES FOR PHI_U, PHI_V, U_PLUS, U_MINUS,
%V_PLUS, V_MINUS, U_CAP AND V_CAP. (MU IS A COEFFICIENT WHICH MEASURES HOW
%TIGHTLY THE LEVEL CURVES ENCLOSE HOMOGENEOUS REGIONS.)

h = 1;
mu = 1000;  %1000;
delT = 10^lev;%1000;
P = 1/(2*thetaOld);
[m,n] = size(phi_u);

%%%%%%%%%%%%%%%%%% EVALUATING COEFFICIENTS FOR PHI_U_NEW %%%%%%%%%%%%%%%%%%
temp = zeros(m+1,n);
temp(1:m,:) = phi_u;
phi1C1 = (temp(2:m+1,:) - phi_u).^2;
temp1 = zeros(m,n+2);
temp2 = temp1;
temp1(:,1:n) = phi_u;
temp2(:,3:n+2) = phi_u;
phi1C2 = ((temp1(:,2:n+1) - temp2(:,2:n+1))./2).^2;
C1 = (phi1C1 + phi1C2  + 0.001).^(-0.5);

temp = zeros(m+1,n);
temp(2:m+1,:) = phi_u;
phi1C1 = (phi_u - temp(1:m,:)).^2;
temp = temp(1:m,:);
temp1 = zeros(m,n+2);
temp2 = temp1;
temp1(:,1:n) = temp;
temp2(:,3:n+2) = temp;
phi1C2 = (temp1(:,2:n+1) - temp2(:,2:n+1)./2).^2;
C2 = (phi1C1 + phi1C2 + 0.001).^(-0.5);

temp1 = zeros(m+2,n);
temp2 = temp1;
temp1(1:m,:) = phi_u;
temp2(3:m+2,:) = phi_u;
phi1C1 = ((temp1(2:m+1,:) - temp2(2:m+1,:))./2).^2;
temp = zeros(m,n+1);
temp(:,1:n) = phi_u;
phi1C2 = (temp(:,2:n+1) - phi_u).^2;
C3 = (phi1C1 + phi1C2 + 0.001).^(-0.5);

tempRow = zeros(m,n+1);
tempRow(:,1:n) = phi_u;
tempRow = tempRow(:,2:n+1);
temp1 = zeros(m+2,n);
temp2 = temp1;
temp1(1:m,:) = tempRow;
temp2(3:m+2,:) = tempRow;
phi1C1 = ((temp1(2:m+1,:) - temp2(2:m+1,:))./2).^2;
temp = zeros(m,n+1);
temp(:,2:n+1) = phi_u;
phi1C2 = (phi_u - temp(:,1:n)).^2;
C4 = (phi1C1 + phi1C2 + 0.001).^(-0.5);

C = C1 + C2 + C3 + C4;
%%%%%%%%%%%%%%%%%% EVALUATING COEFFICIENTS FOR PHI_U_NEW %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% EVALUATING COEFFICIENTS FOR PHI_V_NEW %%%%%%%%%%%%%%%%%%
temp = zeros(m+1,n);
temp(1:m,:) = phi_v;
phi1C1 = (temp(2:m+1,:) - phi_v).^2;
temp1 = zeros(m,n+2);
temp2 = temp1;
temp1(:,1:n) = phi_v;
temp2(:,3:n+2) = phi_v;
phi1C2 = ((temp1(:,2:n+1) - temp2(:,2:n+1))./2).^2;
D1 = (phi1C1 + phi1C2 + 0.001).^(-0.5);

temp = zeros(m+1,n);
temp(2:m+1,:) = phi_v;
phi1C1 = (phi_v - temp(1:m,:)).^2;
temp = temp(1:m,:);
temp1 = zeros(m,n+2);
temp2 = temp1;
temp1(:,1:n) = temp;
temp2(:,3:n+2) = temp;
phi1C2 = (temp1(:,2:n+1) - temp2(:,2:n+1)./2).^2;
D2 = (phi1C1 + phi1C2 + 0.001).^(-0.5);

temp1 = zeros(m+2,n);
temp2 = temp1;
temp1(1:m,:) = phi_v;
temp2(3:m+2,:) = phi_v;
phi1C1 = ((temp1(2:m+1,:) - temp2(2:m+1,:))./2).^2;
temp = zeros(m,n+1);
temp(:,1:n) = phi_v;
phi1C2 = (temp(:,2:n+1) - phi_v).^2;
D3 = (phi1C1 + phi1C2 + 0.001).^(-0.5);

tempRow = zeros(m,n+1);
tempRow(:,1:n) = phi_v;
tempRow = tempRow(:,2:n+1);
temp1 = zeros(m+2,n);
temp2 = temp1;
temp1(1:m,:) = tempRow;
temp2(3:m+2,:) = tempRow;
phi1C1 = ((temp1(2:m+1,:) - temp2(2:m+1,:))./2).^2;
temp = zeros(m,n+1);
temp(:,2:n+1) = phi_v;
phi1C2 = (phi_v - temp(:,1:n)).^2;
D4 = (phi1C1 + phi1C2 + 0.001).^(-0.5);

D = D1 + D2 + D3 + D4;
%%%%%%%%%%%%%%%%%% EVALUATING COEFFICIENTS FOR PHI_V_NEW %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% CALCULATING DIRAC APPROXIMATION OF %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% PHI_U AND PHI_V %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirac_phi_u = 0.3183*ep*(ep^2 + phi_u.^2).^(-1);
dirac_phi_v = 0.3183*ep*(ep^2 + phi_v.^2).^(-1);
%%%%%%%%%%%%%%%%%% CALCULATING DIRAC APPROXIMATION OF %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% PHI_U AND PHI_V %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% ESTIMATING COEFFICIENTS NU_U AND NU_V %%%%%%%%%%%%%%%%%%
K_u = (delT/h^2).*dirac_phi_u;
K_v = (delT/h^2).*dirac_phi_v;
K = delT*mu/h^2;
nu_u = K*dirac_phi_u;
nu_v = K*dirac_phi_v;
%%%%%%%%%%%%%%%%%% ESTIMATED COEFFICIENTS NU_U AND NU_V %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% ESTIMATING C AND D COEFFICIENTS %%%%%%%%%%%%%%%%%%%%%%%%
C = 1 + nu_u.*C;
D = 1 + nu_v.*D;
%%%%%%%%%%%%%%%%%% ESTIMATED C AND D COEFFICIENTS %%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% ESTIMATING FRACTIONAL ORDER DERIVATIVE KERNELS %%%%%%%%%
[M,N] = size(tempmat);
aX = (M+1)/2;
bX = (N+1)/2;
wX = tempmat(aX:end,bX);
wY = tempmat(aX,bX:end);
%%%%%%%%%%%%%%%%%% ESTIMATED FRACTIONAL ORDER DERIVATIVE KERNELS %%%%%%%%%%

%%%%%%%%%%%%%%%%%% ESTIMATING FRACTIONAL ORDER DERIVATIVES %%%%%%%%%%%%%%%%
frac_u_plusX = lambdaFrac*(conv2(u_plus,wX,'same')).^2;
frac_u_plusY = lambdaFrac*(conv2(u_plus,wY,'same')).^2;
frac_u_minusX = lambdaFrac*(conv2(u_minus,wX,'same')).^2;
frac_u_minusY = lambdaFrac*(conv2(u_minus,wY,'same')).^2;

frac_v_plusX = lambdaFrac*(conv2(v_plus,wX,'same')).^2;
frac_v_plusY = lambdaFrac*(conv2(v_plus,wY,'same')).^2;
frac_v_minusX = lambdaFrac*(conv2(v_minus,wX,'same')).^2;
frac_v_minusY = lambdaFrac*(conv2(v_minus,wY,'same')).^2;
%%%%%%%%%%%%%%%%%% ESTIMATED FRACTIONAL ORDER DERIVATIVES %%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% ESTIMATING THE SQUARED DIFFERENCE TERMS %%%%%%%%%%%%%%%%
u_plus_cap = -P.*(u_plus - u_cap).^2;
u_minus_cap = P.*(u_minus - u_cap).^2;

v_plus_cap = -P.*(v_plus - v_cap).^2;
v_minus_cap = P.*(v_minus - v_cap).^2;
%%%%%%%%%%%%%%%%%% ESTIMATED THE SQUARED DIFFERENCE TERMS %%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% CONSTRUCT SHIFTED LEVEL SURFACES %%%%%%%%%%%%%%%%%%%%%%%
tempC1 = zeros(m+1,n);
tempC1(1:m,:) = phi_u;
tempC1 = tempC1(2:m+1,:);

tempC2 = zeros(m+1,n);
tempC2(2:m+1,:) = phi_u;
tempC2 = tempC2(1:m,:);

tempC3 = zeros(m,n+1);
tempC3(:,1:n) = phi_u;
tempC3 = tempC3(:,2:n+1);

tempC4 = zeros(m,n+1);
tempC4(:,2:n+1) = phi_u;
tempC4 = tempC4(:,1:n);

tempD1 = zeros(m+1,n);
tempD1(1:m,:) = phi_v;
tempD1 = tempD1(2:m+1,:);

tempD2 = zeros(m+1,n);
tempD2(2:m+1,:) = phi_v;
tempD2 = tempD2(1:m,:);

tempD3 = zeros(m,n+1);
tempD3(:,1:n) = phi_v;
tempD3 = tempD3(:,2:n+1);

tempD4 = zeros(m,n+1);
tempD4(:,2:n+1) = phi_v;
tempD4 = tempD4(:,1:n);
%%%%%%%%%%%%%%%%%% CONSTRUCTED SHIFTED LEVEL SURFACES %%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% ESTIMATING NEW ITERATIONS FOR PHI_U AND PHI_V %%%%%%%%%%
phi_u = C.^(-1).*(phi_u + nu_u.*(C1.*tempC1+ C2.*tempC2 + C3.*tempC3 + C4.*tempC4) + K_u.*(u_plus_cap + u_minus_cap - frac_u_plusX - frac_u_plusY + frac_u_minusX + frac_u_minusY));
phi_v = D.^(-1).*(phi_v + nu_v.*(D1.*tempD1+ D2.*tempD2 + D3.*tempD3 + D4.*tempD4) + K_v.*(v_plus_cap + v_minus_cap - frac_v_plusX - frac_v_plusY + frac_v_minusX + frac_v_minusY));

%phi_u = C.^(-1).*(phi_u + nu_u.*(C3.*tempC3 + C4.*tempC4) + K_u.*(u_plus_cap + u_minus_cap - frac_u_plusX - frac_u_plusY + frac_u_minusX + frac_u_minusY));
%phi_v = D.^(-1).*(phi_v + nu_v.*(D3.*tempD3 + D4.*tempD4) + K_v.*(v_plus_cap + v_minus_cap - frac_v_plusX - frac_v_plusY + frac_v_minusX + frac_v_minusY));
%%%%%%%%%%%%%%%%%% ESTIMATED NEW ITERATIONS FOR PHI_U AND PHI_V %%%%%%%%%%%
% figure(500),imshow(phi_u);
end

