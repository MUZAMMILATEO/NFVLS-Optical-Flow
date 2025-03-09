function [u_capN, v_capN, thetaOld] = FSFOV_2021_first(fx,fy,ft,fxx,fyy,fxy,fxt,fyt,u,v,u_cap,v_cap,thetaOld)
%THIS FUNCTION EVALUATES U_cap and V_cap, USING U AND V, AND EVALUATES THE
%SYSTEM OF EQUATIONS BASED ON DATA PENALTY.

theta = 1/thetaOld;
lambda =10000; %10000;
delta = 0.001;    %0.4;
delta2 = delta^2;
kappa = lambda*delta;

%%%%%%%%%% Evaluate T term %%%%%%%%%%
T = sqrt(1 + ((ft + fx.*u_cap + fy.*v_cap)./delta).^2);
T = 1./T;
%%%%%%%%%% Evaluated T term %%%%%%%%%%

C = (kappa/delta2).*T;

D = theta.*C.*(fxx + fyy) + theta^2;
D = 1./D;

b1 = u.*theta - C.*fxt;
b2 = v.*theta - C.*fyt;

a = C.*fyy + theta;
c = C.*fxx + theta;
b = -C.*fxy;

u_capN = D.*(a.*b1 + b.*b2);
v_capN = D.*(b.*b1 + c.*b2);





end

