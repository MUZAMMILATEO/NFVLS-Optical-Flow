function [u_plus, u_minus, v_plus, v_minus, ep] = FSFOV_2021_second(u_cap,v_cap,u,v,thetaOld,phiU,phiV,tempmat,wndcff,ep,lambdaFrac)
%THIS FUNCTION USES U_cap AND V_cap EVALUATED FROM FSFOV_2021_FIRST()
%FUNCTION, SECOND ORDER FRACTIONAL ORDER DERIVATIVES OF
%U_plus(U_minus,V_plus,V_minus) AND THE APPROXIMATIONS TO HEAVISIDE'S UNIT
%STEP FUNCTION ALONG WITH LEVEL CURVES AND EVALUATES NEW ITERATIONS FOR
%U_plus(U_minus,V_plus,V_minus).

%%%%%%%%%%%%%%%% EVALUATE HEAVISIDE'S UNIT STEP FUNCTION %%%%%%%%%%%%%%%%%%
heavU = 0.5.*(1 + 0.6366.*atan(phiU/ep));
heavV = 0.5.*(1 + 0.6366.*atan(phiV/ep));
%%%%%%%%%%%%%%%% EVALUATED HEAVISIDE'S UNIT STEP FUNCTION %%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% EVALUATE SEGMENTED U_cap AND SEGMENTED V_cap FUNCTION %%%%
u_cap_P = u_cap.*heavU;
u_cap_N = u_cap.*(1 - heavU);
v_cap_P = v_cap.*heavV;
v_cap_N = v_cap.*(1 - heavV);
%%%%%%%%%%%%%%%% EVALUATED SEGMENTED U_cap AND SEGMENTED V_cap FUNCTION %%%

%%%%%%%%%%%%%%%% EVALUATE SEGMENTED U_plus, U_minus AND SEGMENTED %%%%%%%%%
%%%%%%%%%%%%%%%% U_plus, U_minus FROM THE PREVIOUS ITERATION STEP %%%%%%%%%
u_plus = u.*heavU;
u_minus = u.*(1 - heavU);
v_plus = v.*heavV;
v_minus = v.*(1 - heavV);
%%%%%%%%%%%%%%%% EVALUATE SEGMENTED U_plus, U_minus AND SEGMENTED %%%%%%%%%
%%%%%%%%%%%%%%%% U_plus, U_minus FROM THE PREVIOUS ITERATION STEP %%%%%%%%%


%%%%%%%%%%%%%%%% EVALUATE SECOND ORDER FRACTIONAL DERIVATIVES FOR %%%%%%%%%
%%%%%%%%%%%%%%%% U_plus, U_minus, V_plus, V_minus with respect to x and y %
u_plus = conv2(u_plus,tempmat,'same');
u_minus = conv2(u_minus,tempmat,'same');
v_plus = conv2(v_plus,tempmat,'same');
v_minus = conv2(v_minus,tempmat,'same');
%%%%%%%%%%%%%%%% EVALUATED SECOND ORDER FRACTIONAL DERIVATIVES FOR %%%%%%%%
%%%%%%%%%%%%%%%% U_plus, U_minus, V_plus, V_minus with respect to x and y %

%%%%%%%%%%%%%%%% THIS COEFFICIENT CAN BE PLACED OUTSISE OF THIS %%%%%%%%%%%
%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = 2*thetaOld;
R = 1 + lambdaFrac*theta*wndcff;
R = 1/R;
%%%%%%%%%%%%%%%% THIS COEFFICIENT CAN BE PLACED OUTSISE OF THIS %%%%%%%%%%%
%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% ESTIMATING NEW U_plus, U_minus and V_plus, V_minus %%%%%%%
u_plus = R.*(u_cap_P + lambdaFrac*theta.*u_plus);
u_minus = R.*(u_cap_N + lambdaFrac*theta.*u_minus);
v_plus = R.*(v_cap_P + lambdaFrac*theta.*v_plus);
v_minus = R.*(v_cap_N + lambdaFrac*theta.*v_minus);
%%%%%%%%%%%%%%%% ESTIMATING NEW U_plus, U_minus and V_plus, V_minus %%%%%%%

end

