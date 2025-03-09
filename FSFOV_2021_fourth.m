function [u,v] = FSFOV_2021_fourth(u_plus,u_minus,v_plus,v_minus,phi_u,phi_v,ep)
%THIS FUNCTION EVALUATES NEW ITERATIONS FOR U AND V USING THE NEWLY
%OBTAINED ITERATIONS FOR U_PLUS, U_MINUS, V_PLUS, V_MINUS, PHI_U AND PHI_V.

%%%%%%%%%%%% EVALUATING THE APPROXIMATION OF HEAVISIDE'S %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% UNIT STEP FUNCTION USING PHI_U AND PHI_V %%%%%%%%%%%%%%%%%%%%%
heav_u = 0.5*(1 + 0.6366*atan(phi_u./ep));
heav_v = 0.5*(1 + 0.6366*atan(phi_v./ep));
%%%%%%%%%%%% EVALUATED THE APPROXIMATION OF HEAVISIDE'S %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% UNIT STEP FUNCTION USING PHI_U AND PHI_V %%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% ESTIMATING NEW ITERATIONS FOR U AND V %%%%%%%%%%%%%%%%%%%%%%%%
u = u_plus.*heav_u + u_minus.*(1-heav_u);
v = v_plus.*heav_v + v_minus.*(1-heav_v);
%%%%%%%%%%%% ESTIMATED NEW ITERATIONS FOR U AND V %%%%%%%%%%%%%%%%%%%%%%%%%

end

