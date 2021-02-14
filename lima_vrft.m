% Author: Pedro Rodrigues de Lima
%
% LIMA_VRFT  Virtual Reference Feedback Tuning using instrumental variable
% 
%   [C,p] = lima_vrft(y1,y2,u,Td,C_bar,L,n)
%   
%   Inputs:
%       u,y1,y2: input and output data
%       Td: desired closed-loop transfer function
%       C_bar: controller's class
%       L: mismatched class VRFT filter
%       n: Td(z) input delay
%
%   Outputs:
%       C: resulting controller
%       p: parameter vector

function [C,p] = lima_vrft(y1,y2,u,Td,C_bar,L,n)

    z = zpk('z',Td.Ts);

    % Applying filter
    u_L = lsim(L,u);
    y1_L = lsim(L,y1);
    y2_L = lsim(L,y2);
    
    % Virtual error
    e1 = lsim((1-Td)/(z^(1+n)*Td),y1_L);
    e1 = e1(2+n:end);
    
    e2 = lsim((1-Td)/(z^(1+n)*Td),y2_L);
    e2 = e2(2+n:end);
    
    % Removing extra samples
    u_L = u_L(1:end-(1+n));
    
    % Least squares
    phi = lsim(C_bar,e1);
    zeta = lsim(C_bar,e2);
    p = (zeta'*phi)\(zeta'*u_L);
    
    C = minreal(zpk(p'*C_bar));
end