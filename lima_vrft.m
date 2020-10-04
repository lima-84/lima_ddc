% Author: Pedro Rodrigues de Lima
%
% LIMA_VRFT  Virtual Reference Feedback Tuning using instrumental variable
%    u,y1,y2: input and output data
%    Td: desired closed-loop transfer function
%    C_bar: controller's class
%    L: mismatched class VRFT filter
%    n: Td(z) input delay

function [C,p] = lima_vrft(y1,y2,u,Td,C_bar,L,n)

    z = zpk(0,[],1,Td.Ts);

    % Applying filter
    u_L = lsim(L,u);
    y1_L = lsim(L,y1);
    y2_L = lsim(L,y2);
    
    % Virtual reference
    r1 = lsim(1/(z^(1+n)*Td),y1_L);
    r1 = r1(2+n:end);
    
    r2 = lsim(1/(z^(1+n)*Td),y2_L);
    r2 = r2(2+n:end);
    
    % Removing extra samples
    u_L = u_L(1:end-(1+n));
    y1_L = y1_L(1:end-(1+n)); 
    y2_L = y2_L(1:end-(1+n)); 
    
    % Reference tracking
    e1 = r1-y1_L;
    e2 = r2-y2_L;
    
    % Least squares
    phi = lsim(C_bar,e1);
    zeta = lsim(C_bar,e2);
    p = (zeta'*phi)\(zeta'*u_L);
    
    C = minreal(zpk(p'*C_bar));
end