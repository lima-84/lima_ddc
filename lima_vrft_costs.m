% Author: Pedro Rodrigues de Lima
%
% LIMA_VRFT_COSTS   Computes VRFT J_y and J_VR cost functions
%    u,y: input and output data
%    Td: VRFT's desired closed-loop transfer function
%    C: VRFT's estimated controller C(z,p)
%    L: mismatched class VRFT filter
%    y_hat : estimated output using C(z,p)
%    yd : desired y(t)

function [Jy,Jvr] = lima_vrft_costs(y,u,Td,C,L,y_hat,yd)

    z = zpk(0,[],1,Td.Ts);
    
    % Calculating e_bar
    aux = zpk((1-Td)/Td);
    rd = length(aux.P{1}) - length(aux.Z{1});
    aux = aux*(z^(rd));
    e_bar = lsim(aux,y);
    e_bar = e_bar(1-rd:end);
    
    % Removing extra samples from u
    u = u(1:end+rd);
    
    % Calculating C(z)e_bar(t)
    Ce = lsim(C,e_bar);
    
    % Computing J_VR = E_bar[u(t)-C(z,p)e_bar(t)]^2
    Jaux = lsim(L,u-Ce);    
    Jvr = Jaux'*Jaux/length(u);
    
    % Computing J_y = E_bar[(Td(z)-T(z,p))r(t)]^2
    Jy = (yd-y_hat)'*(yd-y_hat)/length(y_hat);
end

