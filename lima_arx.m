% Author: Pedro Rodrigues de Lima
%
% LIMA_ARX  Estimates ARX parameters using instrumental variable
%
%   Model structure:
%
%    A(z)y(t) = B(z)u(t-nk) + e(t)
%
%    A(z) = 1 + a_1*z^-1 + ... a_na*z^-na
%    B(z) = (b_0 + b_1*z^-1 + ... + b_nb*z^-nb) * z^-nk
%   
%   [G,theta] = lima_arx(y1,y2,u,nz,np,nk,Ts)
%
%   Inputs:
%       u,y1,y2: input and output data
%       nz,np,nk: number of zeros (nb-1), poles (na) and input delays (nk)
%       Ts: sampling time
%
%   Outputs:
%       G: resulting model transfer function
%       theta: parameter vector [a_1 a_2 ... a_na b_0 b_1 ... b_nb]' 

function [G,theta] = lima_arx(y1,y2,u,nz,np,nk,Ts)

    % Defining delayed input vector
    u_delay = zeros(length(u),nz+1);
    for i = 1:nz+1
        u_delay(i+nk:end,i) = u(1:end-nk-i+1);
    end
    
    % Defining delayed output vectors
    y1_delay = zeros(length(y1),np);
    y2_delay = zeros(length(y2),np);
    for i = 1:np
        y1_delay(i+1:end,i) = -y1(1:end-i);
        y2_delay(i+1:end,i) = -y2(1:end-i);
    end
    
    % Regressor vectors
    phi = [u_delay y1_delay];
    zeta = [u_delay y2_delay];
    
    % Removing extra zeros
    xtra = 1 + max(np,nz+nk);
    phi = phi(xtra:end,:);
    zeta = zeta(xtra:end,:);
    y1 = y1(xtra:end,:);
    
    theta = (zeta'*phi)\(zeta'*y1);

    G = zpk(tf(theta(1:nz+1)',[1 theta(nz+2:end)'],Ts));
    G.IODelay = nk;
end