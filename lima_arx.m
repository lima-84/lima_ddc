% Author: Pedro Rodrigues de Lima
%
% LIMA_ARX  Estimates ARX parameters
%
%   Model structure:
%
%    A(z)y(t) = B(z)u(t-nk) + e(t)
%
%    A(z) = 1 + a_1*z^-1 + ... a_na*z^-na
%    B(z) = (b_0 + b_1*z^-1 + ... + b_nb*z^-nb) * z^-nk
%   
%   [G,theta] = lima_arx(y,u,nz,np,nk,Ts)
%
%   Inputs:
%       y,u: input and output data
%       nz,np,nk: number of zeros (nb-1), poles (na) and input delays (nk)
%       Ts: sampling time
%
%   Outputs:
%       G: resulting model transfer function
%       theta: parameter vector [a_1 a_2 ... a_na b_0 b_1 ... b_nb]' 

function [G,theta] = lima_arx(y,u,nz,np,nk,Ts)

    % Defining delayed input vector
    u_delay = zeros(length(u),nz+1);
    for i = 1:nz+1
        u_delay(i+nk:end,i) = u(1:end-nk-i+1);
    end
    
    % Defining delayed output vectors
    y_delay = zeros(length(y),np);
    for i = 1:np
        y_delay(i+1:end,i) = -y(1:end-i);
    end
    
    % Regressor vectors
    phi = [u_delay y_delay];
    
    % Removing extra zeros
    xtra = 1 + max(np,nz+nk);
    phi = phi(xtra:end,:);
    y = y(xtra:end,:);
    
    % Least squares
    theta = phi\y;
    
    G = tf([theta(1:nz+1)' zeros(1,np-nz)],[1 theta(nz+2:end)'],Ts);
    G.IODelay = nk;
end