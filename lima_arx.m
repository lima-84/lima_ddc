% Author: Pedro Rodrigues de Lima
%
% LIMA_ARX   Estimates ARX parameters using instrumental variable
%    u,y1,y2: input and output data
%    nz,np,nk: number of zeros (nb-1), poles (na) and input delays (nk)
%    Ts: sampling time

function [G,theta] = lima_arx(y1,y2,u,nz,np,nk,Ts)

    % Defining input delayed vector
    u_delay = zeros(length(u),nz+1);
    for i = 1:nz+1
        u_delay(i+nk:end,i) = u(1:end-nk-i+1);
    end
    
    % Defining output delayed vectors
    y1_delay = zeros(length(y1),np);
    y2_delay = zeros(length(y2),np);
    for i = 1:np
        y1_delay(i+1:end,i) = y1(1:end-i);
        y2_delay(i+1:end,i) = y2(1:end-i);
    end
    
    % Least squares
    phi = [u_delay y1_delay];
    zeta = [u_delay y2_delay];
    theta = (zeta'*phi)\(zeta'*y1);
	
    G = zpk(tf(theta(1:nz+1)',[1 -theta(nz+2:end)'],Ts));
end