function [theta, Jvr] = lima_vrft_lpv(y1,y2,u,Td,L,Q,n,m,f)

    z = zpk('z',Td.Ts);
   
    % Erro virtual e_bar(t) (invariância no tempo de Td(z))
    aux = (1-Td)/(z*Td);
    e_bar1 = lsim(aux,y1);
    e_bar1 = e_bar1(2:end);

    e_bar2 = lsim(aux,y2);
    e_bar2 = e_bar2(2:end);

    % Removendo amostras extras
    u = u(1:end-1);
    y1 = y1(1:end-1);
    y2 = y2(1:end-1);
    f = f(1:end-1,:);

    % Entrada filtrada u_L(t)
    u_L = lsim(L,u);

    % Vetores a_j(t)
    R = minreal(L/Q);
    a1 = zeros(length(y1),m);
    a2 = zeros(length(y2),m);
    for j = 1:m
        a1(:,j) = lsim(R,f(:,j).*e_bar1);
        a2(:,j) = lsim(R,f(:,j).*e_bar2);
    end

    N = length(y1);
    % Matrizes regressoras psi(t) e zeta(t)
    psi = zeros(n,m);
    Psi = zeros(N,n*m);
    zeta = zeros(n,m);
    Zeta = zeros(N,n*m);
    for i = 1:N
        psi = [a1(i,:); psi(1:end-1,:)];
        zeta = [a2(i,:); zeta(1:end-1,:)];
        % Concatenando as colunas de psi e zeta
        Psi(i,:) = reshape(psi,1,[]);
        Zeta(i,:) = reshape(zeta,1,[]);
    end

    % Identificação de theta
    theta = (Zeta'*Psi)\(Zeta'*u_L);
    Jvr = mean((u_L - Psi*theta).^2);
    
    theta = reshape(theta,n,m);    
    
end