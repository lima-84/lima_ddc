% Author: Pedro Rodrigues de Lima
% 
% LIMA_VDFT_FLEX   Virtual Disturbance Feedback Tuning using the flexible criterion
%    u,y: input and output data
%    F: vector of transfer functions (note that Qd(z) = n'F(z))
%    Qf: fixed part of Qd(z) (note that F(z)=Qf(z)*F_bar(z))
%    C0: initial controller (that provides p0)
%    C_bar: controller's class
%    K: K(z) filter (considering J_VD multiplied by Qd(z))
%    N: number of iterations

function[C,p,Q,n] = lima_vdft_flex(y,u,F,Qf,C0,C_bar,K,N)

    % Fetching num and den as filter() is faster than lsim()
    [num_F,den_F] = tfdata(F,'v');
    [num_Qf,~] = tfdata(Qf,'v'); % den_Qf = 1
    [num_C0,den_C0] = tfdata(C0,'v');
    [num_C_bar,den_C_bar] = tfdata(C_bar,'v');
    [num_K,den_K] = tfdata(K,'v');
    
    F_length = length(F);
    C_length = length(C_bar);
    
    % If length(F) > 1, tfdata returns cell arrays
    if F_length == 1
       den_Q = den_F;
    else
       den_Q = den_F{1}; 
    end

    % Initializing parameter vectors (n,p)
    n = zeros(F_length,1);
    p = zeros(C_length,1);
    
    % Initializing C(z,p) as C(z,p0)
    num_C = num_C0;
    den_C = den_C0;
    
    for i = 1:N
        % w(p_{i-1},t) = K(z)(u(t) + C(z,p_{i-1})y(t))
        w = filter(num_K,den_K,u+filter(num_C,den_C,y));
        
        % Regressor vectors for n_i
        if F_length == 1
            phi = filter(num_F,den_F,w);
        else
            for j = 1:F_length
                phi(:,j) = filter(num_F{j},den_F{j},w);
            end
        end
        Ky = filter(num_K,den_K,y);
        
        % Least squares for n_i
        n(:,i) = (phi'*phi)\(phi'*Ky);
    
        % Multiplying by Qf
        num_Q = conv(n(:,i)',num_Qf);
        num_Q = [zeros(1,length(den_Q)-length(num_Q)) num_Q];

        % Multiplying Q(z,n_i) and K(z)
        num_QK = conv(num_Q,num_K);
        den_QK = conv(den_Q,den_K);

        % v(n_i,t) = Qd(z,n_i)K(z)y(t)
        v = filter(num_Q,den_Q,Ky);

        % Regressor vectors for p_i
        if C_length == 1
           tau = filter(num_C_bar,den_C_bar,v); 
        else
           for j = 1:C_length
              tau(:,j) = filter(num_C_bar{j},den_C_bar{j},v); 
           end
        end
        tau2 = filter(num_QK,den_QK,u) - Ky;

        % Least squares for p
        p(:,i) = -(tau'*tau)\(tau'*tau2);

        % Setting controller for next iteration
        num_C = p(:,i)';
        % In case the classes of C0 and C_bar are different
        if i == 1
           den_C = den_C_bar{1}; 
        end        
    end    
    % Setting final vectors
    n = n(:,end);
    p = p(:,end);
    % Computing transfer functions
    Q = n'*F;
    C = p'*C_bar;
end