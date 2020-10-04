% Author: Pedro Rodrigues de Lima
%
% LIMA_MINREAL  Extended version of minreal()
%   X: vector of transfer functions
%   tol: tolerance for zero-pole cancellation

function [X,X_old] = lima_minreal(X,tol)
    
    % Ensuring X is a zpk
    X = zpk(X);
    X_old = X;
    Ts = X.Ts;
    [ni, no] = size(X);
    
    % Loop through all inputs and outputs (for MIMO systems)
    for i = 1:ni
        for j = 1:no
            % Fetching zeros and poles
            zer = X.Z{i,j};
            pol = X.P{i,j};

            % Comparing to tol and applying a bitmask
            mask = abs(zer) > tol;
            zer = zer.*mask;

            % Same thing for the poles
            mask = abs(pol) > tol;
            pol = pol.*mask;
                        
            % Now minreal() does the rest
            X(i,j) = minreal(zpk(zer,pol,X.K(i,j),Ts),1e-4);
        end
    end
end
