%% This file is distributed under BSD (simplified) license
%% Author: Behrad Soleimani <behrad@umd.edu>

function [A,B,Q] = MStep(G, S, U, h, d, f, A, B, lambda, gamma, T)
    % This function implements the maximization (M-) step.
    
    % Inputs:
    % G,S,U,h,d,f = parameters computed at E-step
    % A       = VAR coefficients; A is a p * 1 cell corresponding to each  
    %           lag such that each cell is an N_x * N_x matrix
    % B       = stimuli coefficients matrix with dimension N_x * N_e
    % lambda  = regularization constant of VAR coefficients
    % gamma   = regularization constant of stimuli coefficients
    % T       = number of samples

    % Outputs:
    % A  =  updated VAR coefficients; A is a p * 1 cell corresponding to 
    %       each lag such that each cell is an N_x * N_x matrix
    % B   = updated stimuli coefficients matrix with dimension N_x * N_e
    % Q   = updated source noise covariance matrix with dimension N_x * N_x

    % ---------------------------------------------------------------------

    p = length(A);
    Nx = length(A{1});
    
    
    a = zeros(p*Nx,Nx);
    
    for i = 1 : Nx
        temp = [];
        for j = 1 : p
            temp = [temp,A{j}(i,:)];
        end
        a(:,i) = temp';
    end
    
    [a,B] = IRLS(G,U,S,h,d,lambda,gamma, a, B);
    
    
    A = cell(p,1);
    for i = 1 : Nx
        for k = 1 : p
            A{k}(i,:) = a((k-1)*Nx + 1 : k*Nx ,i)';
        end
    end

    q = (diag(a'*G*a - 2*h'*a + B*S*B' - 2*B*d + 2*B*U*a)+f)/(T-p);
    q(q<0) = 0.05;

    Q = diag(q);

end