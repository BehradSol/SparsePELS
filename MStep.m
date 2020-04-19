%% This file is distributed under BSD (simplified) license
%% Author: Behrad Soleimani <behrad@umd.edu>

function [A,B,Q] = MStep(G, S, U, h, d, f, A, B,lambda,gamma,T)
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