function [A, B, Q, LL] = EM(y, e, C, R, p, lambda, max_iterations, tol)
    
    if nargin < 7
        max_iterations = 1e3;
    end
    
    if nargin < 8
        tol = 1e-4;    
    end
    
    L = size(y);
    T = L(1,2);
    
    L = size(C);
    Nx = L(1,2);

    L = size(e);
    Ne = L(1,1);

    A = cell(p,1);
    for i = 1 : p
        A{i} = zeros(Nx,Nx);
    end
    
    Q = diag(ones(1,Nx));
    
    B = zeros(Nx, Ne);

    alternating = 1;
    gamma = lambda;

    LLvec = zeros(1,max_iterations);
    LLvec(1) = -0.5*trace(y'*y) -T*log(2*pi) - T/2*log(det(R));
    
    for i = 2 : max_iterations
        fprintf('%d ,',i-1)

        [G, S, U, h, d, f, m, AllCov] = EStep(y, A, Q, C, R, e, B);
        for k = 1 : alternating
            [A,B,Q] = MStep(G, S, U, h, d, f, A, B,lambda,gamma,T);
        end
        LLvec(i) = -(T-p)*sum(log(diag(Q)))/2 -(T-p)*Nx/2 - ...
                    0.5*(trace(C'/R*C*AllCov) - 2*trace(y'/R*C*m)+trace(y'*y)) - T*log(2*pi) - T/2*log(det(R));
        if (abs( (LLvec(i)-LLvec(i-1))/LLvec(i-1) ) < tol)
            LLvec(i:end) = LLvec(i);
            break;
        end

    end
    LL = LLvec(end);
end

