function [A, B, Q, LL] = CrossValidation(LambdaVec, y, e, C, R, p)
    
    L = size(y);
    T = L(1,2);
    
    TA = 0.8*T;
    TB = 0.1*T;
    TC = 0.1*T;

    
    LL = zeros(1,length(LambdaVec));
    
    if (length(LambdaVec) > 1)
        for i = 1 : length(LambdaVec)
            lambda = LambdaVec(i);
            for k = 1 : 2
                fprintf('\n Training (%d) for lambda = %d ... \n', k , lambda)

                yTrain = y(:,1 : TA + (k-1)*TB);
                yTest =  y(:, TA + (k-1)*TB + 1 : TA + TB + (k-1)*TC );

                [A, B, Q, ~] = EM(yTrain, e, C, R, p, lambda);
               
                LL(i) = LL(i) + LLCalculation(yTest, A, Q, C, R, e, B, p);

            end
        end
        LL = real(LL);
%         figure(1)
%         loglog(LambdaVec,(LL/k),'-*','LineWidth',2)
%         xlabel('Regularization Coeff')
%         ylabel('Log-Likelihood')
    
        [~,idx] = max(LL/k);
    else
        idx = 1;
    end
        lambdaStar = LambdaVec(idx);

    fprintf('\n Test with  lambda = %d ... \n', lambdaStar)

    [A, B, Q, LL] =  EM(y, e, C, R, p, lambdaStar);
    fprintf('\n Full Model Done! \n')


end

function [LL] = LLCalculation(y, A, Q, C, R, e, B, p)

    [~, ~, ~, ~, ~, ~, m, AllCov] = EStep(y, A, Q, C, R, e, B);
    
    L = size(y);
    T = L(1,2);
    
    L = size(C);
    Nx = L(1,2);
%     LL = -sum(log(diag(Q)));
    
    LL = -(T-p)*sum(log(diag(Q)))/2 -(T-p)*Nx/2 - ...
         0.5*(trace(C'/R*C*AllCov) - 2*trace(y'/R*C*m)+trace(y'*y)) - T*log(2*pi) - T/2*log(det(R));
end
    
