%% This file is distributed under BSD (simplified) license
%% Author: Behrad Soleimani <behrad@umd.edu>

function [G, S, U, h, d, f, m, AllCov] = EStep(y, Ao, Qo, Co, R, e, Bo)
    
    smallvalue = 1e-15;
    p = length(Ao);
    Nx = length(Ao{1});
    NNx = p*Nx;
    Ny = length(R);
    L = size(y);
    T = L(1,2);
        
    if nargin < 6
        Ne = 1;
        Bo = zeros(Nx , Ne);
        e = zeros(Ne , T);
    end
    
    m0 = zeros(NNx,1);
    
    L = size(Bo);
    Ne = L(1,2);
    
    A = [];
    temp = zeros((p-1)*Nx,p*Nx);
 
    for j = 1 : p
        A = [A, Ao{j}];
        for i = 1 : p-1
            if (i==j)
                temp((i-1)*Nx+1:(i)*Nx, (j-1)*Nx+1:(j)*Nx) = eye(Nx);
            end
        end
    end
    A = [A;temp];
    
   
    Q = Qo;
    for i = 1 : p-1
        Q = blkdiag(Q, smallvalue*eye(Nx));
    end
    
    C = zeros(Ny,NNx);
    C(:,1:Nx) = Co;
    
    B = zeros(NNx,Ne);
    B(1:Nx , :) = Bo;
    
    Sig = idare(A',A'*C',Q,C*Q*C' + R,Q*C',[]);   %Sigma_{k|k}

    if (isempty(Sig))
        disp('Augmented A matrix is not stable!')
    end
    
    SigP = A*Sig*A' + Q;    %Sigma_{k+1|k}

    K = SigP(:,1:Nx)*Co'/( Co*SigP(1:Nx,1:Nx)*Co'+ R);   %Kalman gain
    S = Sig*A'/SigP;    %Smoothing gain
    
 
    mF = zeros(NNx,T);

    mLast = m0;
  
    mTemp = zeros(NNx,1);
    
    for i = 1 : 1 : T
        
        mTemp(Nx+1:end , :) = mLast(1: NNx - Nx , :);
        mTemp(1:Nx , :) = A(1:Nx , :)*mLast + B(1:Nx , :) * e(:,i);
          
        mF(:,i) = mTemp + K*(y(:,i) - Co*mTemp(1:Nx , :));

        mLast = mF(:,i);
             
    end
    
    
    mB = zeros(NNx,T);
    mB(:,end) = mF(:,end);
    
    Cov = zeros(p*Nx , p*Nx , T);
    Cov(: , : , T) = Sig;
    for i = T-1 : -1 : 1         
        mB(:,i) = mF(:,i) + S*(mB(:,i+1) - A*mF(:,i)- B * e(:,i));
        
        Cov(: , : , i) = Sig + S*(Cov(: , : , i+1) - SigP)*S';
    end
 
       
    m = mB( 1 : Nx,:); %Mean vectors with dimension Nx * T
    
    Gain = S;
    
    
    AllCov = sum(Cov(1 :Nx ,1 : Nx , p+1 : T),3) + m*m';

    G = zeros(p*Nx, p*Nx);
    h = zeros(p*Nx,Nx);
    f = zeros(Nx,1);

    [P,D] = eig(Gain);

    for i = 1 : p*Nx
        for j = i : p*Nx
            t1 = ceil(i/Nx);
            t2 = i - (t1-1)*Nx;

            q1 = ceil(j/Nx);
            q2 = j - (q1-1)*Nx;

            tempCov = sum(Cov( : , : , p - t1 + 1 : T - t1),3);
            tempGain = real(P*diag(diag(D .^ (q1 - t1)))*P');
            
            G(i,j) = tempGain(q2 , : ) * tempCov(: , t2) + m(t2 , p - t1 + 1 : T - t1)*m(q2 , p - q1 + 1 : T - q1)';

            G(j,i) = G(i,j);

            
        end
    end


    for i = 1 : Nx
       for k = 1 : p*Nx
            t1 = ceil(k/Nx);
            t2 = k - (t1-1)*Nx;

            tempCov = sum(Cov(:, : , p + 1 : T ),3);
            tempGain = real(P*diag(diag(D .^ (t1)))*P');

            h(k,i) = tempGain(t2 , : ) * tempCov(: , i) + m(i , p + 1 : T)*m(t2 , p + 1 - t1 : T - t1)';

       end
           f(i,1) = tempCov(i , i) + m(i , p + 1 : T)*m(i ,  p + 1 : T)';

    end

    S = zeros(Ne,Ne);
    d = zeros(Ne,Nx);
    U = zeros(Ne,p*Nx);

    S = e(:,p+1:T)*e(:,p+1:T)';
    d = e(:,p+1:T)*m(:,p+1:T)';
    
    X = [];
    for i = 1 : Nx
        for k = 1 : p
            X = [X , m(i, p+1-k : T-1-k+1)'];
        end
    end
    U = e(:,p+1:T)*X;
    

end


