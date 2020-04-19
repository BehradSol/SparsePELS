function  [a,B] = IRLS(G,U,S,h,d,lambda,gamma, a, B)
    L = size(h);
    Nx = L(1,2);
    p = L(1,1) / Nx;

    NumOfIteration = 20;
    
    delta = 0.1;
   
    for m = 1 : NumOfIteration   
        for i = 1 : Nx
            temp = [G + lambda*diag(1./(sqrt(a(:,i).^2 + delta^2))) , U' ; U , S + gamma*diag(1./(sqrt(B(i,:).^2 + delta^2)))]\[h(:,i);d(:,i)];
            a(:,i) = temp(1:p*Nx);
            B(i,:) = temp(p*Nx+1 : end)';
        end      
    end 


end

