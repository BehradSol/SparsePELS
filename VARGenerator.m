function [y, e, R, C, A, Q, B] = VARGenerator(T, p, Nx, Ny, Ne)

    B = zeros(Nx,Ne);
    B(1,1) = -1;
    B(3,1) = 0.7;


    A = cell(p,1);
    for i = 1 : p
        A{i} = zeros(Nx,Nx);
    end
    A{1}(1,1) = 0.8;
    A{1}(1,3) = -0.8;
    A{1}(1,1) = -0.4;
    A{1}(4,3) = -0.7;

    A{1}(3,1) = 0.5;
    A{1}(4,4) = 0.8;

    A{2}(4,4) = -0.9;
    A{2}(3,1) = -0.7;
    A{2}(1,1) = 0.3;

    q = 0.05*ones(1,Nx);
    q(1) = 9;
    q(3) = 8;
    q(4) = 6.5;

    Q = diag(q);
    R = 0.1*eye(Ny);

    C = 2*rand(Ny,Nx) - 1;

    C = 2*(rand(Ny,Nx)>0.5) - 1;
    while (rank(C) ~= min(Ny,Nx))
        C = 2*(rand(Ny,Nx)>0.5) - 1;
    end
    
    e = zeros(Ne,T);
    % t = 2*pi/T:2*pi/T:2*pi;
    % for i = 1 : Ne
    %     e(i,:) = 10*sin(randi(Ne)*t);
    % end
    x = zeros(Nx,T+p);

    for i = 1 : T
        x(:,i+p) = x(:,i+p) + B*e(:,i) + mvnrnd(zeros(1,Nx),Q)';
        for j = 1 : p
            x(:,i+p) = x(:,i+p) + A{j}*x(:,i+p-j);
        end
    end

    y = C*x(:,p+1 : T+p) + mvnrnd(zeros(1,Ny),R,T)';


end

