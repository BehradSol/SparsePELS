%% This file is distributed under BSD (simplified) license
%% Author: Behrad Soleimani <behrad@umd.edu>

function [y, e, R, C, A, Q, B] = VARGenerator(T, p, Nx, Ny, Ne)
    % This function generates a VAR process.
    
    % Inputs:
    % T  = number of samples
    % p  = VAR model order (number of lags)
    % Nx = number of sources
    % Ny = number of observations
    % Ne = number of (external) stimuli

    % Outputs:
    % y   = observations matrix with dimension N_y * T
    % e   = (external) stimuli vector with dimension N_e * 1
    % R   = observation noise covariance matrix with dimension N_y * N_y
    % C   = linear mapping between observations and sources with dimension 
    %       N_y * N_x
    % A   = VAR coefficients; A is a p * 1 cell corresponding to each lag 
    %       such that each cell is an N_x * N_x matrix
    % Q   = source noise covariance matrix with dimension N_x * N_x
    % B   = stimuli coefficients matrix with dimension N_x * N_e

    % ---------------------------------------------------------------------


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

