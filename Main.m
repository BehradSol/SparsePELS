%% This file is distributed under BSD (simplified) license
%% Author: Behrad Soleimani <behrad@umd.edu>

clc
clear 
close all


T = 500;
p = 1;

Nx = 10;
Ny = 5;
Ne = 1;


B = zeros(Nx,Ne);
B(1,1) = 1;
B(4,1) = -0.7;


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


q = 0.05*ones(1,Nx);
q(1) = 9;
q(3) = 6;
q(4) = 7;

Q = diag(q);
R = 0.01*eye(Ny);

C = 2*rand(Ny,Nx) - 1;

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


lambdaVec = [100 300 500];
[Af, Bf, Qf, LL] = CrossValidation(lambdaVec, y, e, C, R, p);



figure
errA = 100*norm(A{1}-Af{1},'F')/norm(A{1},'F')
errQ = 100*norm(Q-Qf,'F')/norm(Q,'F')
errB = 100*norm(B-Bf,'F')/norm(B,'F')

cmap = redblue(200);

subplot(3,2,1)
imagesc(A{1},[-1 1])
ylabel('$$A$$', 'Interpreter', 'LaTeX','FontSize',20,'FontWeight','bold')
colormap(cmap)     

colorbar
axis('square')

subplot(3,2,2)
imagesc(Af{1},[-1 1])
colormap(cmap)     

colorbar
ylabel('$$\tilde{A}$$', 'Interpreter', 'LaTeX','FontSize',20,'FontWeight','bold')
axis('square')


subplot(3,2,3)
imagesc(B,[-1 1])
colormap(cmap)     
colorbar
ylabel('$$B$$', 'Interpreter', 'LaTeX','FontSize',20,'FontWeight','bold')
axis('square')

subplot(3,2,4)
imagesc(Bf,[-1 1])
colormap(cmap)     
colorbar
ylabel('$$\tilde{B}$$', 'Interpreter', 'LaTeX','FontSize',20,'FontWeight','bold')
axis('square')

subplot(3,2,5)
imagesc(Q,[0 10])
colormap(gca,'summer')
colorbar
ylabel('$$Q$$', 'Interpreter', 'LaTeX','FontSize',20,'FontWeight','bold')
axis('square')

subplot(3,2,6)
imagesc(Qf,[0 10])
colormap(gca,'summer')
colorbar
ylabel('$$\tilde{Q}$$', 'Interpreter', 'LaTeX','FontSize',20,'FontWeight','bold')
axis('square')




%        
%         str = sprintf('Iteration = %d, Nx = %d, Ny = %d, p = %d', i, Nx, Ny, p);
%         title({str},'Interpreter', 'TeX','FontSize',15)
%         
%         subplot(3,1,2)
%         imagesc(diag(Qf)',[0 10])
%         colorbar
%         ylabel('diag($$\tilde{Q}$$)', 'Interpreter', 'LaTeX','FontSize',20,'FontWeight','bold')


% [m2, Cov2]  = Filtering(y, A, Q, C, R, e, B);
% 
% plot(1:T , x(1,p+1:T+p) , 'oblack' , 'LineWidth' , 1.5)
% hold on
% plot(1:T , m2(1,:) , 'blue' , 'LineWidth' , 2.5)
% plot(1:T , m(1,:) , 'red', 'LineWidth' , 1.5)
% 
% ylabel('Signal amplitude')
% xlabel('Time index')
% axis('square')
% 
% legend('Gound truth','Estimation via conventional filtering','Estimation via EFBS')
% grid on
% xlim([1 T])
% 
% 
% figure
% 
% plot(1:T , x(1,p+1:T+p) , 'oblack' , 'LineWidth' , 1.5)
% hold on
% plot(1:T , m(1,:) , 'red' , 'LineWidth' , 1.5)
% plot(1:T , m(1,:) + 3*reshape(Cov(1,1,:),1,50) , '--blue', 'LineWidth' , 1)
% plot(1:T , m(1,:) - 3*reshape(Cov(1,1,:),1,50) , '--blue', 'LineWidth' , 1)
% axis('square')
% 
% ylabel('Signal amplitude')
% xlabel('Time index')
% 
% legend('Gound truth','Estimation via EFBS','95% Quantiles')
% grid on
% xlim([1 T])
