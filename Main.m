%% This file is distributed under BSD (simplified) license
%% Author: Behrad Soleimani <behrad@umd.edu>

clc
clear 
close all


T = 200;
p = 2;

Nx = 5;
Ny = 5;
Ne = 1;


[y, e, R, C, A, Q, B] = VARGenerator(T, p, Nx, Ny, Ne);

lambdaVec = [200];
[Af, Bf, Qf, LL] = SparsePELS(lambdaVec, y, e, C, R, p);

cmap = redblue(200);

subplot(4,2,1)
imagesc(A{1},[-1 1])
ylabel('$$A_1$$', 'Interpreter', 'LaTeX','FontSize',20,'FontWeight','bold')
colormap(cmap)     

colorbar
axis('square')

subplot(4,2,2)
imagesc(Af{1},[-1 1])
colormap(cmap)     

colorbar
ylabel('$$\widehat{A}_1$$', 'Interpreter', 'LaTeX','FontSize',20,'FontWeight','bold')
axis('square')

subplot(4,2,3)
imagesc(A{2},[-1 1])
ylabel('$$A_2$$', 'Interpreter', 'LaTeX','FontSize',20,'FontWeight','bold')
colormap(cmap)     

colorbar
axis('square')

subplot(4,2,4)
imagesc(Af{2},[-1 1])
colormap(cmap)     

colorbar
ylabel('$$\widehat{A}_2$$', 'Interpreter', 'LaTeX','FontSize',20,'FontWeight','bold')
axis('square')


subplot(4,2,5)
imagesc([B,zeros(Nx,Nx-Ne)],[-1 1])
colormap(cmap)     
colorbar
ylabel('$$B$$', 'Interpreter', 'LaTeX','FontSize',20,'FontWeight','bold')
axis('square')

subplot(4,2,6)
imagesc([Bf,zeros(Nx,Nx-Ne)],[-1 1])
colormap(cmap)     
colorbar
ylabel('$$\widehat{B}$$', 'Interpreter', 'LaTeX','FontSize',20,'FontWeight','bold')
axis('square')

subplot(4,2,7)
imagesc(Q,[0 10])
colormap(gca,'summer')
colorbar
ylabel('$$Q$$', 'Interpreter', 'LaTeX','FontSize',20,'FontWeight','bold')
axis('square')

subplot(4,2,8)
imagesc(Qf,[0 10])
colormap(gca,'summer')
colorbar
ylabel('$$\widehat{Q}$$', 'Interpreter', 'LaTeX','FontSize',20,'FontWeight','bold')
axis('square')


