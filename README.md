# SparsePELS
Sparse Parameter Estimation for Linear Dynamical Systems

Description: This repository contains implementations of the sparse parameter estimation for linear dynamical systems.

Copyright (c) 2020 Behrad Soleimani All Rights Reserved

Contact: behrad@umd.edu

Citation: If you find these piece of codes helpful in your reserach, please cite the following paper

-Soleimani, B., P. Das, P., J. Kulasingham, J. Z. Simon and B. Babadi (2020) Granger Causal Inference from Indirect Low-Dimensional Measurements with Application to MEG Functional Connectivity Analysis, 2020 54th Annual Conference on Information Sciences and Systems.

Date: May 3, 2020

Requirements: implemented in Matlab R2019a version, but should run on most versions.

Contents: 
> main.m:       **Master script**. 

> SparsePELS.m:       **SparsePELS algorithm**. 

> EM.m:       **Expectation maximization (EM) algorithm**.

> EStep.m:       **Expectation (E-) step**.

> IRLS.m:       **Iteratively re-weighted least square (IRLS) algorithm**.

> MStep.m:       **Maximization (M-) step**.

> VARGenerator.m:       **VAR process generator**.

> redblue.m:  **Red-Blue color map**.

> SparsePELS.pdf: **Derivation and details of the algorithm**.

Instructions: Simple and easy. Download all the codes in a directory and run main.m, that will generate one example described below. To use the functions individually, please look at the function descriptions. The derivations and details are also explained in .pdf file.

Example:

In this example, we assume that there are N<sub>y</sub>=5 observations, N<sub>x</sub>=10 sources, and T=200 time samples. The underlying source dynamic is considered as a vector auto-regressive process with 2 lags, i.e. VAR(2). There are 3 active sources with unknown indeces.

In Fig.1, the estimation and ground truth versions of the VAR coefficients along with the noise covariance matrix is demonstrated. The overall relative error between estimation and ground truth is less than 8%. 

| ![](Figs/FullModel.png) | 
|:--:| 
| Fig 1. Comparison of the conventional filtering and EfficientFFBS |
