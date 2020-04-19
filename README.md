# SparsePELS
Sparse Parameter Estimation for Linear Dynamical Systems

Description: This repository contains implementations of the sparse parameter estimation for linear dynamical systems.

Copyright (c) 2020 Behrad Soleimani All Rights Reserved

Contact: behrad@umd.edu

Citation: If you find these piece of codes helpful in your reserach, please cite the following paper

-Soleimani, B., P. Das, P., J. Kulasingham, J. Z. Simon and B. Babadi (2020) Granger Causal Inference from Indirect Low-Dimensional Measurements with Application to MEG Functional Connectivity Analysis, 2020 54th Annual Conference on Information Sciences and Systems.

Date: April 20, 2020

Requirements: implemented in Matlab R2019a version, but should run on most versions.

Contents: 
> main.m:       **Master script**. 

> EFBS.m:       **Efficient forward filtering and backward smoothing function**.

> Filtering.m:  **Conventional forward filtering and backward smoothing function**.

> EfficientFFBS.pdf: **Derivation and details of the algorithm**.

Instructions: Simple and easy. Download all the codes in a directory and run main.m, that will generate one example described below. To use the functions individually, please look at the function descriptions. The derivations and details are also explained in .pdf file.

Example:

| ![]<span style="display:block;text-align:center">(Figs/FullModel.png) </span> | 
|:--:| 
| Fig 1. Comparison of the conventional filtering and EfficientFFBS |

