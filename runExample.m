
clear
close all
clc

%% Setup data settings (dS)
% settings_ds105_01;
settings_HCP;
% settings_br_tum;

subject = str2double(dS.subjStr);

%% Preprocessing

% % Head motion correction
% preprocessing(dS);

%% Run 3D

% % Run SVB
% VBMethod = 'SVB3D';
% runSVB(dS,VBMethod);

% % Run MCMC
% MCMCMethod = 'MCMC3D';
% samplingMethod = 'PCG';
% runMCMC(dS,MCMCMethod,VBMethod,samplingMethod);

% % Compute marginal and joint PPMs
% computePPMs(dS.outputPath,subject,VBMethod);
% computePPMs(dS.outputPath,subject,MCMCMethod);

%% Run 2D

% Run SVB
VBMethod = 'SVB2D';
% a = tic;
% runSVB(dS,VBMethod);
% b = toc(a);

% Run MCMC
MCMCMethod = 'MCMC2D';
samplingMethod = 'PCG';
% c = tic;
% runMCMC_Per(dS,MCMCMethod,VBMethod,samplingMethod);
% d = toc(c);

% Requires Tools for NIfTI and ANALYZE image Matlab-package
computePPMs(dS.outputPath,dS.subjStr,VBMethod);
computePPMs(dS.outputPath,dS.subjStr,MCMCMethod);
