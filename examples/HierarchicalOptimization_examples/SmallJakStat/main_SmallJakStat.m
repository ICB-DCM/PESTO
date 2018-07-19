% Main file of the JakStat signaling I example.
%
% Demonstrates the use of:
% * logLikelhoodHierarchical()
%
% Performs parameter estimation for the standard and hierarchical approach
% and Gaussian and Laplace noise.

clear all 
close all
clc

%% Compilation of simulation files using AMICI
addpath('models/')
compilation_JakStat

%% Optimization using fmincon
runEstimation_JakStat('hierarchical','normal')
runEstimation_JakStat('hierarchical','laplace')

runEstimation_JakStat('standard','normal')
runEstimation_JakStat('standard','laplace')

%% Profile calculation 
% requires runEstimation_JakStat to be called before
runProfiles_JakStat('standard','normal')
runProfiles_JakStat('standard','laplace')

runProfiles_JakStat('hierarchical','normal')
runProfiles_JakStat('hierarchical','laplace')

%% Optimization using PSwarm
runEstimation_JakStat('hierarchical','normal','pswarm')
runEstimation_JakStat('hierarchical','laplace','pswarm')

runEstimation_JakStat('standard','normal','pswarm')
runEstimation_JakStat('standard','laplace','pswarm')
