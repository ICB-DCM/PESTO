% Main file of the JakStat signaling model I  example.
%
% Demonstrates the use of:
% * logLikelhoodHierarchical()
%
% Performs parameter estimation for the standard and hierarchical approach
% and Gaussian and Laplace noise.

clear all 
close all
clc

compilation_JakStat
%%
runEstimation_JakStat('hierarchical','normal')
runEstimation_JakStat('hierarchical','laplace')

runEstimation_JakStat('standard','normal')
runEstimation_JakStat('standard','laplace')
%%
runEstimation_JakStat('hierarchical','normal','pswarm')
runEstimation_JakStat('hierarchical','laplace','pswarm')

runEstimation_JakStat('standard','normal','pswarm')
runEstimation_JakStat('standard','laplace','pswarm')
