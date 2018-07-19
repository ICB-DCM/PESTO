% Main file of the RAF/MEK/ERK example.
%
% Demonstrates the use of:
% * logLikelhoodHierarchical()
%
% Performs parameter estimation for the standard and hierarchical approach
% and Gaussian and Laplace noise.

clear all 
close all
clc

addpath('models/')
compilation_RafMekErk

%% Optimization
runEstimation_RafMekErk('hierarchical','normal')
runEstimation_RafMekErk('hierarchical','laplace')

runEstimation_RafMekErk('standard','normal')
runEstimation_RafMekErk('standard','laplace')

%% Profile calculation 
% requires runEstimation_JakStat to be called before
runProfiles_RafMekErk('hierarchical','normal')
runProfiles_RafMekErk('hierarchical','laplace')

runProfiles_RafMekErk('standard','normal')
runProfiles_RafMekErk('standard','laplace')

