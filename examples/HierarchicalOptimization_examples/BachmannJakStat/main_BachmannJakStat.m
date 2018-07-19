% Main file of the Bachmann JAK-STAT signaling example.
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
compilation_BachmannJakStat

runEstimation_BachmannJakStat('hierarchical','normal')
runEstimation_BachmannJakStat('hierarchical','laplace')

runEstimation_BachmannJakStat('standard','normal')
runEstimation_BachmannJakStat('standard','laplace')
