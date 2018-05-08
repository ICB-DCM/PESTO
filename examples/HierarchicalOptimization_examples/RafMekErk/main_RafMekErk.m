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

compilation_RafMekErk

runEstimation_RafMekErk('hierarchical','normal')
runEstimation_RafMekErk('hierarchical','laplace')

runEstimation_RafMekErk('standard','normal')
runEstimation_RafMekErk('standard','laplace')
