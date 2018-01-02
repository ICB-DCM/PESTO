% This file shows how to execute several PESTO examples using the
% RAmPART-sampler. The options are set in selectRun.m and the results are
% getting stored in [pwd filesep 'results'].

% Add pathes of PESTO and the examples
addpath(['..' filesep '..']);
addpath([pwd filesep 'Banana']);
addpath([pwd filesep 'Ring']);
addpath([pwd filesep 'Gauss']);
addpath([pwd filesep 'Bachmann_JakStat']);
addpath([pwd filesep 'JakStat']);
addpath([pwd filesep 'RafMekErk']);
addpath([pwd filesep 'mRNA_Transfection']);

% Add you AMICI path here (necessary for some of the examples)
if exist('amiwrap.m')~=2
    error('Please add the AMICI toolbox to your path environment.')
end

% Add results folder for SAVE files
if exist([pwd filesep 'results'],'dir')~=7
    mkdir([pwd filesep 'results'])
end

%% Run RAmPART in RING example
% The results are getting stored in '/results'. For detailed options, see
% selectRun.m
selectRun(101,[pwd filesep 'results']);

%% Run RAmPART in the RING example with another random seed
selectRun(102,[pwd filesep 'results']);

%% Run RAmPART in the GAUSS example
selectRun(301,[pwd filesep 'results']);

%% Run Parallel Tempering in the GAUSS example
selectRun(201,[pwd filesep 'results']);

%% Run RAmPART in the JakStat example 
% (has to be compiled via AMICI)
amiwrap('jakstat_pesto','jakstat_pesto_syms',[],false);
selectRun(701,[pwd filesep 'results']);

%% Run RAmPART in mRNA-transfection example with experimental data 
% (has to be compiled via AMICI)
selectRun(1301,[pwd filesep 'results']);

%% Run RB-AM in mRNA-transfection example with experimental data
selectRun(1701,[pwd filesep 'results']);
























