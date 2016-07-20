close all;

options.nsimu_warmup = 1e3;
options.nsimu_run    = 1e4;

options.SCMC.n_proposals = 1;
%options.sampling_scheme = 'single-chain multi-core'; options.proposal_scheme = 'AM'; options.AM.adaption_scheme = 'difference'; options.AM.memory_length = 1;
options.sampling_scheme = 'single-chain'; options.proposal_scheme = 'MH'; options.AM.adaption_scheme = 'difference'; options.AM.memory_length = 1;
%options.sampling_scheme = 'single-chain'; options.proposal_scheme = 'MALA';
%options.sampling_scheme = 'DRAM';
options.plot_options.A.plot_type = 1;
options.plot_options.S.bins = 15;

tic
parameters = getParameterSamples(parameters,logL,options);
toc