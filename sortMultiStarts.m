% sortMultiStarts.m sorts the multi-start results.
%
% USAGE:
% ======
% [parameters] = sortMultiStart(parameters)
%
% INPUTS:
% =======
% parameters ... parameters struct
%
% Outputs:
% ========
% parameters ... parameter struct with sorted MS
%
% 2014/06/12 Jan Hasenauer

function [parameters] = sortMultiStarts(parameters)

%% Sort
[~,ind] = sort(parameters.MS.logPost,1,'descend');
ind = ind([find(~isnan(parameters.MS.logPost(ind)));...
           find( isnan(parameters.MS.logPost(ind)))]);

%% Assignment of variables which are always contained in the struct
if isfield(parameters.MS,'par0')
    parameters.MS.par0 = parameters.MS.par0(:,ind);
end

if isfield(parameters.MS,'par')
    parameters.MS.par = parameters.MS.par(:,ind);
end

if isfield(parameters.MS,'logPost0')
    parameters.MS.logPost0 = parameters.MS.logPost0(ind);
end

if isfield(parameters.MS,'logPost')
    parameters.MS.logPost = parameters.MS.logPost(ind);
end

if isfield(parameters.MS,'gradient')
    parameters.MS.gradient = parameters.MS.gradient(:,ind);
end

if isfield(parameters.MS,'hessian')
    parameters.MS.hessian = parameters.MS.hessian(:,:,ind);
end

if isfield(parameters.MS,'n_objfun')
    parameters.MS.n_objfun = parameters.MS.n_objfun(ind);
end

if isfield(parameters.MS,'n_iter')
    parameters.MS.n_iter = parameters.MS.n_iter(ind);
end

if isfield(parameters.MS,'t_cpu')
    parameters.MS.t_cpu = parameters.MS.t_cpu(ind);
end

if isfield(parameters.MS,'exitflag')
    parameters.MS.exitflag = parameters.MS.exitflag(ind);
end

%% Assignment of variables which are not always contained in the struct
if isfield(parameters.MS,'par_trace')
    parameters.MS.par_trace = parameters.MS.par_trace(:,:,ind);
end

if isfield(parameters.MS,'fval_trace')
    parameters.MS.fval_trace = parameters.MS.fval_trace(:,ind);
end

