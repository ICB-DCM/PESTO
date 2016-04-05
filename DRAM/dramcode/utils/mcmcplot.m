function mcmcplot(chain,inds,names,plottype,varargin)
%MCMCPLOT Plot mcmc chain
% mcmcplot(chain,inds,names,plottype, ...)
% plot several plots from mcmc chain
%  chain    - mcmc chain matrix
%  inds     - columns to use, default all columns
%  names    - names of the columns of chain as char array
%  plottype - 'chain', 'pairs', 'dens', 'acf', 'hist'
% example: mcmcplot(chain(2000:10:10000,:),1:4,[],'pairs')

% varargin is smo for density*

% ML, 2001
% $Revision: 1.15 $  $Date: 2006/03/22 09:28:15 $

if exist('lowess') ==2
  smooth=1; % lowess smooth with 'chain' plot
else
  smooth=0;
end
maxpoints = 500; % max n:o of points in the chain plot

if nargin  < 4 | isempty(plottype)
  plottype='chain';
end

if nargin==1 & isstruct(chain) & isfield(chain,'class') & strcmp(chain.class,'MCMC')
  % calling with results struct only
  results = chain;
  chain = GetArray('chain');
  names = results.names;
  plottype='chainpanel';
end

[nsimu,npar2]=size(chain);

if nargin < 2 | isempty(inds)
  inds = 1:npar2;
end
if nargin < 3 | isempty(names)
  for i=1:length(inds);
    names{i} = num2str(inds(i));
  end
else 
  if isstruct(names) % results struct from mcmcrun
    results=names;
    names=results.names;
  end
  names = names(inds);
end

skip=1;

switch plottype
 case('pairs')
  pairs(chain(:,inds),'panellims',names,skip,varargin{:})
 case('dens')
  for i=1:length(inds)
    [y,x]=density(chain(:,inds(i)),[],varargin{:});
    plot(x,y)
    title(sprintf('Parameter %g: %s',i,names{i}))
    drawnow
    if i~=length(inds)  % not on last time
      disp('Pause, press enter to continue ...');
      pause
    end
  end
 case('denspanel')
  np  = length(inds);
  ns1 = ceil(sqrt(np));
  ns2 = round(sqrt(np));
  for i=1:np
    h=subplot(ns1,ns2,i);
    [y,x]=density(chain(:,inds(i)),[],varargin{:});
    plot(x,y)
    set(h,'ytick',[]);
    title(sprintf('%s',names{i}))
    % add prior density
    if exist('results','var')
      ii  = inds(i);
      mus = results.prior(ii,:);
      mu  = mus(1);
      sig = mus(2);
      mi  = results.limits(ii,1);
      ma  = results.limits(ii,2);
      if ~isinf(sig)
	xdr=x(end)-x(1);
	hold on
	xp = linspace(max(mi,mu-3*sig),min([ma,mu+3*sig,x(end)+xdr]));
	yp = norpf(xp,mu,sig^2);
	yn = nordf((mi-mu)/sig);
	plot([xp(1),xp],[0,yp./(1-yn)],'--k')
	hold off
	set(gca,'ytick',[])
      end
      % force axis to min max region
      xlim = get(h,'xlim');
      if xlim(1) < mi
        xlim(1) = mi-(min(xlim(2),ma)-mi)*0.05;
      end
      if xlim(2) > ma
        xlim(2) = ma+(ma-max(xlim(1),mi))*0.05;
      end
      if ~isinf(sig)
        xlim(2) = min(xlim(2),xp(end));
        xlim(1) = max(xlim(1),xp(1));
      end
      set(h,'xlim',xlim);
    end
  end
 case('hist')
  np  = length(inds);
  ns1 = ceil(sqrt(np));
  ns2 = round(sqrt(np));
  for i=1:np
    h=subplot(ns1,ns2,i);
    histp(chain(:,inds(i)),varargin{:});
    set(h,'ytick',[]);
    title(sprintf('%s',names{i}))
  end
 case('chain')
  % max points in the plot
  if nsimu>maxpoints
    skip = floor(nsimu/maxpoints);
  end
  
  for i=1:length(inds)
    plot(1:skip:nsimu,chain(1:skip:nsimu,inds(i)),'.')
    if length(inds)>1
      title(sprintf('Parameter %g: %s',i,names{i}))
    else
      title(sprintf('%s',names{i}))
    end
    if smooth
      ys = lowess(1:skip:nsimu,chain(1:skip:nsimu,inds(i)));
      hold on
      plot(1:skip:nsimu,ys,'-k');
      hold off
    end
    drawnow
    if i~=length(inds)  % not on last time
      disp('Pause, press enter to continue ...');
      pause
    end
  end
 case('chainpanel')
  % max points in the plot
  if nsimu>maxpoints
    skip = floor(nsimu/maxpoints);
  end

  np  = length(inds);
  ns1 = ceil(sqrt(np));
  ns2 = round(sqrt(np));
  for i=1:np
    h=subplot(ns1,ns2,i);
    plot(1:skip:nsimu,chain(1:skip:nsimu,inds(i)),'.')
%    set(h,'ytick',[]);
    title(sprintf('%s',names{i}))
    if smooth
      ys = lowess(1:skip:nsimu,chain(1:skip:nsimu,inds(i)));
      hold on
      plot(1:skip:nsimu,ys,'-k');
      hold off
    end
    set(h,'xlim',[1 nsimu]);
    if i<=ns1*ns2-ns2
      set(h,'xticklabel',' ');
    end
  end
 case('acf')
  np = length(inds);
  if np>1;clf;end
  for i=1:np
    if np>1; h=subplot(np,1,i); else; h=gca; end
    plot(acf(chain(1:skip:nsimu,inds(i)),varargin{:}),'-'); 
    if i~=np; set(h,'xtick',[]); end
%    axis([1 100 -1 1]);
%    set(h,'ylim',[-1 1]);
    ylim = get(h,'ylim');
    ylim(1) = min(ylim(1),0);
    set(h,'ylim',ylim);

    title(sprintf('Parameter %g: %s',i,names{i}))
  end
 otherwise
  disp('unknown plot type');
end % switch
