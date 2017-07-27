% plotting script, plot in log10 space to match visualization in d2d

close all
%% get analytical values
if strcmp(options.llh.approach,'analytical')
    options.llh.save_analytical = 1;
    [ll] = logLikelihood_Bachmann(xi,D,options);
    load analytical_results c
end

%% simulation
options.ami.sensi = 0;
for cond = 1:numel(D)
    if numel(D(cond).t) > 1
        tsim = linspace(D(cond).t(1),D(cond).t(end));
    else
        tsim = D(cond).t;
    end
    options.ami.x0 = D(cond).init(xi,D(cond).u);
    sol_temp(cond) = simulate_Bachmann_JAKSTAT_red(tsim,xi,D(cond).u,[],options.ami);
    
end
sol = getSimulation_Bachmann_JAKSTAT_offsetscaling(xi,sol_temp,D,options);

if strcmp(options.llh.approach,'analytical')
    for cond = 1:numel(D)
        sol(cond).y = bsxfun(@times,c(:,:,:,cond),sol(cond).y);
    end
end
%% plot conditions
for cond = 1:14
    figure('name',D(cond).name);
    if numel(D(cond).t) > 1
        tsim = linspace(D(cond).t(1),D(cond).t(end));
    else
        tsim = D(cond).t;
    end
    obsplot = [];
    for iobs = 1:20
        if sum(sum((~isnan(squeeze(D(cond).my(:,iobs,:))))))>0
            obsplot = [obsplot,iobs];
        end
    end
    for iobs = 1:numel(obsplot)
        sx = round(sqrt(numel(obsplot)));
        sy = ceil(numel(obsplot)/sx);
        subplot(sx,sy,iobs);
        if obsplot(iobs) == 11 %pSTAT5B_rel
            plot(D(cond).t,(squeeze(D(cond).my(:,obsplot(iobs),:))),'k.'); hold on;
        else
            plot(D(cond).t,log10(squeeze(D(cond).my(:,obsplot(iobs),:))),'k.'); hold on;
        end
        if numel(D(cond).t) > 1
            if obsplot(iobs) == 11 %pSTAT5B_rel
                plot(tsim,((sol(cond).y(:,obsplot(iobs)))),'r'); hold on;
                xlim([D(cond).t(1),D(cond).t(end)]);
            else
                plot(tsim,log10((sol(cond).y(:,obsplot(iobs)))),'r'); hold on;
                xlim([D(cond).t(1),D(cond).t(end)]);
            end
        else
            plot(tsim,log10((sol(cond).y(:,obsplot(iobs)))),'rd'); hold on;
        end
        title(observable_names{obsplot(iobs)})
    end
    
end
%% plot dose responses in one figure each
condition_groups = {[15:19],[20:25],[26:31],[32:36]};
for group = 1:numel(condition_groups)
    figure('name',D(condition_groups{group}(1)).name);
    for cond = condition_groups{group}
        tsim = D(cond).t;
        obsplot = [];
        for iobs = 1:20
            if sum(~isnan(squeeze(D(cond).my(:,iobs,:))))>1
                obsplot = [obsplot,iobs];
            end
        end
        for iobs = 1:numel(obsplot)
            sx = round(sqrt(numel(obsplot)));
            sy = ceil(numel(obsplot)/sx);
            subplot(sx,sy,iobs);
            plot(log10(D(cond).u(end)),(squeeze(D(cond).my(:,obsplot(iobs),:))),'k.'); hold on;
            plot(log10(D(cond).u(end)),((sol(cond).y(:,obsplot(iobs)))),'rd'); hold on;
            title(observable_names{obsplot(iobs)})
            xlabel('dose');
            ylabel('conc [au]');
        end
    end
end