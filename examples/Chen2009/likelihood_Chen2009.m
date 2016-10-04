function [logL,dlogL,ddlogL] = likelihood_villab(theta,D,simulate_model,options)

f_plot = 0;

nderiv = nargout-1;

logL = 0;
if(isfield(options,'sens_ind'))
    if(nderiv>=1)
        dlogL = zeros(length(options.sens_ind),1);
    end
    if(nderiv>=2)
        ddlogL = zeros(length(options.sens_ind),length(options.sens_ind));
    end
else
    if(nderiv>=1)
        dlogL = zeros(length(theta),1);
    end
    if(nderiv>=2)
        ddlogL = zeros(length(theta),length(theta));
    end
end

for idata = 1:length(D)
    if(nderiv>=1)
        options.sensi = 1;
        options.sensi_meth = 'adjoint';
        sol = simulate_model(D(idata).t,theta,D(idata).condition,D(idata),options);
        if(sol.status<0)
            error('integration error');
        else
            logL = logL + sol.llh;
            dlogL = dlogL + sol.sllh;
        end
    else
        options.sensi = 0;
        sol = simulate_model(D(idata).t,theta,D(idata).condition,D(idata),options);
        if(sol.status<0)
            error('integration error');
        else
            logL = logL + sol.llh;
        end
    end
end

if(f_plot == 1)
    errorbar(repmat(D.t',[1,size(D.Y,2)]),D.Y,D.Sigma_Y,'LineWidth',2);
    hold on
    plot(D.t,sol.y)
    drawnow
    hold off
end

