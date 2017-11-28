function sol_ret = getSimulation_Bachmann_JAKSTAT_offsetscaling(xi,sol,D,options)

% check dimension of parameters
if size(xi,1) == 1
    xi = xi';
end

for cond=1:numel(D)
    if options.llh.original % (observable = scaling*simulation + offset)
        offset = zeros(20,1);
        offset(~isnan(D(cond).offset)) = 10.^xi(D(cond).offset(~isnan(D(cond).offset)));
        scaling = ones(20,1);
        scaling(~isnan(D(cond).scaling)) = 10.^xi(D(cond).scaling(~isnan(D(cond).scaling)));
        c(1,:,1,cond) = scaling;
        sol_ret(cond).y = bsxfun(@plus,offset,bsxfun(@times,scaling,sol(cond).y'))';
        
        if options.ami.sensi
            sol(cond).sy(:,:,end+1:numel(xi)) = 0;
            sol_ret(cond).sy = zeros(size(sol(cond).sy,1),size(sol(cond).sy,2),numel(xi));
            
            doffset = zeros(1,20,numel(xi));
            doffset(1,(~isnan(D(cond).offset)),D(cond).offset(~isnan(D(cond).offset))) = ...
                diag(10.^xi(D(cond).offset(~isnan(D(cond).offset)))*log(10));
            
            dscaling = zeros(1,20,numel(xi));
            dscaling(1,(~isnan(D(cond).scaling)),D(cond).scaling(~isnan(D(cond).scaling))) = ...
                diag(10.^xi(D(cond).scaling(~isnan(D(cond).scaling)))'*log(10));
            
            for iobs = 1:size(sol(cond).y,2)
                sol_ret(cond).sy(:,iobs,:) = bsxfun(@plus,bsxfun(@plus,...
                    bsxfun(@times,sol(cond).y(:,iobs),dscaling(1,iobs,:)), ...
                    (sol(cond).sy(:,iobs,:))*scaling(iobs)),...
                    doffset(1,iobs,:));
            end
        end
    else % modified version (observable = scaling*(simulation + offset)
        offset = zeros(20,1);
            offset(~isnan(D(cond).offset)) = 10.^xi(D(cond).offset(~isnan(D(cond).offset)));
            scaling = ones(20,1);
            scaling(~isnan(D(cond).scaling)) = 10.^xi(D(cond).scaling(~isnan(D(cond).scaling)));
            c(1,:,1,cond) = scaling;
            sol_ret(cond).y = bsxfun(@times,scaling,bsxfun(@plus,offset,sol(cond).y'))';
            if options.ami.sensi
                sol(cond).sy(:,:,end+1:numel(xi)) = 0;
                sol_ret(cond).sy = zeros(size(sol(cond).sy,1),size(sol(cond).sy,2),numel(xi));
             
                doffset = zeros(1,20,numel(xi));
                doffset(1,(~isnan(D(cond).offset)),D(cond).offset(~isnan(D(cond).offset))) = ...
                    diag(10.^xi(D(cond).offset(~isnan(D(cond).offset)))*log(10));
                dscaling = zeros(1,20,numel(xi));
                dscaling(1,(~isnan(D(cond).scaling)),D(cond).scaling(~isnan(D(cond).scaling))) = ...
                    diag(10.^xi(D(cond).scaling(~isnan(D(cond).scaling)))'*log(10));
                
                for iobs = 1:size(sol(cond).y,2)
                    sol_ret(cond).sy(:,iobs,:) = bsxfun(@plus,(sol_ret(cond).sy(:,iobs,:)) + ...
                        bsxfun(@times,sol(cond).y(:,iobs),dscaling(1,iobs,:)) + (sol(cond).sy(:,iobs,:))*scaling(iobs),...
                        scaling(iobs)*doffset(1,iobs,:) + ...
                        dscaling(1,iobs,:)*offset(iobs));
                end
            end
    end
end