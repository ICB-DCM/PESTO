function sol_ret = getSimulation_Bachmann_JAKSTAT_offsetscaling(xi,sol,D,options)
if size(xi,1) == 1
    xi = xi';
end
for cond=1:numel(D)
    switch options.llh.approach
        case 'standard'
            offset = zeros(20,1);
            offset(~isnan(D(cond).offset)) = 10.^xi(D(cond).offset(~isnan(D(cond).offset)));
            scaling = ones(20,1);
            scaling(~isnan(D(cond).scaling)) = 10.^xi(D(cond).scaling(~isnan(D(cond).scaling)));
            c(1,:,1,cond) = scaling;
            if isfield(options.llh,'original') && options.llh.original % y = offset + scaling*x
                % lines below removed and 1 substracted from the data!
                %                 if cond == 3
                %                     offset([12:17]) = 1;
                %                 end
                sol_ret(cond).y = bsxfun(@plus,offset,bsxfun(@times,scaling,sol(cond).y'))';
            else % y = scaling*(offset + x)
                % lines below removed and 1 substracted from the data!
                % if cond == 3
                %     offset([12:17]) = 1;%10.^(-xi(D(cond).scaling(~isnan(D(cond).scaling))));
                % end
                sol_ret(cond).y = bsxfun(@times,scaling,bsxfun(@plus,offset,sol(cond).y'))';
            end
            if options.ami.sensi
                sol(cond).sy(:,:,end+1:numel(xi)) = 0;
                sol_ret(cond).sy = zeros(size(sol(cond).sy,1),size(sol(cond).sy,2),numel(xi));
                %                 if cond == 3
                %                     doffset = zeros(1,20,numel(xi));
                %                     doffset(1,(~isnan(D(cond).scaling)),D(cond).scaling(~isnan(D(cond).scaling))) = ...
                %                         diag(-10.^(-xi(D(cond).scaling(~isnan(D(cond).scaling))))'*log(10));
                %                 else
                doffset = zeros(1,20,numel(xi));
                doffset(1,(~isnan(D(cond).offset)),D(cond).offset(~isnan(D(cond).offset))) = ...
                    diag(10.^xi(D(cond).offset(~isnan(D(cond).offset)))*log(10));
                %                 end
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
            % if options.llh.save_analytical
            %     save c_standard c
            % end
        case 'analytical'
            offset = zeros(20,1);
            offset(~isnan(D(cond).offset)) = 10.^xi(D(cond).offset(~isnan(D(cond).offset)));
            %             if cond == 3
            %                 offset([12:17]) = 1;
            %             end
            sol_ret(cond).y = (offset+sol(cond).y')';
            if options.ami.sensi
                sol(cond).sy(:,:,end+1:numel(xi)) = 0;
                sol_ret(cond).sy = sol(cond).sy;%zeros(size(sol(cond).sy,1),size(sol(cond).sy,2),numel(xi));
                doffset = zeros(1,20,numel(xi));
                doffset(1,(~isnan(D(cond).offset)),D(cond).offset(~isnan(D(cond).offset))) = ...
                    diag(10.^xi(D(cond).offset(~isnan(D(cond).offset)))*log(10));
                for iobs = 1:size(sol(cond).y,2)
                    sol_ret(cond).sy(:,iobs,:) = (sol_ret(cond).sy(:,iobs,:)) + doffset(1,iobs,:);
                end
            end
    end
    
end
