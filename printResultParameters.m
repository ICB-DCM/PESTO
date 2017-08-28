function printResultParameters( parameters )
    fprintf('No.\t|\texitflag\t|\tt_cpu\t|\tlog_post\t|\tpar\n');
    for j=1:length(parameters.MS.logPost)
        fprintf(strcat('%d\t|\t%d\t|\t%0.12f\t|\t%0.12f\t|\t',mat2str(parameters.MS.par(:,j)),'\n'),...
            j,...
            parameters.MS.exitflag(j),...
            parameters.MS.t_cpu(j),...
            parameters.MS.logPost(j));
    end
end

