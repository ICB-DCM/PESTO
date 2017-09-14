function cell_results = testOnToyFunctions(mode,varargin)
% mode:
%     arbdim-local
%     arbdim-global
%     fixeddim-local
%     fixeddim-global

    
    %% prepare exercises

    rng(1);

    cell_exercises = createExercises(mode,varargin{:});

    %% solve problems

    nExercises = length(cell_exercises);

    cell_results = cell(nExercises,1);

    parfor jExercise = 1:nExercises
        ex = cell_exercises{jExercise};
        disp([ex.name ' ' ex.alg]);

        result = doExercise(ex);

        cell_results{jExercise} = result;
    end

    tab_results = Result.cell_to_table(cell_results);
    if(nargin > 1)
        varg = ['-',varargin{1}];
    else
        varg = '';
    end
    time = '';%['-',datestr(datetime('now'),'yyyymmddHHMMSS')];
    save(['cell_results_test-',mode,varg,time,'.mat'],'cell_results');
    save(['tab_results_test-',mode,varg,time,'.mat'],'tab_results');

end
