function cell_results = testOnToyFunctions(mode)
% mode:
%   arbdim-local
%   arbdim-global
%   fixeddim-local
%   fixeddim-global

    
%% prepare exercises
    
cell_exercises = createExercises(mode);

%% solve problems

nExercises = length(cell_exercises);

cell_results = cell(nExercises,1);

parfor jExercise = 1:nExercises
    ex = cell_exercises{jExercise};
    
    result = doExercise(ex);
    
    cell_results{jExercise} = result;
end

tab_results = Result.cell_to_table(cell_results);
save(['cell_results_test-',mode,'.mat'],'cell_results');
save(['tab_results_test-',mode,'.mat'],'tab_results');

end
