function testToys(mode)
% flag:
%   0: arbdim local
%   1: arbdim global
%   2: fixeddim local
%   3: fixeddim global

    
%% prepare exercises
    
cell_exercises = createExercises(mode);

%% solve problems

nExercises = length(cell_exercises);

cell_results = cell(nExercises,1);

for jExercise = 1:nExercises
    ex = cell_exercises{jExercise};
    
    result = doExercise(ex);
    
    cell_results{jExercise} = result;
end

tabTest = Result.cell_to_table(cell_results);

save(['tabTest-',mode,'.mat'],'tabTest');
