function printResult( x, fval, exitflag, output )
% print result of optimization algorithm

    disp('Best guess:');
    disp(['fval: ', num2str(fval)]);
    disp(['x:    ', mat2str(x(:)')]);
    disp('Algorithm data:');
    if (isfield(output,'algorithm'))
        disp(['algorithm:            ', output.algorithm]);
    end
    disp(['exitflag:             ', num2str(exitflag)]);
    if (isfield(output,'iterations'))
        disp(['iterations:           ', num2str(output.iterations)]);
    end
    if (isfield(output,'funcCount'))
        disp(['function evaluations: ', num2str(output.funcCount)]);
    end
    if (isfield(output,'message'))
        disp(['message:              ', output.message]);
    end
    if (isfield(output,'t_cpu'))
        disp(['t_cpu:                ', num2str(output.t_cpu)]);
    end

end

