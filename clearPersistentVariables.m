function clearPersistentVariables()
% clearPersistentVariables() clears variables which are used by PESTO in
% its algorithms to improve performance. Can be called manually by the user
% between running different simulations.

    % Clear variables in getFiniteDifferences()
    clear getFiniteDifferences;
    % alternatively:
    % getFiniteDifferences([], [], 4);
end