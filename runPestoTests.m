function runPestoTests()
    % runPestoTests Run a set of PESTO unit tests
    import matlab.unittest.TestSuite
    run(TestSuite.fromFolder('tests'))
end

