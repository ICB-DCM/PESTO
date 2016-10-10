function runTestExamples()
% runTestExamples runs some examples for the PESTO toolbox, to test if new
% implementations cause problems

run conversion_reaction/mainConversionReaction;
run enzymatic_catalysis/mainEnzymaticCatalysis;
run mRNA_transfection/main;
run jakstat_signaling/mainJakstatSignaling;

end