Examples
========

# Examples {#examples}

## Overview {#overview}

PESTO comes with a number of ready-to-run examples to demonstrate its usage and functionality. The examples mostly stem from computational biology and 
comprise ordinary and partial differential equation (ODE / PDE), and Gaussian mixture models. 
More background information is provided in the respective example folder in the `main*.m` script.

The following examples are included: 
* *Conversion reaction*, the simplest example and demonstrating all PESTO features (`examples/conversion_reaction/mainConversionReaction.m`, ODE)
* *Enzymatic catalysis* (`examples/enzymatic_catalysis/mainEnzymaticCatalysis.m`, ODE)
* *Gaussian mixture* (`examples/GaussExample/mainExampleGauss.m`)
* *Hyperring* (`examples/RingExample/mainExampleRing.m`)
* *mRNA transfection* § (`examples/mRNA_transfection/mainTransfection.m`, ODE)
* *Pom1p gradient formation* § (`examples/Pom1p_gradient_formation/mainPom1.m`, PDE)
* *Jak-Stat-signaling* § (`examples/jakstat_signaling/mainJakstatSignaling.m`, ODE)
* *ERBB signaling*, largest model among the examples § (`examples/erbb_signaling/mainErbBSignaling.m`, ODE)

§ These models require the freely available AMICI toolbox to run (http://icb-dcm.github.io/AMICI/).

The following table provides an overview of which of the PESTO functions are demonstrated in the different examples: 
 
|                           | conversion reaction | enzymatic catalysis | Gaussian mixture | Hyperring | transfection | Pom1p | JakStat | ErbB |
|:-------------------------:|:-------------------:|:-------------------:|:----------------:|:---------:|:------------:|:-----:|:-------:|:----:|
| getMultiStarts()          |          X          |           X         |                  |     X     |      X       |   X   |    X    |  X   |
| plotMultiStarts()         |          X          |           X         |                  |     X     |      X       |   X   |    X    |  X   |
| getParameterConfidenceIntervals() |  X          |           X         |                  |     X     |      X       |       |         |      |
| plotConfidenceIntervals() |          X          |           X         |                  |     X     |      X       |       |         |      |
| getParameterProfiles()    |          X          |           X         |                  |     X     |      X       |       |         |      |
| plotParameterProfiles()   |          X          |           X         |          X       |     X     |      X       |   X   |         |      |
| getParameterSamples()     |          X          |           X         |          X       |     X     |      X       |       |         |      |
| plotParameterSamples()    |          X          |           X         |          X       |     X     |      X       |       |         |      |
| getPropertyMultiStarts()  |          X          |                     |                  |           |      X       |       |         |      |
| plotPropertyMultiStarts() |          X          |                     |                  |           |      X       |       |         |      |
| getPropertyConfidenceIntervals() |   X          |                     |                  |           |      X       |       |         |      |
| getPropertyProfiles()     |          X          |                     |                  |           |      X       |       |         |      |
| plotPropertyProfiles()    |          X          |                     |                  |           |      X       |       |         |      |
| getPropertySamples()      |          X          |                     |                  |           |      X       |       |         |      |
| plotPropertySamples()     |          X          |                     |                  |           |      X       |       |         |      |
| testGradient()            |                     |                     |                  |           |              |       |         |      |
| collectResults()          |                     |                     |                  |           |              |       |         |      |
