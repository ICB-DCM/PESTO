Examples
========

# Examples {#examples}

## Overview {#overview}

PESTO comes with a number of ready-to-run examples to demonstrate its usage and functionality. All examples come from the life science field. The are examples based on ordinary and partial differential equation (ODE / PDE), and Gaussian mixture models. More background information is provided in the respective example folder in the `main.m` script.

The following examples are included: 
* Conversion reaction (`examples/conversion_reaction/runEstimation__CR.m`, ODE)
* Enzymatic catalysis (`examples/enzymatic_catalysis/`, ODE)
* Gaussian mixture (`examples/Gaussian_mixture/`)
* mRNA transfection § (`examples/mRNA_transfection/`, ODE)
* Pom1p gradient formation § (`examples/Pom1p_gradient_formation/`, PDE)
* TODO: JAK/STAT signaling (`...`, ODE)

§ These models require the freely available AMICI toolbox to run (http://icb-dcm.github.io/AMICI/).

The following table provides an overview of which of the PESTO functions are demonstrated in the different examples: 
 
|                           |conversion_reaction  | enzymatic_catalysis | Gaussian_mixture | mRNA transfection | Pom1p_gradient_formation |
|:-------------------------:|:-------------------:|:-------------------:|:----------------:|:-----------------:|:------------------------:|
| getMultiStarts()          |          X          |           X         |          X       |         X         |            X             |
| getParameterConfidenceIntervals() |  X          |           X         |          X       |         X         |            X             |
| getParameterProfiles()    |          X          |           X         |          X       |         X         |            X             |
| getPropertyConfidenceIntervals() |   X          |           X         |          X       |         X         |            X             |
| getPropertyMultiStarts()  |          X          |           X         |          X       |         X         |            X             |
| getPropertyProfiles()     |          X          |           X         |          X       |         X         |            X             |
| getPropertySamples()      |          X          |           X         |          X       |         X         |            X             |
| plotMultiStarts()         |          X          |           X         |          X       |         X         |            X             |
| plotParameterProfiles()   |          X          |           X         |          X       |         X         |            X             |
| plotParameterSamples()    |          X          |           X         |          X       |         X         |            X             |
| plotParameterUncertainty()|          X          |           X         |          X       |         X         |            X             |
| plotPropertyMultiStarts() |          X          |           X         |          X       |         X         |            X             |
| plotPropertyProfiles()    |          X          |           X         |          X       |         X         |            X             |
| plotPropertySamples()     |          X          |           X         |          X       |         X         |            X             |
| plotPropertyUncertainty() |          X          |           X         |          X       |         X         |            X             |
| testGradient()            |          X          |           X         |          X       |         X         |            X             |
| collectResults()          |          X          |           X         |          X       |         X         |            X             |
