PESTO 1.0 Documentation {#doc}
=======================

# Introduction         {#intro}

Computational models are commonly used in diverse disciplines such as computational biology, engineering, or meteorology. The parameterization of the these models is usually based on measurements or observations. The process of inferring model parameters from such data is called model calibration or parameter estimation. This parameter estimation is often not straightforward due to non-linearities in the model equations or due to the mere size and the resulting computational challenges. Therefore, efficient algorithms are required to provide robust results within acceptable time.

PESTO is a freely available Parameter EStimation TOolbox for MATLAB (MathWorks) implementing a number of state-of-the-art algorithms for parameter estimation. It provides the following features, which are explained in more detail [below](@ref features):

* Global optimization based on multi-start optimization
* Visualization routines for ... 
* Sampling routines for ....
* Profile-likelihood analysis ...
* Parallel processing (requires MATLAB Parallel Toolbox) 
* ...

PESTO functions can be applied to any user-provided model with an objective function that can be evaluated in MATLAB. The user needs to provide such an objective function as well as bounds for the function parameters and may specify a number of additional options.


# Availability         {#availability}

PESTO can be freely obtained from https://github.com/ICB-DCM/PESTO/ by downloading the an .... archive or cloning the `git` repository.

# Installation         {#installation}

If the packed archive was downloaded, it needs to be unzipped and the main folder has to added to the matlab path. 

If the repository was cloned, only the main folder needs to be added to the MATLAB path

# Licensing         {#licensing}

# How to cite         {#citation}

# Code organisation         {#org}

The end-user interface is provided by the MATLAB functions and classes in the top-level directory. PESTO examples application are provided in `/examples/`. All other folders contain files used internally in PESTO.

# Features ## {#features}

## Global optimization ### 

Non-linear optimization problems can have multiple optima. Their number, location and values are usually not known in advance. Usually the user is interested in obtaining the global optimum. However, this would mean, that the whole parameter space has to be searched, which is computationally prohibitive.

One strategy for obtaining global optima is the multi-start local optimization: Here random points from across the parameter space are chosen as starting points for local optimization. If an adequate number of starting points spanning the whole parameter space has been chosen, the lowest/highest minimum/maximum is the global minimum/maximum.

This functionality is provided in getMultiStarts.m, getPropertyMultiStarts.m and the respective plotting routines plotMultiStarts.m and plotPropertyMultiStarts.m.

## Sensitivity analysis / Confidence Intervals ### 

TODO: Explain Sensitivity analysis 

## Profile likelihood ### 

TODO

## Plotting ### 

TODO

## Properties ### 

Explain concept of "properties"
A property is any function g(theta) ... could also be expressed as state variable, but... 

The following functions are available to analyze properties:
* getPropertyConfidenceIntervals.m
* getPropertyMultiStarts.m
* getPropertyProfiles.m
* getPropertySamples.m
* plotPropertyMultiStarts.m
* plotPropertyProfiles.m
* plotPropertySamples.m
* plotPropertyUncertainty.m


# TODO ##

As opposed to other toolboxes ... 


The following table show functionality of PESTO in comparison to other available toolboxes:

| Feature \ Toolbox  | *PESTO*    |  A           | B         | 
|:-------------------|:----------:|:------------:|:---------:|
| Global optimization|     X      |              |           |
| Profile likelihood |     X      |              |           |
| Feature 1          |            |              |           |
| Feature 2          |            |              |           |

# References 

... 