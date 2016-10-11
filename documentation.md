PESTO Documentation {#doc}
===================

# Introduction         {#intro}

Computational models are commonly used in diverse disciplines such as computational biology, engineering, or meteorology. The parameterization of the these models is usually based on measurements or observations. The process of inferring model parameters from such data is called model calibration or parameter estimation. This parameter estimation is often not straightforward due to non-linearities in the model equations or due to the mere size and the resulting computational challenges. Therefore, efficient algorithms are required to provide robust results within acceptable time.

PESTO is a freely available Parameter EStimation TOolbox for MATLAB (MathWorks) implementing a number of state-of-the-art algorithms for parameter estimation. It provides the following features, which are explained in more detail [below](@ref features):

* Global optimization based on multi-start local optimization
* Parameter sampling
* Profile-likelihood analysis ...
* Visualization routines for all above analysis
* Parallel processing (requires [MATLAB Parallel Computing Toolbox](https://mathworks.com/products/parallel-computing/)) 
* ...

PESTO functions can be applied to any user-provided formulation of an optimization problem with an objective function that can be evaluated in MATLAB. Besides the objective function, upper and lower bounds for the function parameters need to be specified.

# Availability         {#availability}

PESTO can be freely obtained from https://github.com/ICB-DCM/PESTO/ by downloading the zip archive at https://github.com/ICB-DCM/PESTO/archive/master.zip or cloning the `git` repository via `
```
git clone git@github.com:ICB-DCM/PESTO.git
```

# Installation         {#installation}

If the zip archive was downloaded, it needs to be unzipped and the main folder has to added to the MATLAB search path. 

If the repository was cloned, the main folder needs to be added to the MATLAB search path.

Detailed instructions on how to modify your MATLAB search path are provided here: https://de.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html

# Licensing         {#licensing}

**TODO**

# How to cite         {#citation}

This section will be updated upon publication of PESTO.

# Code organisation         {#org}

The end-user interface is provided by the MATLAB functions and classes in the top-level directory. PESTO examples application are provided in `/examples/`. All other folders only contain files used internally in PESTO.

# Features ## {#features}

PESTO implements a number of state-of-the-art algorithms related to parameter estimations. The main features are described below. Various [Examples](@ref examples) demonstrate their application.

## Global optimization ### {#global-optimization}

Non-linear optimization problems can have multiple optima. Their number, location and values are usually not known in advance. Usually the user is interested in obtaining the global optimum. However, this would mean, that the whole parameter space has to be searched, which is computationally prohibitive.

One strategy for obtaining the global optimum is multi-start local optimization: Here, random points from across the parameter space are chosen as starting points for local optimization. If an adequate number of starting points spanning the whole parameter space has been chosen, the lowest/highest minimum/maximum is the global minimum/maximum.

This functionality is provided in getMultiStarts.m, getPropertyMultiStarts.m and the respective plotting routines plotMultiStarts.m and plotPropertyMultiStarts.m.

## Sensitivity analysis / Confidence Intervals ### {#confidence-intervals}

TODO: Explain Sensitivity analysis 

getParameterConfidenceIntervals()
* Threshold based confidence intervals using profile likelihood or local approximation and the probability mass
* Mass based confidence intervals from local approximation.
* Bayesian confidence interval Confidence intervals computed using [profile likelihood](@ref profile-likelihood)


## Profile likelihood ### {#profile-likelihood}

TODO

## Parameter sampling ### {#parameter-sampling}

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