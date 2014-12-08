Description
=============

SpatialCondition

* Is an R package for estimating spatial and spatiotemporal variation in individual condition using Gaussian random fields and potentially including covariates
* Is intended to allow rapid comparison of condition factor among multiple species.

Instructions
=============
First, please install JAGS (Just Another Gibbs Sampler) here: http://mcmc-jags.sourceforge.net/

Next, please use R version >=3.1.1 and install the package:


    # Install package
    install.packages("devtools")
    library("devtools")
    install_github("James-Thorson/spatial_condition_factor")
    # Load package
    library(SpatialDeltaGLMM)

Please see examples folder for an example of how to run the model:
https://github.com/James-Thorson/spatial_condition_factor/blob/master/examples/example_for_NWFSC_shelf-slope_data.R

Known installation/usage issues
=============
none

Further reading
=============

For more details regarding development and testing of this delta-GLMM software please see:
* Thorson, J. T., A. O. Shelton, E. J. Ward, and H. Skaug. In press. Geostatistical delta-generalized linear mixed models improve precision for estimated abundance indices for West Coast groundfishes. ICES Journal of Marine Science.
* Thorson, J. T., H. Skaug, K. Kristensen, A. O. Shelton, E. J. Ward, J. Harms, and J. Benante. In press. The importance of spatial models for estimating the strength of density dependence. Ecology.

