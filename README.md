## PhenoPath

PhenoPath learns genomic trajectories (pseudotimes) in the presence of heterogenous environmental and genetic backgrounds encoded as additional covariates and identifies interactions between the trajectories and covariates. Scalable variational Bayesian inference allows the trajectory and interactions to be inferred for thousands of samples and genes quickly.

[![Build Status](https://travis-ci.org/kieranrcampbell/phenopath.svg?branch=master)](https://travis-ci.org/kieranrcampbell/phenopath)

To quickly get started, see either of the vignettes:

* [Introduction to PhenoPath](https://kieranrcampbell.github.io/phenopath/introduction_to_phenopath.html)
* [End-to-end example on Shalek et al. (2014) from raw data to stimulant-pseudotime interactions](https://kieranrcampbell.github.io/phenopath/phenopath_shalek_vignette.html)

### Installation

```r
install.packages("devtools") # If not already installed
devtools::install_github("kieranrcampbell/phenopath", build_vignettes = TRUE)
```

### Overview

PhenoPath models the observed expression *y* in terms of a latent pathway score (pseudotime) *z*. Uniquely, the evolution of genes along the trajectory isn't common to each gene but can be perturbed by an additional sample-specific covariate. For example, this could be the mutational status of each sample or a drug that each sample was exposed to.

<img src="https://user-images.githubusercontent.com/2039489/26843141-1deacabc-4ae7-11e7-8630-a6a6e425ca0b.png" width="500"/>

#### Inference

Inference is performed using co-ordinate ascent variational inference (CAVI) implemented using the Rcpp library for increased performance. This minimises the KL-divergence between a set of approximating distributions and the true posterior by making a full-factorized mean-field approximation.


### Getting started

See the vignette for details.



### Authors

Kieran R Campbell & Christopher Yau