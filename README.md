# CoordinatedMicrobiome
Analysis for [Neonatal gut and respiratory microbiota: coordinated development through time and space](https://www.biorxiv.org/content/early/2018/01/19/247122).

Requires R 3.4.3 and Bioconductor 3.5.

# Contents

  - CCA.Rmd : Canonical correlation analysis
  - metacommunity_post_analysis.Rmd : Further posthoc analysis of CSTs including Generalized additive models
  - data/ : OTU tables and metadata
  - Composition modeling/ : regression models for each OTU and CST-OTU prediction models
  - DMN Cluster Selection/ : Dirichlet Multionomial model search
  - SingleIndex/ : R package used to fit single index models

# To run
1.  Clone this repo
2.  From the repo root, run
```r
install("CoordinatedMBDeps", dependencies = TRUE)
```
This will install all needed dependencies. **Caution: this might downgrade your bioconductor installation to version 3.5.**
Before starting R you may wish to set the environment variable `R_LIBS_USER` -- or within R `.libPaths()` -- to an unused directory before proceeding.  Or use the bioconductor 3.5 AMI or docker images...

3.  `rmarkdown::render('Composition Modeling/pairwise_composition_modeling.Rmd')`
4.  `rmarkdown::render('Composition Modeling/get_model_pvals_betas.Rmd')` to generate univariate associations for each OTU and CSTs
5.  `rmarkdown::render('DMN Cluster Selection/normalized_cluster_selection.Rmd')` to run DMN resampling procedure
5.  `rmarkdown::render('CCA.Rmd')` for canonical correlation analysis
6.  `rmarkdown::render('metacommunity_post_analysis.Rmd')` for remaining plots

# Updates
 - Added pseudopackage to collect dependencies
 - Removed some obsolete code
 - Added some missing data files
