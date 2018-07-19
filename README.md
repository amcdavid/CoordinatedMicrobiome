# CoordinatedMicrobiome
Analysis for "Neonatal gut and respiratory microbiota: coordinated development through time and space"

# Todo

  - ensure no absolute paths
  - CCA.Rmd should use data in data/ rather than it's own .biom files

# Contents

  - CCA.Rmd : Canonical correlation analysis
  - metacommunity_post_analysis.Rmd : Further posthoc analysis of CSTs including Generalized additive models
  - data/ : OTU tables and metadata (need to make sure scripts are self-contained)
  - Composition modeling/ : regression models for each OTU and CST-OTU prediction models
  - DMN Cluster Selection/ : Dirichlet Multionomial model search
  - SingleIndex/ : R package used to fit single index models

# Required libraries/software

R 3.4.2 and the following:

  library(ComplexHeatmap)
  library(DT)
  library(DirichletMultinomial)
  library(GGally)
  library(SummarizedExperiment)
  library(biomformat)
  library(broom)
  library(circlize)
  library(data.table)
  library(dplyr)
  library(dtplyr)
  library(forcats)
  library(ggplot2)
  library(ggrepel)
  library(ggthemes)
  library(mgcv)
  library(phyloseq)
  library(readr)
  library(readxl)
  library(reshape2)
  library(stringr)
  library(tidyr)
  library(tidyverse)
  library(whoami)
  library(knitr)
  library(plyr)
  library(purrr)
  library(splines)
  library(cowplot)
  library(tibble)
  library(viridis)
  library(circlize)
  library(devtools)