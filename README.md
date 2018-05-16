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
