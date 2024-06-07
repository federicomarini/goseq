# goseq

The `goseq` package provides a method to detect Gene Ontology and/or other user 
defined categories which are over/under represented in RNA-seq data.

`goseq` can be found on Bioconductor
(<https://www.bioconductor.org/packages/goseq>).

## Installation

You can install the version of `goseq` which is on Bioconductor with these commands:

``` r
if (!require("BiocManager")) {
  install.packages("BiocManager")
}
BiocManager::install("goseq")
```

Alternatively, you can install the development version of `goseq` from GitHub with:

``` r
library("remotes")
remotes::install_github("federicomarini/goseq")
```

