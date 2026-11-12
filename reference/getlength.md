# Retrieves Gene length data

Gets the length of each gene in a vector.

## Usage

``` r
getlength(genes, genome, id)
```

## Arguments

- genes:

  A vector or list of the genes for which length information is
  required.

- genome:

  A string identifying the genome that `genes` refer to. For a list of
  supported organisms run
  [`supportedGenomes`](https://rdrr.io/pkg/geneLenDataBase/man/supportedGenomes.html).

- id:

  A string identifying the gene identifier used by `genes`. For a list
  of supported gene IDs run
  [`supportedGeneIDs`](https://rdrr.io/pkg/geneLenDataBase/man/supportedGeneIDs.html).

## Value

Returns a vector of the gene lengths, in the same order as `genes`. If
length data is unavailable for a particular gene NA is returned in that
position. The returned vector is intended for use with the `bias.data`
option of the [`nullp`](nullp.md) function.

## Details

Length data is obtained from data obtained from the UCSC genome browser
for each combination of `genome` and `id`. As fetching this data at
runtime is time consuming, a local copy of the length information for
common genomes and gene ID are included in the geneLenDataBase package.
This function uses this package to fetch the required data.

The length of a gene is taken to be the median length of all its mature,
mRNA, transcripts. It is always preferable to obtain length information
directly for the gene ID used to summarize your count data, rather than
converting IDs and then using the supplied databases. Even when two
genes have a one-to-one mapping between different identifier conventions
(which is often not the case), they frequently refer to slightly
different regions of the genome with different lengths. It is therefore
recommended that the user perform the full analysis in terms of only one
gene ID, or manually obtain their own length data for the identifier
used to bin reads by gene.

## See also

[`supportedGenomes`](https://rdrr.io/pkg/geneLenDataBase/man/supportedGenomes.html),
[`supportedGeneIDs`](https://rdrr.io/pkg/geneLenDataBase/man/supportedGeneIDs.html),
[`nullp`](nullp.md), geneLenDataBase

## Author

Matthew D. Young <myoung@wehi.edu.au>

## Examples

``` r

genes <- c("ENSG00000124208", 
           "ENSG00000182463", 
           "ENSG00000124201", 
           "ENSG00000124205", 
           "ENSG00000124207")
getlength(genes,'hg19','ensGene')
#> Loading hg19 length data...
#> [1] 1978.0 3133.0 2973.0 2593.0 3036.5
```
