# Androgen stimulation of prostate cancer Cell lines.

This data set gives the RNA-seq data from an experiment measuring the
effects of androgen stimulation on prostate cancer. Information is given
about all (ENSEMBL) genes for which there was at least one mapping read
in either the treated or untreated RNA-seq experiment. The `edgeR`
package was used to determine which genes were differentially expressed.
The details of the analysis can be found in the `goseq` vignette.

## Format

A named vector of ENSEMBL genes, with 1 representing differential
expression.

## Source

Determination of tag density required for digital transcriptome
analysis: application to an androgen-sensitive prostate cancer model,
2008, Li et. al.

## References

Li, H., Lovci, M. T., Kwon, Y. S., Rosenfeld, M. G., Fu, X. D., Yeo, G.
W. (2008) *Determination of tag density required for digital
transcriptome analysis: application to an androgen-sensitive prostate
cancer model* Proceedings of the National Academy of Sciences of the
United States of America Date: Dec 23 Vol: 105 Issue: 51 Pages: 20179-84

## Examples

``` r

data(genes)
head(genes)
#> ENSG00000230758 ENSG00000182463 ENSG00000124208 ENSG00000230753 ENSG00000224628 
#>               0               0               0               0               0 
#> ENSG00000125835 
#>               0 
```
