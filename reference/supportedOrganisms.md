# Supported Organisms

Lists which genomes and gene ids are supported for gene lengths and for
gene ontology

## Usage

``` r
supportedOrganisms()
```

## Value

A data.frame containing supported genomes and gene ids

## Details

Goseq allows a user to provide their own bias data (usually gene
lengths) and/or gene categories (usually gene ontologies), but goseq
also provides this data automatically for many commonly used species.
This function lists which genome and gene ids are automatically
supported by goseq. The first to third columns list the genomes, gene
ids, and gene id descriptions respectively. The fourth column indicates
whether this combination of genome and id are available in the
geneLengthDataBase. If a particular combination is absent, goseq may
still automatically fetch the gene lengths from either a TxDB annotation
package (must be installed) or download the data from UCSC. For example
gene lengths for `hg38` are not supported in `geneLengthDataBase` but
may be fetched by these other means. However, this is not always the
case. The final column indicates if GO terms will be automatically
fetched for the genome and id combination. Goseq relies on an `org`
annotation package (e.g. `org.Hs.eg.db`) existing for the organism. In
general, if either the lengths or GO terms are not supported, the user
must enter this information manually.

## Author

Nadia Davidson <nadia.davidson@mcri.edu.au>

## Examples

``` r
supportedOrganisms()
#> Warning: 'getMethods' is deprecated.
#> Use 'getMethodsForDispatch(f, TRUE)' instead.
#> See help("Deprecated")
#>            Genome         Id       Id Description Lengths in geneLeneDataBase
#> 139   Arabidopsis                                                       FALSE
#> 140   E. coli K12                                                       FALSE
#> 141 E. coli Sakai                                                       FALSE
#> 142       Malaria                                                       FALSE
#> 1         anoCar1    ensGene      Ensembl gene ID                        TRUE
#> 2         anoGam1    ensGene      Ensembl gene ID                        TRUE
#> 98        anoGam3                                                       FALSE
#> 3         apiMel2    ensGene      Ensembl gene ID                        TRUE
#> 56        bosTau2 geneSymbol          Gene Symbol                        TRUE
#> 4         bosTau3    ensGene      Ensembl gene ID                        TRUE
#> 57        bosTau3 geneSymbol          Gene Symbol                        TRUE
#> 5         bosTau4    ensGene      Ensembl gene ID                        TRUE
#> 58        bosTau4 geneSymbol          Gene Symbol                        TRUE
#> 99        bosTau6                                                       FALSE
#> 100       bosTau7                                                       FALSE
#> 101       bosTau8                                                       FALSE
#> 102       bosTau9                                                       FALSE
#> 6         canFam1    ensGene      Ensembl gene ID                        TRUE
#> 59        canFam1 geneSymbol          Gene Symbol                        TRUE
#> 7         canFam2    ensGene      Ensembl gene ID                        TRUE
#> 60        canFam2 geneSymbol          Gene Symbol                        TRUE
#> 103       canFam3                                                       FALSE
#> 104       canFam4                                                       FALSE
#> 105       canFam5                                                       FALSE
#> 106       canFam6                                                       FALSE
#> 8         cavPor3    ensGene      Ensembl gene ID                        TRUE
#> 107          ce10                                                       FALSE
#> 108          ce11                                                       FALSE
#> 61            ce2 geneSymbol          Gene Symbol                        TRUE
#> 62            ce4 geneSymbol          Gene Symbol                        TRUE
#> 9             ce6    ensGene      Ensembl gene ID                        TRUE
#> 63            ce6 geneSymbol          Gene Symbol                        TRUE
#> 64            ci1 geneSymbol          Gene Symbol                        TRUE
#> 10            ci2    ensGene      Ensembl gene ID                        TRUE
#> 65            ci2 geneSymbol          Gene Symbol                        TRUE
#> 109      danRer10                                                       FALSE
#> 110      danRer11                                                       FALSE
#> 11        danRer3    ensGene      Ensembl gene ID                        TRUE
#> 66        danRer3 geneSymbol          Gene Symbol                        TRUE
#> 12        danRer4    ensGene      Ensembl gene ID                        TRUE
#> 67        danRer4 geneSymbol          Gene Symbol                        TRUE
#> 13        danRer5    ensGene      Ensembl gene ID                        TRUE
#> 54        danRer5   vegaGene HAVANA Pseudogene ID                        TRUE
#> 68        danRer5 geneSymbol          Gene Symbol                        TRUE
#> 14        danRer6    ensGene      Ensembl gene ID                        TRUE
#> 69        danRer6 geneSymbol          Gene Symbol                        TRUE
#> 111       danRer7                                                       FALSE
#> 70            dm1 geneSymbol          Gene Symbol                        TRUE
#> 71            dm2 geneSymbol          Gene Symbol                        TRUE
#> 72            dm3 geneSymbol          Gene Symbol                        TRUE
#> 112           dm6                                                       FALSE
#> 73        equCab1 geneSymbol          Gene Symbol                        TRUE
#> 15        equCab2    ensGene      Ensembl gene ID                        TRUE
#> 74        equCab2 geneSymbol          Gene Symbol                        TRUE
#> 16        felCat3    ensGene      Ensembl gene ID                        TRUE
#> 75        felCat3 geneSymbol          Gene Symbol                        TRUE
#> 17            fr1    ensGene      Ensembl gene ID                        TRUE
#> 18            fr2    ensGene      Ensembl gene ID                        TRUE
#> 19        galGal2    ensGene      Ensembl gene ID                        TRUE
#> 76        galGal2 geneSymbol          Gene Symbol                        TRUE
#> 20        galGal3    ensGene      Ensembl gene ID                        TRUE
#> 77        galGal3 geneSymbol          Gene Symbol                        TRUE
#> 113       galGal4                                                       FALSE
#> 114       galGal5                                                       FALSE
#> 115       galGal6                                                       FALSE
#> 21        gasAcu1    ensGene      Ensembl gene ID                        TRUE
#> 22           hg16    ensGene      Ensembl gene ID                        TRUE
#> 45           hg16  knownGene       Entrez Gene ID                        TRUE
#> 78           hg16 geneSymbol          Gene Symbol                        TRUE
#> 23           hg17    ensGene      Ensembl gene ID                        TRUE
#> 46           hg17  knownGene       Entrez Gene ID                        TRUE
#> 55           hg17   vegaGene HAVANA Pseudogene ID                        TRUE
#> 79           hg17 geneSymbol          Gene Symbol                        TRUE
#> 24           hg18    ensGene      Ensembl gene ID                        TRUE
#> 47           hg18  knownGene       Entrez Gene ID                        TRUE
#> 80           hg18 geneSymbol          Gene Symbol                        TRUE
#> 25           hg19    ensGene      Ensembl gene ID                        TRUE
#> 48           hg19  knownGene       Entrez Gene ID                        TRUE
#> 81           hg19 geneSymbol          Gene Symbol                        TRUE
#> 116          hg38                                                       FALSE
#> 117          mm10                                                       FALSE
#> 118          mm39                                                       FALSE
#> 26            mm7    ensGene      Ensembl gene ID                        TRUE
#> 49            mm7  knownGene       Entrez Gene ID                        TRUE
#> 82            mm7 geneSymbol          Gene Symbol                        TRUE
#> 27            mm8    ensGene      Ensembl gene ID                        TRUE
#> 50            mm8  knownGene       Entrez Gene ID                        TRUE
#> 83            mm8 geneSymbol          Gene Symbol                        TRUE
#> 28            mm9    ensGene      Ensembl gene ID                        TRUE
#> 51            mm9  knownGene       Entrez Gene ID                        TRUE
#> 84            mm9 geneSymbol          Gene Symbol                        TRUE
#> 29        monDom4    ensGene      Ensembl gene ID                        TRUE
#> 85        monDom4 geneSymbol          Gene Symbol                        TRUE
#> 30        monDom5    ensGene      Ensembl gene ID                        TRUE
#> 86        monDom5 geneSymbol          Gene Symbol                        TRUE
#> 31        ornAna1    ensGene      Ensembl gene ID                        TRUE
#> 87        ornAna1 geneSymbol          Gene Symbol                        TRUE
#> 32        oryLat2    ensGene      Ensembl gene ID                        TRUE
#> 88        oryLat2 geneSymbol          Gene Symbol                        TRUE
#> 33        panTro1    ensGene      Ensembl gene ID                        TRUE
#> 34        panTro2    ensGene      Ensembl gene ID                        TRUE
#> 89        panTro2 geneSymbol          Gene Symbol                        TRUE
#> 119       panTro3                                                       FALSE
#> 120       panTro4                                                       FALSE
#> 121       panTro5                                                       FALSE
#> 122       panTro6                                                       FALSE
#> 35        ponAbe2    ensGene      Ensembl gene ID                        TRUE
#> 90        ponAbe2 geneSymbol          Gene Symbol                        TRUE
#> 123      rheMac10                                                       FALSE
#> 36        rheMac2    ensGene      Ensembl gene ID                        TRUE
#> 91        rheMac2 geneSymbol          Gene Symbol                        TRUE
#> 124       rheMac3                                                       FALSE
#> 125       rheMac8                                                       FALSE
#> 37            rn3    ensGene      Ensembl gene ID                        TRUE
#> 52            rn3  knownGene       Entrez Gene ID                        TRUE
#> 92            rn3 geneSymbol          Gene Symbol                        TRUE
#> 38            rn4    ensGene      Ensembl gene ID                        TRUE
#> 53            rn4  knownGene       Entrez Gene ID                        TRUE
#> 93            rn4 geneSymbol          Gene Symbol                        TRUE
#> 126           rn5                                                       FALSE
#> 127           rn6                                                       FALSE
#> 128           rn7                                                       FALSE
#> 129           rn8                                                       FALSE
#> 39        sacCer1    ensGene      Ensembl gene ID                        TRUE
#> 40        sacCer2    ensGene      Ensembl gene ID                        TRUE
#> 130       sacCer3                                                       FALSE
#> 94        strPur1 geneSymbol          Gene Symbol                        TRUE
#> 95        strPur2 geneSymbol          Gene Symbol                        TRUE
#> 131      susScr11                                                       FALSE
#> 132       susScr2                                                       FALSE
#> 133       susScr3                                                       FALSE
#> 41        taeGut1    ensGene      Ensembl gene ID                        TRUE
#> 96        taeGut1 geneSymbol          Gene Symbol                        TRUE
#> 42        tetNig1    ensGene      Ensembl gene ID                        TRUE
#> 43        tetNig2    ensGene      Ensembl gene ID                        TRUE
#> 134       xenTro1                                                       FALSE
#> 135      xenTro10                                                       FALSE
#> 44        xenTro2    ensGene      Ensembl gene ID                        TRUE
#> 97        xenTro2 geneSymbol          Gene Symbol                        TRUE
#> 136       xenTro3                                                       FALSE
#> 137       xenTro7                                                       FALSE
#> 138       xenTro9                                                       FALSE
#>     GO Annotation Available
#> 139                    TRUE
#> 140                    TRUE
#> 141                    TRUE
#> 142                    TRUE
#> 1                     FALSE
#> 2                      TRUE
#> 98                     TRUE
#> 3                     FALSE
#> 56                     TRUE
#> 4                      TRUE
#> 57                     TRUE
#> 5                      TRUE
#> 58                     TRUE
#> 99                     TRUE
#> 100                    TRUE
#> 101                    TRUE
#> 102                    TRUE
#> 6                      TRUE
#> 59                     TRUE
#> 7                      TRUE
#> 60                     TRUE
#> 103                    TRUE
#> 104                    TRUE
#> 105                    TRUE
#> 106                    TRUE
#> 8                     FALSE
#> 107                    TRUE
#> 108                    TRUE
#> 61                     TRUE
#> 62                     TRUE
#> 9                      TRUE
#> 63                     TRUE
#> 64                    FALSE
#> 10                    FALSE
#> 65                    FALSE
#> 109                    TRUE
#> 110                    TRUE
#> 11                     TRUE
#> 66                     TRUE
#> 12                     TRUE
#> 67                     TRUE
#> 13                     TRUE
#> 54                    FALSE
#> 68                     TRUE
#> 14                     TRUE
#> 69                     TRUE
#> 111                    TRUE
#> 70                     TRUE
#> 71                     TRUE
#> 72                     TRUE
#> 112                    TRUE
#> 73                    FALSE
#> 15                    FALSE
#> 74                    FALSE
#> 16                    FALSE
#> 75                    FALSE
#> 17                    FALSE
#> 18                    FALSE
#> 19                     TRUE
#> 76                     TRUE
#> 20                     TRUE
#> 77                     TRUE
#> 113                    TRUE
#> 114                    TRUE
#> 115                    TRUE
#> 21                    FALSE
#> 22                     TRUE
#> 45                     TRUE
#> 78                     TRUE
#> 23                     TRUE
#> 46                     TRUE
#> 55                    FALSE
#> 79                     TRUE
#> 24                     TRUE
#> 47                     TRUE
#> 80                     TRUE
#> 25                     TRUE
#> 48                     TRUE
#> 81                     TRUE
#> 116                    TRUE
#> 117                    TRUE
#> 118                    TRUE
#> 26                     TRUE
#> 49                     TRUE
#> 82                     TRUE
#> 27                     TRUE
#> 50                     TRUE
#> 83                     TRUE
#> 28                     TRUE
#> 51                     TRUE
#> 84                     TRUE
#> 29                    FALSE
#> 85                    FALSE
#> 30                    FALSE
#> 86                    FALSE
#> 31                    FALSE
#> 87                    FALSE
#> 32                    FALSE
#> 88                    FALSE
#> 33                     TRUE
#> 34                     TRUE
#> 89                     TRUE
#> 119                    TRUE
#> 120                    TRUE
#> 121                    TRUE
#> 122                    TRUE
#> 35                    FALSE
#> 90                    FALSE
#> 123                    TRUE
#> 36                     TRUE
#> 91                     TRUE
#> 124                    TRUE
#> 125                    TRUE
#> 37                     TRUE
#> 52                     TRUE
#> 92                     TRUE
#> 38                     TRUE
#> 53                     TRUE
#> 93                     TRUE
#> 126                    TRUE
#> 127                    TRUE
#> 128                    TRUE
#> 129                    TRUE
#> 39                     TRUE
#> 40                     TRUE
#> 130                    TRUE
#> 94                    FALSE
#> 95                    FALSE
#> 131                    TRUE
#> 132                    TRUE
#> 133                    TRUE
#> 41                    FALSE
#> 96                    FALSE
#> 42                    FALSE
#> 43                    FALSE
#> 134                    TRUE
#> 135                    TRUE
#> 44                     TRUE
#> 97                     TRUE
#> 136                    TRUE
#> 137                    TRUE
#> 138                    TRUE
```
