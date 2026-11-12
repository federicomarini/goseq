# Fetch GO categories

Obtains all gene ontology (GO) categories associated with a set of genes
using the relevant organism package.

## Usage

``` r
getgo(genes, genome, id, fetch.cats = c("GO:CC", "GO:BP", "GO:MF"))
```

## Arguments

- genes:

  A vector or list of genes to get the associated GO categories.

- genome:

  A string identifying the genome that `genes` refer to. For a list of
  supported organisms run
  [`supportedGenomes`](https://rdrr.io/pkg/geneLenDataBase/man/supportedGenomes.html).

- id:

  A string identifying the gene identifier used by `genes`. For a list
  of supported gene IDs run
  [`supportedGeneIDs`](https://rdrr.io/pkg/geneLenDataBase/man/supportedGeneIDs.html).

- fetch.cats:

  A vector specifying which categories to fetch the mapping between
  category names and genes for. See details for valid options.

## Value

A list where each entry is named by a gene and contains a vector of all
the associated GO categories. This can be used directly with the
`gene2cat` option in [`goseq`](goseq.md).

## Details

This function attempts to make use of the organism packages
(org.\<Genome\>.\<GeneID\>.db) to obtain the mapping between gene ID and
GO categories. As with [`getlength`](getlength.md) it is preferable that
the same gene identifier system is used for both summarization and
retrieving GO categories.

Valid options for the `fetch.cats` argument are any combination of
"GO:CC", "GO:BP", "GO:MF" & "KEGG". The three GO terms refer to the
Cellular Component, Biological Process and Molecular Function
respectively. "KEGG" refers to KEGG pathways.

Note that `getgo` is a convenience function, designed to make extracting
mappings between GO categories and Gene ID easy. For less common
organisms and/or gene ID `getgo` may fail to return a mapping even when
a legitimate mapping exists in the relevant organism package. If `getgo`
fails, you should always try to build the mapping yourself from the
organism package (if one exists) before deciding that the information is
unavailable. Further information and examples of this can be found in
the package Vignette.

## See also

[`supportedGenomes`](https://rdrr.io/pkg/geneLenDataBase/man/supportedGenomes.html),
[`supportedGeneIDs`](https://rdrr.io/pkg/geneLenDataBase/man/supportedGeneIDs.html),
[`goseq`](goseq.md)

## Author

Matthew D. Young <myoung@wehi.edu.au>

## Examples

``` r

genes <- c("ENSG00000124208", 
           "ENSG00000182463", 
           "ENSG00000124201", 
           "ENSG00000124205", 
           "ENSG00000124207")
getgo(genes,'hg19','ensGene')
#> Loading required package: AnnotationDbi
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following object is masked from ‘package:utils’:
#> 
#>     data
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, scale, sequence, table,
#>     tapply, transform, unique, unsplit, which.max, which.min
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> Loading required package: IRanges
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:geneLenDataBase’:
#> 
#>     unfactor
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> 
#> $<NA>
#> NULL
#> 
#> $ENSG00000182463
#>  [1] "GO:0006139" "GO:0006351" "GO:0006355" "GO:0006357" "GO:0006366"
#>  [6] "GO:0008150" "GO:0008152" "GO:0009058" "GO:0009059" "GO:0009889"
#> [11] "GO:0009987" "GO:0010467" "GO:0010468" "GO:0010556" "GO:0016070"
#> [16] "GO:0019219" "GO:0019222" "GO:0032774" "GO:0034654" "GO:0043170"
#> [21] "GO:0044238" "GO:0050789" "GO:0050794" "GO:0051252" "GO:0060255"
#> [26] "GO:0065007" "GO:0080090" "GO:0090304" "GO:0141187" "GO:2001141"
#> [31] "GO:0000785" "GO:0005575" "GO:0005622" "GO:0005634" "GO:0005694"
#> [36] "GO:0043226" "GO:0043227" "GO:0043228" "GO:0043229" "GO:0043231"
#> [41] "GO:0043232" "GO:0110165" "GO:0000981" "GO:0003674" "GO:0003676"
#> [46] "GO:0003677" "GO:0003700" "GO:0005488" "GO:0005515" "GO:0140110"
#> 
#> $ENSG00000124201
#>   [1] "GO:0002218" "GO:0002221" "GO:0002253" "GO:0002376" "GO:0002682"
#>   [6] "GO:0002683" "GO:0002684" "GO:0002753" "GO:0002757" "GO:0002758"
#>  [11] "GO:0002764" "GO:0002831" "GO:0002833" "GO:0006325" "GO:0006338"
#>  [16] "GO:0006950" "GO:0006952" "GO:0006954" "GO:0006955" "GO:0007154"
#>  [21] "GO:0007165" "GO:0008150" "GO:0008152" "GO:0009058" "GO:0009059"
#>  [26] "GO:0009605" "GO:0009607" "GO:0009615" "GO:0009617" "GO:0009889"
#>  [31] "GO:0009890" "GO:0009892" "GO:0009966" "GO:0009968" "GO:0009987"
#>  [36] "GO:0010467" "GO:0010468" "GO:0010556" "GO:0010558" "GO:0010605"
#>  [41] "GO:0010629" "GO:0010646" "GO:0010648" "GO:0016032" "GO:0016043"
#>  [46] "GO:0019058" "GO:0019079" "GO:0019222" "GO:0022607" "GO:0023051"
#>  [51] "GO:0023052" "GO:0023057" "GO:0030522" "GO:0031047" "GO:0031048"
#>  [56] "GO:0031333" "GO:0031347" "GO:0031349" "GO:0031507" "GO:0032101"
#>  [61] "GO:0032103" "GO:0035556" "GO:0039531" "GO:0039532" "GO:0040029"
#>  [66] "GO:0042742" "GO:0043170" "GO:0043207" "GO:0043254" "GO:0043933"
#>  [71] "GO:0044085" "GO:0044087" "GO:0044419" "GO:0044546" "GO:0045069"
#>  [76] "GO:0045071" "GO:0045087" "GO:0045088" "GO:0045089" "GO:0045814"
#>  [81] "GO:0048518" "GO:0048519" "GO:0048523" "GO:0048525" "GO:0048583"
#>  [86] "GO:0048584" "GO:0048585" "GO:0050776" "GO:0050778" "GO:0050789"
#>  [91] "GO:0050792" "GO:0050794" "GO:0050896" "GO:0051128" "GO:0051129"
#>  [96] "GO:0051607" "GO:0051707" "GO:0051716" "GO:0060255" "GO:0062207"
#> [101] "GO:0065003" "GO:0065007" "GO:0070828" "GO:0071840" "GO:0080134"
#> [106] "GO:0098542" "GO:0140546" "GO:0140632" "GO:0141084" "GO:0141085"
#> [111] "GO:0141086" "GO:1900225" "GO:1900226" "GO:1902531" "GO:1902532"
#> [116] "GO:1903900" "GO:0005575" "GO:0005622" "GO:0005634" "GO:0005737"
#> [121] "GO:0005739" "GO:0005740" "GO:0005741" "GO:0005829" "GO:0010494"
#> [126] "GO:0016020" "GO:0019867" "GO:0030880" "GO:0031090" "GO:0031379"
#> [131] "GO:0031380" "GO:0031966" "GO:0031967" "GO:0031968" "GO:0032991"
#> [136] "GO:0035770" "GO:0036464" "GO:0043226" "GO:0043227" "GO:0043228"
#> [141] "GO:0043229" "GO:0043231" "GO:0043232" "GO:0061695" "GO:0098588"
#> [146] "GO:0099080" "GO:0110165" "GO:0140513" "GO:0140535" "GO:1902494"
#> [151] "GO:1990234" "GO:0003674" "GO:0003676" "GO:0003723" "GO:0003824"
#> [156] "GO:0004386" "GO:0004842" "GO:0005488" "GO:0008270" "GO:0016740"
#> [161] "GO:0016746" "GO:0016755" "GO:0016853" "GO:0019787" "GO:0036094"
#> [166] "GO:0043167" "GO:0043169" "GO:0046872" "GO:0046914" "GO:0061630"
#> [171] "GO:0061659" "GO:0120543" "GO:0120545" "GO:0140096" "GO:0140311"
#> [176] "GO:0140313" "GO:0140640" "GO:0140657"
#> 
#> $ENSG00000124205
#>   [1] "GO:0000165" "GO:0000278" "GO:0000280" "GO:0000902" "GO:0001990"
#>   [6] "GO:0002027" "GO:0002376" "GO:0002682" "GO:0002684" "GO:0002685"
#>  [11] "GO:0002687" "GO:0002688" "GO:0002690" "GO:0002790" "GO:0003008"
#>  [16] "GO:0003012" "GO:0003013" "GO:0003015" "GO:0003018" "GO:0003044"
#>  [21] "GO:0003073" "GO:0003100" "GO:0006810" "GO:0006873" "GO:0006874"
#>  [26] "GO:0006935" "GO:0006936" "GO:0006937" "GO:0006939" "GO:0006940"
#>  [31] "GO:0006996" "GO:0007049" "GO:0007088" "GO:0007154" "GO:0007165"
#>  [36] "GO:0007166" "GO:0007267" "GO:0007275" "GO:0007346" "GO:0007399"
#>  [41] "GO:0007409" "GO:0007411" "GO:0008015" "GO:0008016" "GO:0008150"
#>  [46] "GO:0008152" "GO:0008217" "GO:0008283" "GO:0008284" "GO:0009058"
#>  [51] "GO:0009059" "GO:0009605" "GO:0009653" "GO:0009889" "GO:0009914"
#>  [56] "GO:0009966" "GO:0009967" "GO:0009987" "GO:0010460" "GO:0010467"
#>  [61] "GO:0010468" "GO:0010556" "GO:0010564" "GO:0010638" "GO:0010646"
#>  [66] "GO:0010647" "GO:0010817" "GO:0014821" "GO:0014826" "GO:0014829"
#>  [71] "GO:0015833" "GO:0016043" "GO:0016049" "GO:0016477" "GO:0019222"
#>  [76] "GO:0019229" "GO:0019725" "GO:0022008" "GO:0022402" "GO:0023051"
#>  [81] "GO:0023052" "GO:0023056" "GO:0023061" "GO:0030003" "GO:0030030"
#>  [86] "GO:0030072" "GO:0030154" "GO:0030182" "GO:0030334" "GO:0030335"
#>  [91] "GO:0030593" "GO:0030595" "GO:0031175" "GO:0032101" "GO:0032103"
#>  [96] "GO:0032501" "GO:0032502" "GO:0032879" "GO:0032940" "GO:0033043"
#> [101] "GO:0035150" "GO:0035296" "GO:0035556" "GO:0040007" "GO:0040011"
#> [106] "GO:0040012" "GO:0040017" "GO:0042127" "GO:0042221" "GO:0042310"
#> [111] "GO:0042330" "GO:0042592" "GO:0043170" "GO:0043408" "GO:0043410"
#> [116] "GO:0044057" "GO:0045595" "GO:0045597" "GO:0045787" "GO:0045823"
#> [121] "GO:0045840" "GO:0045933" "GO:0045987" "GO:0046879" "GO:0046883"
#> [126] "GO:0046887" "GO:0046903" "GO:0048285" "GO:0048468" "GO:0048518"
#> [131] "GO:0048522" "GO:0048583" "GO:0048584" "GO:0048588" "GO:0048589"
#> [136] "GO:0048666" "GO:0048667" "GO:0048675" "GO:0048699" "GO:0048731"
#> [141] "GO:0048812" "GO:0048856" "GO:0048858" "GO:0048869" "GO:0048870"
#> [146] "GO:0048878" "GO:0050789" "GO:0050793" "GO:0050794" "GO:0050801"
#> [151] "GO:0050886" "GO:0050896" "GO:0050900" "GO:0050920" "GO:0050921"
#> [156] "GO:0051046" "GO:0051047" "GO:0051049" "GO:0051050" "GO:0051094"
#> [161] "GO:0051128" "GO:0051130" "GO:0051179" "GO:0051234" "GO:0051239"
#> [166] "GO:0051240" "GO:0051716" "GO:0051726" "GO:0051783" "GO:0051785"
#> [171] "GO:0055074" "GO:0055080" "GO:0055082" "GO:0060047" "GO:0060255"
#> [176] "GO:0060326" "GO:0060560" "GO:0061564" "GO:0065007" "GO:0065008"
#> [181] "GO:0070887" "GO:0071621" "GO:0071705" "GO:0071840" "GO:0090066"
#> [186] "GO:0090068" "GO:0090257" "GO:0097485" "GO:0097529" "GO:0097530"
#> [191] "GO:0097746" "GO:0098771" "GO:0120036" "GO:0120039" "GO:0140014"
#> [196] "GO:0140352" "GO:0141124" "GO:0160108" "GO:1902531" "GO:1902533"
#> [201] "GO:1903047" "GO:1903522" "GO:1903524" "GO:1903530" "GO:1903532"
#> [206] "GO:1990138" "GO:1990266" "GO:2000145" "GO:2000147" "GO:0005575"
#> [211] "GO:0005576" "GO:0110165" "GO:0001664" "GO:0003674" "GO:0005102"
#> [216] "GO:0005179" "GO:0005488" "GO:0005515" "GO:0030545" "GO:0030546"
#> [221] "GO:0031705" "GO:0031708" "GO:0048018" "GO:0071855" "GO:0098772"
#> [226] "GO:0140677"
#> 
#> $ENSG00000124207
#>  [1] "GO:0006611" "GO:0006810" "GO:0006886" "GO:0006913" "GO:0008104"
#>  [6] "GO:0008150" "GO:0009987" "GO:0015031" "GO:0033036" "GO:0045184"
#> [11] "GO:0046907" "GO:0051168" "GO:0051169" "GO:0051179" "GO:0051234"
#> [16] "GO:0051641" "GO:0051649" "GO:0005575" "GO:0005576" "GO:0005622"
#> [21] "GO:0005634" "GO:0005635" "GO:0005654" "GO:0005737" "GO:0005829"
#> [26] "GO:0012505" "GO:0016020" "GO:0031967" "GO:0031974" "GO:0031981"
#> [31] "GO:0031982" "GO:0043226" "GO:0043227" "GO:0043229" "GO:0043230"
#> [36] "GO:0043231" "GO:0043233" "GO:0065010" "GO:0070013" "GO:0070062"
#> [41] "GO:0110165" "GO:1903561" "GO:0003674" "GO:0005049" "GO:0005488"
#> [46] "GO:0005515" "GO:0019899" "GO:0031267" "GO:0051020" "GO:0140104"
#> [51] "GO:0140142"
#> 
```
