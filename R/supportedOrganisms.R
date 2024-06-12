#' Supported Organisms
#' 
#' Lists which genomes and gene ids are supported for gene lengths and for gene
#' ontology
#' 
#' Goseq allows a user to provide their own bias data (usually gene lengths)
#' and/or gene categories (usually gene ontologies), but goseq also provides
#' this data automatically for many commonly used species. This function lists
#' which genome and gene ids are automatically supported by goseq. The first to
#' third columns list the genomes, gene ids, and gene id descriptions
#' respectively. The fourth column indicates whether this combination of genome
#' and id are available in the geneLengthDataBase. If a particular combination
#' is absent, goseq may still automatically fetch the gene lengths from either
#' a TxDB annotation package (must be installed) or download the data from
#' UCSC. For example gene lengths for \code{hg38} are not supported in
#' \code{geneLengthDataBase} but may be fetched by these other means. However,
#' this is not always the case. The final column indicates if GO terms will be
#' automatically fetched for the genome and id combination.  Goseq relies on an
#' \code{org} annotation package (e.g. \code{org.Hs.eg.db}) existing for the
#' organism. In general, if either the lengths or GO terms are not supported,
#' the user must enter this information manually.
#' 
#' @return A data.frame containing supported genomes and gene ids
#' 
#' @author Nadia Davidson \email{nadia.davidson@@mcri.edu.au}
#' 
#' @export
#' 
#' @examples
#' supportedOrganisms()
#' 
supportedOrganisms <- function() {
  requireNamespace("rtracklayer")

  # Start by getting the geneLenDataBase supported genomes and geneIDs
  geneIDs <- supportedGeneIDs()

  # The names of the supported geneIDs is hard-coded so we can remove things like
  # "refGene" which actually correspond to an Entrez gene ID in geneLenDataBase
  uniqueIDs <- c("knownGene", "vegaGene", "ensGene", "geneSymbol")

  geneIDs <- geneIDs[geneIDs$tablename %in% uniqueIDs, ]

  # extract the corresponding genome IDs and rearrange
  genomes <- strsplit(geneIDs$AvailableGenomes, ",")
  id <- rep(geneIDs$tablename, sapply(genomes, length))
  description <- (rep(geneIDs$GeneID, sapply(genomes, length)))
  length_supported <- rep(TRUE, length(id))
  genomes <- unlist(genomes)
  tab_temp <- data.frame(genomes, id, description, length_supported)

  ## Is there GO annotation package for that genome and gene ID?
  genome_supported_go <- gsub("[0-9]+", "", genomes) %in% names(.ORG_PACKAGES)
  id_supported_go <- id %in% names(.ID_MAP)

  go_supported <- genome_supported_go & id_supported_go
  tab_temp <- data.frame(tab_temp, go_supported)
  colnames(tab_temp) <- c(
    "Genome", "Id", "Id Description", "Lengths in geneLeneDataBase",
    "GO Annotation Available"
  )

  ### Add in the newer genomes which aren't supported for length, but are for GO terms
  base <- unfactor(ucscGenomes())$db
  other_genomes <- base[!base %in% genomes]
  new_genomes <- other_genomes[gsub("[0-9]+", "", other_genomes) %in% names(.ORG_PACKAGES)]
  new_genomes <- c(new_genomes, names(.ORG_PACKAGES)[!names(.ORG_PACKAGES) %in% gsub("[0-9]+", "", c(new_genomes, genomes))])
  tab_new <- data.frame(new_genomes, "", "", FALSE, TRUE)
  colnames(tab_new) <- c(
    "Genome", "Id", "Id Description", "Lengths in geneLeneDataBase",
    "GO Annotation Available"
  )
  tab_supported <- rbind(tab_temp, tab_new)
  tab_supported[order(tab_supported$Genome), ]
}
