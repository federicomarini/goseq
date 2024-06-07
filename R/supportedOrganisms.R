#############################################################################
# Description: Lists which genomes and gene IDs are supported by goseq (length and GO terms)
# Notes:
# Author: Nadia Davidson
# Date Modified: 2016/10/04

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
