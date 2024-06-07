#  ZZZ.R

.onAttach <- function(libname, pkgname) {
  if (.Platform$OS.type == "windows" && .Platform$GUI == "Rgui") {
    winMenuAddItem("Vignettes", "goseq", "shell.exec(system.file(\"doc\",\"goseq.pdf\",package=\"goseq\"))")
  }
}

# These two variables are required for automatic fetching of categories to function.  Their purpose is to take the UCSC genome and gene ID values given when looking up length data and convert them to the names used for the same organism and gene identifier in the organism packages.

# Mappings that are primarily required by getgo, the purpose of this is to convert the UCSC genome IDs, to the bioconductor organism names, e.g. "mm"->"org.Mm."
.ORG_PACKAGES <- paste("org.", c("Ag.eg", "At.tair", "Bt.eg", "Ce.eg", "Cf.eg", "Dm.eg", "Dr.eg", "EcK12.eg", "EcSakai.eg", "Gg.eg", "Hs.eg", "Mm.eg", "Mmu.eg", "Pf.plasmo", "Pt.eg", "Rn.eg", "Sc.sgd", "Ss.eg", "Xl.eg"), sep = "")
names(.ORG_PACKAGES) <- c("anoGam", "Arabidopsis", "bosTau", "ce", "canFam", "dm", "danRer", "E. coli K12", "E. coli Sakai", "galGal", "hg", "mm", "rheMac", "Malaria", "panTro", "rn", "sacCer", "susScr", "xenTro")

# These are the only formats supported by getgo at the moment, the purpose is to convert the USCC gene
# ID formats, to the shorthand used by the bioconductor organism packages, .e.g. "refGene"->"ENSEMBL"
.ID_MAP <- c("eg", "eg", "ENSEMBL", "SYMBOL", "sgd", "plasmo", "tair")
names(.ID_MAP) <- c("knownGene", "refGene", "ensGene", "geneSymbol", "sgd", "plasmo", "tair")

# Below are the exceptions to the function name for gene to go term mappings
.ORG_GOMAP_FUNCTION <- c("GO2ALLEGS", "GO2ALLTAIRS", "GO2ALLORFS", "GO2ALLORFS")
names(.ORG_GOMAP_FUNCTION) <- c("default", "org.At.tair", "org.Pf.plasmo", "org.Sc.sgd")

# TxDb Length databases
.TXDB_ORGS <- c("ce6", "dm3", "hg18", "hg19", "hg38", "mm10", "mm9", "rn4", "rn5", "sacCer2", "sacCer3")
