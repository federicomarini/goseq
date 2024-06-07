#############################################################################
# Description: Fetches gene length data for the genome and id specified.
# Notes:
# Author: Matthew Young, Nadia Davidson
# Date Modified: 7/4/2015



#' Retrieves Gene length data
#' 
#' Gets the length of each gene in a vector.
#' 
#' Length data is obtained from data obtained from the UCSC genome browser for
#' each combination of \code{genome} and \code{id}.  As fetching this data at
#' runtime is time consuming, a local copy of the length information for common
#' genomes and gene ID are included in the \pkg{geneLenDataBase} package.  This
#' function uses this package to fetch the required data.
#' 
#' The length of a gene is taken to be the median length of all its mature,
#' mRNA, transcripts.  It is always preferable to obtain length information
#' directly for the gene ID used to summarize your count data, rather than
#' converting IDs and then using the supplied databases.  Even when two genes
#' have a one-to-one mapping between different identifier conventions (which is
#' often not the case), they frequently refer to slightly different regions of
#' the genome with different lengths.  It is therefore recommended that the
#' user perform the full analysis in terms of only one gene ID, or manually
#' obtain their own length data for the identifier used to bin reads by gene.
#' 
#' @param genes A vector or list of the genes for which length information is
#' required.
#' @param genome A string identifying the genome that \code{genes} refer to.
#' For a list of supported organisms run \code{\link{supportedGenomes}}.
#' @param id A string identifying the gene identifier used by \code{genes}.
#' For a list of supported gene IDs run \code{\link{supportedGeneIDs}}.
#' @return Returns a vector of the gene lengths, in the same order as
#' \code{genes}.  If length data is unavailable for a particular gene NA is
#' returned in that position.  The returned vector is intended for use with the
#' \code{bias.data} option of the \code{\link{nullp}} function.
#' 
#' @author Matthew D. Young \email{myoung@@wehi.edu.au}
#' 
#' @seealso \code{\link{supportedGenomes}}, \code{\link{supportedGeneIDs}},
#' \code{\link{nullp}}, \pkg{geneLenDataBase}
#' 
#' @export
#' 
#' @examples
#' 
#' genes <- c("ENSG00000124208", 
#'            "ENSG00000182463", 
#'            "ENSG00000124201", 
#'            "ENSG00000124205", 
#'            "ENSG00000124207")
#' getlength(genes,'hg19','ensGene')
#' 
getlength <- function(genes, genome, id) {
  # geneLenDataBase contains all the gene length data, check it is installed
  if (!any(library()$results[, "Package"] == "geneLenDataBase")) {
    stop("Package geneLenDataBase is required for automatic gene length lookup.")
  }
  # Get a list of the length data files that are in the package
  repo <- grep(".*\\..*\\.LENGTH", as.data.frame(data(package = "geneLenDataBase")$results,
    stringsAsFactors = FALSE
  )$Item, ignore.case = TRUE, value = TRUE)
  # Is this genome/id in the local repository?
  if (paste(genome, id, "LENGTH", sep = ".") %in% repo) {
    # Yes it is, so load it and group transcript lengths into genes
    message("Loading ", genome, " length data...")
    data(list = paste(genome, id, "LENGTH", sep = "."), package = "geneLenDataBase")
    data <- split(get(paste(genome, id, "LENGTH", sep = "."))$Length, get(paste(genome, id, "LENGTH", sep = "."))$Gene)

    # Get the names of the genes used in the genome package
    gene_names <- get(paste(genome, id, "LENGTH", sep = "."))$Gene
  } else { ### Try to find the gene lengths from TXDB annotation packages
    message("Can't find ", genome, "/", id, " length data in genLenDataBase...")

    ## Is there a TxDB package installed for the organism?
    installedPackages <- rownames(installed.packages())

    # first try looking for a match with the same type of gene IDs
    txdbPattern <- paste("TxDb", "*", genome, id, sep = ".")
    txdbPack <- installedPackages[grep(txdbPattern, installedPackages)]

    # if not look for the package with generic gene IDs and convert later
    if (length(txdbPack) == 0) {
      orgstring <- as.character(.ORG_PACKAGES[match(gsub("[0-9]+", "", genome), names(.ORG_PACKAGES))])
      packageExtension <- tail(strsplit(orgstring, ".", fixed = TRUE)[[1]], n = 1)
      genericID <- names(.ID_MAP[which(.ID_MAP == packageExtension)])
      txdbPattern <- paste("TxDb", "*", genome, genericID, sep = ".")
      txdbPack <- installedPackages[grep(txdbPattern, installedPackages)]
      if (length(txdbPack) > 0) {
        library(paste(orgstring, "db", sep = "."), character.only = TRUE)
        userid <- as.character(.ID_MAP[match(id, names(.ID_MAP))])
        user2core <- toTable(get(paste(orgstring, userid, sep = "")))
      }
    }
    if (length(txdbPack) > 0) {
      message("Found the annotation package, ", txdbPack)
      message("Trying to get the gene lengths from it.")
      library(txdbPack, character.only = TRUE)
      txdb <- get(txdbPack)
      txlens <- transcriptLengths(txdb)
      core2len <- split(txlens$tx_len, txlens$gene_id)
      if (length(grep(id, txdbPattern)) != 0) {
        data <- core2len
        gene_names <- names(data)
      } else { # We need to conver between gene IDs. Just use Matt's code from getgo.R
        # Now we need to replicate the core IDs that need replicating
        core2len <- core2len[match(user2core[, 1], names(core2len))]
        # Now we can replace the labels on this list with the user ones from user2core,
        # but there will be duplicates, so we have to unlist, label, then relist
        data <- split(unlist(core2len, FALSE, FALSE), rep(user2core[, 2], sapply(core2len, length)))
        gene_names <- names(data)
      }
    } else { ## Can't use TxDB
      if (length(txdbPack) == 0 & any(.TXDB_ORGS == genome)) {
        message("A TxDb annotation package exists for ", genome, ". Consider installing it to get the gene lengths.", sep = "")
      }
      # If it's not in geneLenDataBase, try to load the data from UCSC using rtracklayer
      if (!requireNamespace("rtracklayer")) {
        stop("Can't load the rtracklayer library. Is it installed?...")
      }
      message(
        "Trying to download from UCSC. ",
        "This might take a couple of minutes. "
      )

      tryCatch(
        {
          # If not, then the user will have to supply it themselves
          # download the relevant table from UCSC
          session <- browserSession()
          genome(session) <- genome
          if (id == "geneSymbol") { # take the geneSymbols from refSeq genes
            query <- ucscTableQuery(session, "refGene")
          } else {
            query <- ucscTableQuery(session, id)
          }
          transTable <- getTable(query)

          # all this below is to get the length information of each gene
          starts <- strsplit(as.character(transTable$exonStart), ",")
          ends <- strsplit(as.character(transTable$exonEnd), ",")
          exons <- mapply(data.frame, starts, ends, SIMPLIFY = FALSE)
          get_length <- function(x) {
            sum(as.integer(as.character(x[, 2])) - as.integer(as.character(x[, 1])) + 1)
          }
          lengths <- sapply(exons, get_length)

          # now we need to get the correct list of gene names.
          names2 <- as.character(transTable$name2)
          gene_names <- names2
          if (!any(toupper(genes) %in% toupper(names2))) {
            # check which column has the gene names
            names1 <- as.character(transTable$name1)
            if (any(toupper(genes) %in% toupper(names1))) {
              gene_names <- names1
            }
          }
          data <- split(lengths, gene_names)
        },
        error = function(e) {
          stop(
            "Length information for genome ", genome, " and gene ID ", id,
            " is not available. You will have to specify bias.data manually."
          )
        }
      )
    }
  }
  # In order to assign just one length to each gene, we store the Median, Min Max
  # and count of transcripts for each gene for later use
  len <- data.frame(
    Gene = names(data), Median = sapply(data, median), Min = sapply(data, min),
    Max = sapply(data, max), Count = sapply(data, length)
  )

  # How many of the genes match the names the user specified with "genes"
  matched_frac <- sum(toupper(genes) %in% toupper(gene_names)) / length(genes)
  # None? Throw an error.
  if (matched_frac == 0) {
    stop(
      "The gene names specified do not match the gene names for genome ",
      genome, " and ID ", id, ".\n\tGene names given were: ", paste(genes[1:10], collapse = ", "),
      "\n\tRequired gene names are: ", paste(gene_names[1:10], collapse = ", ")
    )
  }
  # Less than 60%? Give a warning and continue
  if (matched_frac < .6) {
    warning(
      "More than 40% of gene names specified did not match the gene names for genome ",
      genome, " and ID ", id, ".  No length data will be available for these genes.\n\tGene names which failed to match were: ",
      paste(genes[!genes %in% gene_names][1:10], collapse = ", "), "\n\tRequired gene names are: ",
      paste(gene_names[1:10], collapse = ", ")
    )
  }
  # Finally, return the median transcript length for each matching gene
  return(len$Median[match(toupper(genes), toupper(len$Gene))])
}
