#' Fetch GO categories
#' 
#' Obtains all gene ontology (GO) categories associated with a set of genes
#' using the relevant organism package.
#' 
#' This function attempts to make use of the organism packages
#' (org.<Genome>.<GeneID>.db) to obtain the mapping between gene ID and GO
#' categories.  As with \code{\link{getlength}} it is preferable that the same
#' gene identifier system is used for both summarization and retrieving GO
#' categories.
#' 
#' Valid options for the \code{fetch.cats} argument are any combination of
#' "GO:CC", "GO:BP", "GO:MF" & "KEGG".  The three GO terms refer to the
#' Cellular Component, Biological Process and Molecular Function respectively.
#' "KEGG" refers to KEGG pathways.
#' 
#' Note that \code{getgo} is a convenience function, designed to make
#' extracting mappings between GO categories and Gene ID easy.  For less common
#' organisms and/or gene ID \code{getgo} may fail to return a mapping even when
#' a legitimate mapping exists in the relevant organism package.  If
#' \code{getgo} fails, you should always try to build the mapping yourself from
#' the organism package (if one exists) before deciding that the information is
#' unavailable.  Further information and examples of this can be found in the
#' package Vignette.
#' 
#' @param genes A vector or list of genes to get the associated GO categories.
#' @param genome A string identifying the genome that \code{genes} refer to.
#' For a list of supported organisms run \code{\link{supportedGenomes}}.
#' @param id A string identifying the gene identifier used by \code{genes}.
#' For a list of supported gene IDs run \code{\link{supportedGeneIDs}}.
#' @param fetch.cats A vector specifying which categories to fetch the mapping
#' between category names and genes for.  See details for valid options.
#' @return A list where each entry is named by a gene and contains a vector of
#' all the associated GO categories.  This can be used directly with the
#' \code{gene2cat} option in \code{\link{goseq}}.
#' 
#' @author Matthew D. Young \email{myoung@@wehi.edu.au}
#' 
#' @seealso \code{\link{supportedGenomes}}, \code{\link{supportedGeneIDs}},
#' \code{\link{goseq}}
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
#' getgo(genes,'hg19','ensGene')
#' 
getgo <- function(genes, genome, id, fetch.cats = c("GO:CC", "GO:BP", "GO:MF")) {
  # Check for valid input
  if (any(!fetch.cats %in% c("GO:CC", "GO:BP", "GO:MF", "KEGG"))) {
    stop("Invaled category specified.  Categories can only be GO:CC, GO:BP, GO:MF or KEGG")
  }
  # Convert from genome ID to org.__.__.db format
  orgstring <- as.character(.ORG_PACKAGES[match(gsub("[0-9]+", "", genome), names(.ORG_PACKAGES))])
  # Multimatch or no match
  if (length(orgstring) != 1) {
    stop("Couldn't grab GO categories automatically.  Please manually specify.")
  }
  # Load the library
  library(paste(orgstring, "db", sep = "."), character.only = TRUE)
  # What is the default ID that the organism package uses?
  coreid <- strsplit(orgstring, "\\.")[[1]][3]

  # Now we need to convert it into the naming convention used by the organism packages
  userid <- as.character(.ID_MAP[match(id, names(.ID_MAP))])
  # Multimatch or no match
  if (is.na(userid) | (length(userid) != 1)) {
    stop("Couldn't grab GO categories automatically.  Please manually specify.")
  }
  # The (now loaded) organism package contains a mapping between the internal ID and whatever
  # the default is (usually eg), the rest of this function is about changing that mapping to
  # point from categories to the ID specified
  # Fetch the mapping in its current format
  # Because GO is a directed graph, we need to get not just the genes associated with each ID,
  # but also those associated with its children.  GO2ALLEGS does this.
  core2cat <- NULL
  if (length(grep("^GO", fetch.cats)) != 0) {
    # Get the name of the function which maps gene ids to go terms
    # usually this will be "GO2ALLEG"
    gomapFunction <- .ORG_GOMAP_FUNCTION[orgstring]
    if (is.na(gomapFunction)) gomapFunction <- .ORG_GOMAP_FUNCTION["default"]
    x <- toTable(get(paste(orgstring, gomapFunction, sep = "")))
    # Keep only those ones that we specified and keep only the names
    # 		core2cat=x[x$Ontology%in%gsub("^GO:",'',fetch.cats),1:2]
    x[!x$Ontology %in% gsub("^GO:", "", fetch.cats), 2] <- "Other"
    core2cat <- x[, 1:2]
    colnames(core2cat) <- c("gene_id", "category")
  }
  if (length(grep("^KEGG", fetch.cats)) != 0) {
    x <- toTable(get(paste(orgstring, "PATH", sep = "")))
    # Either add it to existing table or create a new one
    colnames(x) <- c("gene_id", "category")
    if (!is.null(core2cat)) {
      core2cat <- rbind(core2cat, x)
    } else {
      core2cat <- x
    }
  }

  # Now we MAY have to convert the "gene_id" column to the format we are using
  if (coreid != userid) {
    # The mapping between user id and core id, don't use the <USER_ID>2<CORE_ID> object as the naming is not always consistent
    user2core <- toTable(get(paste(orgstring, userid, sep = "")))
    # Throw away any user ID that doesn't appear in core2cat
    user2core <- user2core[user2core[, 1] %in% core2cat[, 1], ]
    # Make a list version of core2cat, we'll need it
    list_core2cat <- split(core2cat[, 2], core2cat[, 1])
    # Now we need to replicate the core IDs that need replicating
    list_core2cat <- list_core2cat[match(user2core[, 1], names(list_core2cat))]
    # Now we can replace the labels on this list with the user ones from user2core,
    # but there will be duplicates, so we have to unlist, label, then relist
    user2cat <- split(unlist(list_core2cat, FALSE, FALSE), rep(user2core[, 2], sapply(list_core2cat, length)))
    # Now we only want each category listed once for each entry...
    user2cat <- sapply(user2cat, unique)
    ### In case you don't believe that this works as it should, here is the slow as all hell way for comparison...
    ### Make first list
    ## list_user2core=split(user2core[,1],user2core[,2])
    ### Make the second
    ## list_core2cat=split(core2cat[,2],core2cat[,1])
    ### Go through each entry in first list and expand using second...
    ## user2cat=sapply(list_user2core,function(u){unique(unlist(list_core2cat[u],FALSE,FALSE))})
  } else {
    # We don't need to convert anything (WOO!), so just make it into a list
    user2cat <- split(core2cat[, 2], core2cat[, 1])
    user2cat <- sapply(user2cat, unique)
  }
  # remove any empty strings
  user2cat <- lapply(user2cat, function(x) {
    if (length(x) > 1) x <- x[x != "Other"]
    x
  })

  ## we don't like case sensitivity
  names(user2cat) <- toupper(names(user2cat))

  # Now look them up
  return(user2cat[toupper(genes)])
}
