############################################################################/**
# @RdocDefault installDNAcopy
#
# @title "Install the DNAcopy package"
#
# @synopsis
#
# \description{
#   @get "title", if missing.
# }
#
# \arguments{
#   \item{...}{Arguments passed to the install function.}
#   \item{force}{If @FALSE and the \pkg{DNAcopy} package is already
#     installed, then it will not be re-install.
#     If @TRUE, it will be installed.}
# }
#
# \value{
#   Returns nothing.
# }
#
# \details{
#   This function is will download and call the \code{biocLite()}
#   installation function from the Bioconductor Project website.
#   This function will also make sure that \pkg{DNAcopy} is loaded so
#   that it is reported by @see "utils::sessionInfo".
# }
#
# @author "HB"
#
# @keyword internal
#*/############################################################################
setMethodS3("installDNAcopy", "default", function(..., force=FALSE) {
  # Argument 'force':
  force <- Arguments$getLogical(force)

  # Package to be installed
  pkgName <- "DNAcopy"

  # Is DNAcopy already available?
  if (!force && isPackageInstalled(pkgName)) {
    library(pkgName, character.only=TRUE)
    return(invisible())
  }

  # If not, install it...
  if (getRversion() >= "3.5.0") {
    if (!isPackageInstalled("BiocManager")) install.packages("BiocManager")
    ns <- getNamespace("BiocManager")
    BiocManager_install <- get("install", envir = ns)
    BiocManager_install(pkgName, ...)
  } else {
    # To please R CMD check
    biocLite <- NULL; rm(list="biocLite")
    source("http://www.bioconductor.org/biocLite.R")
    biocLite(pkgName, ...)
  }
  
  # ...and load it
  library(pkgName, character.only=TRUE)

  return(invisible())
}) # installDNAcopy()
