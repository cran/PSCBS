###########################################################################/**
# @set "class=AbstractCBS"
# @RdocMethod updateMeansTogether
# @alias updateMeansTogether.CBS
# @alias updateMeansTogether.PairedPSCBS
#
# @title "Updates the CN mean levels jointly in sets of segments"
#
# \description{
#  @get "title" as if they were one large segment.
#  The locus-level data is not updated/modified.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an object of the same class.
# }
#
# @author "HB"
#
# \seealso{
#   This method is utilized by @seemethod "pruneByHClust".
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("updateMeansTogether", "AbstractCBS", abstract=TRUE, private=TRUE)


###########################################################################/**
# @set "class=AbstractCBS"
# @RdocMethod hclustCNs
#
# @title "Performs a hierarchical clustering of the CN mean levels"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{size}{Argument passed to @seemethod "sampleCNs".}
#  \item{distMethod, hclustMethod}{Argument \code{method} for
#    @see "stats::dist" and "stats::hclust", respectively.}
#  \item{...}{Not used.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a \code{hclust} object as returned by @see "stats::hclust".
# }
#
# @author
#
# \seealso{
#   This method is utilized by @seemethod "pruneByHClust".
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("hclustCNs", "AbstractCBS", function(fit, size=NULL, distMethod="euclidean", hclustMethod="ward.D", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Hierarchical clustering of segmented copy numbers")
  verbose && enter(verbose, "Extracting/sampling CNs")
  C <- sampleCNs(fit, size=size, splitters=FALSE)
  verbose && str(verbose, C)

  # Drop also segments with no data points
  ok <- !is.na(C)
  ok <- rowAlls(ok, useNames=FALSE)
  C <- C[ok,,drop=FALSE]
  verbose && str(verbose, C)
  verbose && exit(verbose)

  verbose && enter(verbose, "Calculating distance matrix")
  D <- stats::dist(C, method=distMethod)
  verbose && str(verbose, D)
  verbose && exit(verbose)

  verbose && enter(verbose, "Clustering")

  # TODO: Do *weighted* hierarchical clustering
  tree <- stats::hclust(D, method=hclustMethod)
  verbose && str(verbose, tree)
  verbose && exit(verbose)

  verbose && exit(verbose)

  tree
}, private=TRUE) # hclustCNs()



###########################################################################/**
# @RdocMethod pruneByHClust
#
# @title "Prunes the CN profile by pruning and merging through hierarchical clustering"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @see "stats::cutree",
#    particularly either of thresholds \code{h} or \code{k}.}
#  \item{size, distMethod, hclustMethod}{Arguments (as well as
#    some of \code{...}) passed to @seemethod "hclustCNs".}
#  \item{merge}{If @TRUE, consecutive segments that belong to the
#    same PSCN cluster will be merged into one large segment.}
#  \item{update}{If @TRUE, segment means are updated afterwards, otherwise not.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a pruned object of the same class.
# }
#
# \examples{\dontrun{
#  fitP <- pruneByHClust(fit, h=0.25)
# }}
#
# @author
#
# @keyword internal
#*/###########################################################################
setMethodS3("pruneByHClust", "AbstractCBS", function(fit, ..., size=NULL, distMethod="euclidean", hclustMethod="ward.D", merge=TRUE, update=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Prune segments by hierarchical clustering")
  verbose && cat(verbose, "Clustering arguments:")
  verbose && str(verbose, c(list(size=size, distMethod=distMethod, hclustMethod=hclustMethod), list(...)))

  verbose && enter(verbose, "Clustering")
  tree <- hclustCNs(fit, size=size, distMethod=distMethod,
                    hclustMethod=hclustMethod, ..., verbose=less(verbose,5))
  verbose && print(verbose, tree)
  verbose && exit(verbose)

  verbose && enter(verbose, "Cutting tree")
  verbose && cat(verbose, "Cutting arguments:")
  verbose && str(verbose, c(list(tree=tree), list(...)))
  p <- cutree(tree, ...)
  verbose && str(verbose, p)

  # Group segments
  idxList <- by(names(p), p, FUN=function(x) {
    idxs <- as.integer(as.character(x))
    idxs <- sort(unique(idxs))
    list(idxs)
  })
  verbose && str(verbose, idxList)

  verbose && exit(verbose)


  # Dropping previous segment calls and quantile mean-level estimates.
  fit <- resetSegments(fit)

  verbose && enter(verbose, "Merging mean levels of clustered segments")
  fit <- updateMeansTogether(fit, idxList=idxList, verbose=less(verbose, 10))
  verbose && exit(verbose)

  if (merge) {
    verbose && enter(verbose, "Merging neighboring segments within each cluster")
    lefts <- c()
    for (ii in seq_along(idxList)) {
      verbose && enter(verbose, sprintf("Cluster #%d of %d", ii, length(idxList)))
      idxs <- idxList[[ii]]
      verbose && cat(verbose, "Segments in cluster:")
      verbose && str(verbose, idxs)

      # Indices to segments to merge
      leftsII <- idxs[which(diff(idxs) == 1L)]
      verbose && cat(verbose, "Left indices of neighboring segments:")
      verbose && str(verbose, leftsII)

      lefts <- c(lefts, leftsII)
      verbose && exit(verbose)
    } # for (ii ...)

    lefts <- sort(unique(lefts))
    verbose && cat(verbose, "Left indices of segments to be merged:")
    verbose && str(verbose, lefts)
    verbose && exit(verbose)

    verbose && enter(verbose, "Merging segments")
    lefts <- rev(lefts)
    for (ii in seq_along(lefts)) {
      fit <- mergeTwoSegments(fit, left=lefts[ii], update=FALSE)
    } # for (ii ...)
    verbose && exit(verbose)
  } # if (merge)

  if (update) {
    verbose && enter(verbose, "Updating segment means")
##  fit <- updateBoundaries(fit, verbose=less(verbose, 50))
    fit <- updateMeans(fit, verbose=less(verbose, 50))
    verbose && exit(verbose)
  }

  verbose && exit(verbose)

  fit
}, protected=TRUE) # pruneByHClust()
