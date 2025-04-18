###########################################################################/**
# @RdocDefault segmentByCBS
# @alias segmentByCBS.data.frame
# @alias segmentByCBS.CBS
# @alias segmentByCBS.CNA
# @alias segmentByCBS
#
# @title "Segment genomic signals using the CBS method"
#
# \description{
#  @get "title" of the \pkg{DNAcopy} package.
#  This is a convenient low-level wrapper for the \code{DNAcopy::segment()}
#  method.  It is intended to be applied to a sample at the time.
#  For more details on the Circular Binary Segmentation (CBS) method
#  see [1,2].
# }
#
# @synopsis
#
# \arguments{
#   \item{y}{A @numeric @vector of J genomic signals to be segmented.}
#   \item{chromosome}{Optional @numeric @vector of length J, specifying
#       the chromosome of each loci.  If a scalar, it is expanded to
#       a vector of length J.}
#   \item{x}{Optional @numeric @vector of J genomic locations.
#            If @NULL, index locations \code{1:J} are used.}
#   \item{index}{An optional @integer @vector of length J specifying
#     the genomewide indices of the loci.}
#   \item{w}{Optional @numeric @vector in [0,1] of J weights.}
#   \item{undo}{A non-negative @numeric.  If greater than zero, then
#       arguments \code{undo.splits="sdundo"} and \code{undo.SD=undo}
#       are passed to \code{DNAcopy::segment()}.
#       In the special case when \code{undo} is +@Inf, the segmentation
#       result will not contain any changepoints (in addition to what
#       is specified by argument \code{knownSegments}).}
#   \item{avg}{A @character string specifying how to calculating
#         segment mean levels \emph{after} change points have been
#         identified.}
#   \item{...}{Additional arguments passed to the \code{DNAcopy::segment()}
#       segmentation function.}
#   \item{joinSegments}{If @TRUE, there are no gaps between neighboring
#     segments.
#     If @FALSE, the boundaries of a segment are defined by the support
#     that the loci in the segments provides, i.e. there exist a locus
#     at each end point of each segment.  This also means that there
#     is a gap between any neighboring segments, unless the change point
#     is in the middle of multiple loci with the same position.
#     The latter is what \code{DNAcopy::segment()} returns.
#   }
#   \item{knownSegments}{Optional @data.frame specifying
#     \emph{non-overlapping} known segments.  These segments must
#     not share loci.  See @see "findLargeGaps" and @see "gapsToSegments".}
#   \item{seed}{An (optional) @integer specifying the random seed to be
#     set before calling the segmentation method.  The random seed is
#     set to its original state when exiting.  If @NULL, it is not set.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @see "CBS" object.
# }
#
# \details{
#   Internally @see "DNAcopy::segment" of \pkg{DNAcopy} is used to
#   segment the signals.
#   This segmentation method support weighted segmentation.
# }
#
# \section{Reproducibility}{
#   The \code{DNAcopy::segment()} implementation of CBS uses approximation
#   through random sampling for some estimates.  Because of this,
#   repeated calls using the same signals may result in slightly
#   different results, unless the random seed is set/fixed.
# }
#
# \section{Missing and non-finite values}{
#   Signals may contain missing values (@NA or @NaN), but not
#   infinite values (+/-@Inf).  Loci with missing-value signals
#   are preserved and keep in the result.
#
#   Likewise, genomic positions may contain missing values.
#   However, if they do, such loci are silently excluded before
#   performing the segmentation, and are not kept in the results.
#   The mapping between the input locus-level data and ditto of
#   the result can be inferred from the \code{index} column of
#   the locus-level data of the result.
#
#   None of the input data may have infinite values,
#   i.e. -@Inf or +@Inf. If so, an informative error is thrown.
# }
#
# \examples{
#   @include "../incl/segmentByCBS.Rex"
#   @include "../incl/segmentByCBS,plot.Rex"
#   @include "../incl/segmentByCBS,tests.Rex"
# }
#
# @author "HB"
#
# \references{
#  [1] @include "../incl/OlshenVenkatraman_2004.Rd" \cr
#  [2] @include "../incl/VenkatramanOlshen_2007.Rd" \cr
# }
#
# \seealso{
#   To segment allele-specific tumor copy-number signals from a tumor
#   \emph{with} a matched normal, see @see "segmentByPairedPSCBS".
#   For the same \emph{without} a matched normal,
#   see @see "segmentByNonPairedPSCBS".
#
#   It is also possible to prune change points after segmentation (with
#   identical results) using
#   \code{\link[PSCBS:pruneBySdUndo.CBS]{pruneBySdUndo}()}.
# }
# @keyword IO
#*/###########################################################################
setMethodS3("segmentByCBS", "default", function(y, chromosome=0L, x=NULL, index=seq_along(y), w=NULL, undo=0, avg=c("mean", "median"), ..., joinSegments=TRUE, knownSegments=NULL, seed=NULL, verbose=FALSE) {
  R_SANITY_CHECK <- TRUE

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # DNAcopy::getbdry() is slow for now default settings.  Below we
  # implement a memoized version of this function.
  getbdry2 <- function(eta, nperm, alpha, tol=0.01, verbose=FALSE) {
    # Explictly setup cache root here, since it's only done by 'R.cache'
    # if that package is attached.  Here we only load it. /HB 2013-09-27
    .setupCacheRootPath()

    key <- list(method="segmentByCBS",
                eta=eta, nperm=as.integer(nperm), alpha=alpha, tol=tol,
                version="0.16.1")
    dirs <- c("PSCBS", "segmentByCBS", "sbdry")
    bdry <- loadCache(key=key, dirs=dirs)
    if (!is.null(bdry)) return(bdry)

    max.ones <- floor(nperm * alpha) + 1L
    bdry <- getbdry(eta=eta, nperm=nperm, max.ones=max.ones, tol=tol)

    saveCache(bdry, key=key, dirs=dirs)

    bdry
  } # getbdry2()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'y':
  disallow <- c("Inf")
  y <- Arguments$getDoubles(y, disallow=disallow)
  nbrOfLoci <- length(y)

  length2 <- rep(nbrOfLoci, times=2)

  # Argument 'chromosome':
  if (is.null(chromosome)) {
    chromosome <- 0L
  } else {
    disallow <- c("Inf")
    chromosome <- Arguments$getIntegers(chromosome, range=c(0,Inf), disallow=disallow)
    if (length(chromosome) > 1) {
      chromosome <- Arguments$getIntegers(chromosome, length=length2, disallow=disallow)
  ##    # If 'chromosome' is a vector of length J, then it must contain
  ##    # a unique chromosome.
  ##    chromosomes <- sort(unique(chromosome))
  ##    if (length(chromosomes) > 1) {
  ##      stop("Argument 'chromosome' specifies more than one unique chromosome: ", paste(seqToHumanReadable(chromosomes), collapse=", "))
  ##    }
  ##    chromosome <- chromosomes
    }
  }

  # For future usage
  chrom <- rep(chromosome, length.out=nbrOfLoci)

  # Argument 'x':
  if (is.null(x)) {
    x <- seq_len(nbrOfLoci)
  } else {
    disallow <- c("Inf")
    x <- Arguments$getDoubles(x, length=length2, disallow=disallow)
  }

  # Argument 'index':
  if (is.null(index)) {
    index <- seq_along(y)
  } else {
    index <- Arguments$getIndices(index)
  }

  # Argument 'w':
  hasWeights <- !is.null(w)
  if (hasWeights) {
    disallow <- c("NA", "NaN", "Inf")
    w <- Arguments$getDoubles(w, range=c(0,1), length=length2, disallow=disallow)
  }

  # Argument 'undo':
  undo <- Arguments$getDouble(undo, range=c(0,Inf))

  # Argument 'avg':
  avg <- match.arg(avg)

  # Argument 'cpFlavor':
  joinSegments <- Arguments$getLogical(joinSegments)

  # Argument 'knownSegments':
  if (is.null(knownSegments)) {
    knownSegments <- data.frame(chromosome=integer(0), start=integer(0), end=integer(0))
  } else {
#    if (!joinSegments) {
#      stop("Argument 'knownSegments' should only be specified if argument 'joinSegments' is TRUE.")
#    }
  }

  if (!is.data.frame(knownSegments)) {
    stop("Argument 'knownSegments' is not a data.frame: ", class(knownSegments)[1])
  }

  if (!all(is.element(c("chromosome", "start", "end"), colnames(knownSegments)))) {
    stop("Argument 'knownSegments' does not have the required column names: ", hpaste(colnames(knownSegments)))
  }

  # Detailed validation of 'knownSegments'.
  for (chr in sort(unique(knownSegments$chromosome))) {
    dd <- subset(knownSegments, chromosome == chr)

    # Order segments by 'start'.
    o <- order(dd$start)
    dd <- dd[o,]

    # Known segments must not share 'start' or 'end' loci
    for (field in c("start", "end")) {
      xs <- dd[[field]]
      xs <- xs[!is.na(xs)]
      if (anyDuplicated(xs) > 0) {
        print(knownSegments)
        stop(sprintf("Detected segments on chromosome %s with non-unique '%s' positions in argument 'knownSegments'", chr, field))
      }
    } # for (field ...)

    # Known segments must not overlap
    if (!all(dd$start[-1] >= dd$end[-nrow(dd)], na.rm=TRUE)) {
      print(knownSegments)
      stop("Detected overlapping segments on chromosome ", chr, " in argument 'knownSegments'.")
    }
  }

  # Argument 'seed':
  if (!is.null(seed)) {
    seed <- Arguments$getIntegers(seed)
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Segmenting by CBS")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setup up data", level=-10)
  data <- data.frame(chrom=chrom, x=x, y=y, index=index)
  if (hasWeights) {
    verbose && cat(verbose, "Adding locus-specific weights", level=-10)
    data$w <- w
  }
  verbose && str(verbose, data, level=-10)
  # Not needed anymore
  chrom <- x <- index <- y <- w <- NULL
  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Drop data points without known genomic positions, because that
  # is what DNAcopy::CNA() will do otherwise.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ok <- (!is.na(data$chrom) & !is.na(data$x))
  if (any(!ok)) {
    verbose && enter(verbose, "Dropping loci with unknown locations", level=-10)
    verbose && cat(verbose, "Number of loci dropped: ", sum(!ok), level=-10)
    data <- data[ok,,drop=FALSE]
    verbose && exit(verbose)
  }
  ok <- NULL # Not needed anymore


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reorder data points along the genome, because that is what
  # DNAcopy::segment() will return.  At the end, we will undo
  # the sort such that the returned 'data' object is always in
  # the same order and number of loci as the input data.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Ordering data along genome", level=-50)
  o <- order(data$chrom, data$x, decreasing=FALSE, na.last=TRUE)
  # Any change?
  if (any(o != seq_along(o))) {
    data <- data[o,,drop=FALSE]
  }
  o <- NULL # Not needed anymore
  verbose && str(verbose, data, level=-50)
  verbose && exit(verbose)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Multiple chromosomes?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify all chromosomes, excluding missing values
  chromosomes <- sort(unique(data$chrom), na.last=NA)
  nbrOfChromosomes <- length(chromosomes)
  if (nbrOfChromosomes > 1) {
    verbose && enter(verbose, "Segmenting multiple chromosomes")
    verbose && cat(verbose, "Number of chromosomes: ", nbrOfChromosomes)

    # Generate random seeds?
    seeds <- NULL
    if (!is.null(seed)) {
      randomSeed("set", seed=seed, kind="L'Ecuyer-CMRG")
      on.exit(randomSeed("reset"), add=TRUE)
      verbose && printf(verbose, "Random seed temporarily set (seed=c(%s), kind=\"L'Ecuyer-CMRG\")\n", paste(seed, collapse=", "))
      seeds <- randomSeed("advance", n=nbrOfChromosomes)
      verbose && printf(verbose, "Produced %d seeds from this stream for future usage\n", length(seeds))
    }

    fitList <- list()
    for (kk in seq_len(nbrOfChromosomes)) {
      chromosomeKK <- chromosomes[kk]
      chrTag <- sprintf("Chr%02d", chromosomeKK)
      verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", kk, chrTag, nbrOfChromosomes))

      seedKK <- seeds[[kk]]

      # Extract subset of data and parameters for this chromosome
      dataKK <- subset(data, chrom == chromosomeKK)
      verbose && str(verbose, dataKK, level=-10)
      chrom <- x <- index <- y <- w <- NULL
      fields <- attachLocally(dataKK, fields=c("chrom", "x", "index", "y", "w"))
      dataKK <- NULL # Not needed anymore

      knownSegmentsKK <- NULL
      if (!is.null(knownSegments)) {
        knownSegmentsKK <- subset(knownSegments, chromosome == chromosomeKK)
        if (nrow(knownSegmentsKK) == 0L) {
          knownSegmentsKK <- data.frame(chromosome=chromosomeKK, start=-Inf, end=+Inf)
        }
        verbose && cat(verbose, "Known segments:", level=-5)
        verbose && print(verbose, knownSegmentsKK, level=-5)
      }

      fitList[[chrTag]] <- future({
        fit <- segmentByCBS(y=y,
                  chromosome=chrom, x=x,
                  w=w,
                  index=index,
                  undo=undo,
                  avg=avg,
                  joinSegments=joinSegments,
                  knownSegments=knownSegmentsKK,
                  ...,
                  seed=seedKK,
                  verbose=verbose)

        # Sanity checks
        if (R_SANITY_CHECK) {
          if (nrow(knownSegmentsKK) == 0) {
            # Since all missing data have been dropped...
            .stop_if_not(nrow(fit$data) == length(y))
            # ...and ordered along the genome already.
            .stop_if_not(all.equal(fit$data$y, y))
          }

          # Assert weights were used
          .stop_if_not(!hasWeights || !is.null(fit$data$w))
        } # if (R_SANITY_CHECK)

        verbose && print(verbose, head(as.data.frame(fit)), level=-10)
        verbose && print(verbose, tail(as.data.frame(fit)), level=-10)

        fit
      }, seed = TRUE, label = sprintf("segmentByCBS-%s", chrTag))

      rm(list=fields) # Not needed anymore
      verbose && exit(verbose)
    } # for (kk ...)

    verbose && enter(verbose, "Merging (independently) segmented chromosome", level=-50)
    fitList <- value(fitList)
    ## former Reduce() w/ append(..., addSplit = TRUE)
    fit <- do.call(c, args = c(fitList, addSplit = TRUE))
    # Not needed anymore
    fitList <- NULL

    # Update parameters that otherwise may be incorrect
    fit$params$seed <- seed

    verbose && str(verbose, fit, level=-10)
    verbose && exit(verbose)

    segs <- as.data.frame(fit)
    if (nrow(segs) < 6) {
      verbose && print(verbose, segs, level=-10)
    } else {
      verbose && print(verbose, head(segs), level=-10)
      verbose && print(verbose, tail(segs), level=-10)
    }

    verbose && exit(verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Return results
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return(fit)
  } # if (nbrOfChromosomes > 1)

  verbose && cat(verbose, "Chromosome: ", data$chrom[1L])

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset 'knownSegments'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Keeping only current chromosome for 'knownSegments'", level=-10)

  # Assume no missing values
  currChromosome <- data$chrom[1]
  verbose && cat(verbose, "Chromosome: ", currChromosome, level=-10)

  knownSegments <- subset(knownSegments, chromosome == currChromosome)
  if (nrow(knownSegments) == 0L) {
    knownSegments <- data.frame(chromosome=currChromosome, start=-Inf, end=+Inf)
  }
  nbrOfSegments <- nrow(knownSegments)

  verbose && cat(verbose, "Known segments for this chromosome:", level=-10)
  verbose && print(verbose, knownSegments, level=-10)

  verbose && exit(verbose)


  # Sanity checks
  if (R_SANITY_CHECK) {
    # Here 'knownSegments' should specify at most a single chromosome
    uChromosomes <- sort(unique(knownSegments$chromosome))
    if (length(uChromosomes) > 1) {
      stop("INTERNAL ERROR: Argument 'knownSegments' specifies more than one chromosome: ", hpaste(uChromosomes))
    }
  } # if (R_SANITY_CHECK)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Multiple segments?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity check of limitation  /HB 2011-10-19
  if (nbrOfSegments > 1) {
    verbose && enter(verbose, "Segmenting multiple segments on current chromosome", level=-5)
    verbose && cat(verbose, "Number of segments: ", nbrOfSegments, level=-5)

    # Create a splitter-only CBS object
    dataS <- data.frame(y=c(0,0), chromosome=c(1,2), x=c(0,0))
    if (hasWeights) dataS$w <- 1
    splitter <- segmentByCBS(dataS)
    dataS <- NULL
    suppressWarnings({
      splitter <- extractSegment(splitter, 2)
      # Sanity check
      if (R_SANITY_CHECK) {
        .stop_if_not(nbrOfSegments(splitter, splitters=TRUE) == 1)
      } # if (R_SANITY_CHECK)
    })


    # Generate random seeds?
    seeds <- NULL
    if (!is.null(seed)) {
      randomSeed("set", seed=seed, kind="L'Ecuyer-CMRG")
      on.exit(randomSeed("reset"), add=TRUE)
      verbose && printf(verbose, "Random seed temporarily set (seed=c(%s), kind=\"L'Ecuyer-CMRG\")\n", paste(seed, collapse=", "))
      seeds <- randomSeed("advance", n=nbrOfSegments)
      verbose && printf(verbose, "Produced %d seeds from this stream for future usage\n", length(seeds))
    }

    fitList <- list()
    for (jj in seq_len(nbrOfSegments)) {
      seg <- knownSegments[jj,]
      chromosomeJJ <- seg$chromosome
      xStart <- seg$start
      xEnd <- seg$end
      segTag <- sprintf("chr%s:(%s,%s)", chromosomeJJ, xStart, xEnd)
      verbose && enter(verbose, sprintf("Segment #%d ('%s') of %d", jj, segTag, nbrOfSegments), level=-10)

      ## Nothing to do?
      isSplitter <- (is.na(xStart) && is.na(xEnd))
      if (isSplitter) {
        fit <- splitter
        verbose && cat(verbose, "Nothing to segment. Inserting an explicit splitter.", level=-10)
        fitList[[segTag]] <- fit
        fit <- NULL
        verbose && exit(verbose)
        next
      }

      # Extract subset of data and parameters for this segment
      dataJJ <- subset(data, chrom == chromosomeJJ & xStart <= x & x <= xEnd)
      verbose && str(verbose, dataJJ, level=-50)
      chrom <- x <- index <- y <- w <- NULL
      fields <- attachLocally(dataJJ, fields=c("chrom", "x", "index", "y", "w"))
      dataJJ <- NULL # Not needed anymore
      nbrOfLoci <- length(y)

      # Empty segment?
      # [AD HOC. Should be done by segmentCBS(). /HB 2011-10-21]
      if(nbrOfLoci == 0) {
        fit <- splitter
        fit$output$chromosome <- chromosomeJJ
        fit$output$start <- xStart
        fit$output$end <- xEnd
        fit$output$nbrOfLoci <- nbrOfLoci
        fitList[[segTag]] <- fit
        fit <- NULL
        verbose && exit(verbose)
        next
      }

      seedJJ <- seeds[[jj]]

      fitList[[segTag]] <- future({
        fit <- segmentByCBS(y=y,
                  chromosome=chrom, x=x,
                  w=w,
                  index=index,
                  undo=undo,
                  avg=avg,
                  joinSegments=joinSegments,
                  knownSegments=seg,
                  ...,
                  seed=seedJJ,
                  verbose=less(verbose,1))

        # Sanity checks
        if (R_SANITY_CHECK) {
          .stop_if_not(nrow(fit$data) == nbrOfLoci)
          .stop_if_not(all.equal(fit$data$y, y))

          # Assert weights were used
          .stop_if_not(!hasWeights || !is.null(fit$data$w))
        } # if (R_SANITY_CHECK)

        segs <- as.data.frame(fit)
        if (nrow(segs) < 6) {
          verbose && print(verbose, segs, level=-10)
        } else {
          verbose && print(verbose, head(segs), level=-10)
          verbose && print(verbose, tail(segs), level=-10)
        }

        # Sanity check
        if (R_SANITY_CHECK) {
          ## NOTE: Need to use PSCBS::nbrOfSegments() for this to work
          ## when 'PSCBS' is not attached. The reason is that the future
          ## picks up the local integer variable 'nbrOfSegments' instead
          ## of the nbrOfSegments() function of 'PSCBS' as a global
          ## variable.
          ## https://github.com/HenrikBengtsson/PSCBS/issues/67
          .stop_if_not(TRUE && PSCBS::nbrOfSegments(fit, splitters=TRUE) > 0)
        } # if (R_SANITY_CHECK)

        fit
      }, seed = TRUE, label = sprintf("segmentByCBS-%s", segTag))

      rm(list=fields) # Not needed anymore

      verbose && exit(verbose)
    } # for (jj ...)


    verbose && enter(verbose, "Merging (independently) segmented known segments", level=-10)
    verbose && cat(verbose, "Number of segments: ", length(fitList), level=-10)
    fitList <- value(fitList)
    verbose && str(verbose, fitList, level=-50)
    ## former Reduce() w/ append(..., addSplit = FALSE)
    fit <- do.call(c, args = c(fitList, addSplit = FALSE))
    # Not needed anymore
    fitList <- NULL

    # Update parameters that otherwise may be incorrect
    fit$params$seed <- seed

    verbose && str(verbose, fit, level=-10)
    verbose && exit(verbose)

    segs <- getSegments(fit)
    if (nrow(segs) > 6) {
      verbose && print(verbose, head(segs), level=-10)
      verbose && print(verbose, tail(segs), level=-10)
    } else {
      verbose && print(verbose, segs, level=-10)
    }

    # Sanity checks
    if (R_SANITY_CHECK) {
      segs <- getSegments(fit)
      .stop_if_not(all(segs$start[-1] >= segs$end[-nrow(segs)], na.rm=TRUE))
      .stop_if_not(all(diff(segs$start) >= 0, na.rm=TRUE)) ## FIXME: > 0
      .stop_if_not(all(diff(segs$end) >= 0, na.rm=TRUE)) ## FIXME: > 0

  #    if (nrow(fit$data) != length(y)) {
  #      print(c(nrow(fit$data), nrow(data)))
  #    }
  #    .stop_if_not(nrow(fit$data) == nrow(data))
  #    .stop_if_not(all(fit$data$chromosome == data$chromosome))
  #    .stop_if_not(all(fit$data$x == data$x))
  #    .stop_if_not(all(fit$data$index == data$index))
  #    .stop_if_not(all.equal(fit$data$y, data$y))
    } # if (R_SANITY_CHECK)

    verbose && exit(verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Return results
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return(fit)
  } # if (nbrOfSegments > 1)

  nbrOfSegments <- nrow(knownSegments)

  # Sanity check
  if (R_SANITY_CHECK) {
    .stop_if_not(nbrOfSegments <= 1)
  } # if (R_SANITY_CHECK)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Specific segment?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (nbrOfSegments > 0) {
    knownSegments <- subset(knownSegments, chromosome == chromosome)
    nbrOfSegments <- nrow(knownSegments)
    # Sanity check
    if (R_SANITY_CHECK) {
      .stop_if_not(nbrOfSegments <= 1)
    } # if (R_SANITY_CHECK)
  }

  if (nbrOfSegments == 1) {
    seg <- knownSegments[1,]
    chromosomeJJ <- seg$chromosome
    xStart <- seg$start
    xEnd <- seg$end
    segTag <- sprintf("chr%s:(%s,%s)", chromosomeJJ, xStart, xEnd)
    verbose && printf(verbose, "Extracting segment '%s'", segTag, level=-50)

    # Extract subset of data and parameters for this segment
    data <- subset(data, chrom == chromosomeJJ & xStart <= x & x <= xEnd)
    verbose && str(verbose, data, level=-50)
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving the fit function", level=-50)
  # We need to attach the 'DNAcopy' package
  pkgName <- "DNAcopy"
  use(pkgName)
  pkg <- packageDescription(pkgName)
  pkgVer <- pkg$Version
  pkgDetails <- sprintf("%s v%s", pkgName, pkgVer)

  methodName <- "segment"
  verbose && cat(verbose, "Method: ", methodName, level=-50)
  verbose && cat(verbose, "Package: ", pkgDetails, level=-50)

  # Get the fit function for the segmentation method
#  fitFcn <- getExportedValue(pkgName, methodName)
  fitFcn <- getFromNamespace(methodName, pkgName)
  verbose && str(verbose, "Function: ", fitFcn, level=-50)
  formals <- formals(fitFcn)
  verbose && cat(verbose, "Formals:", level=-50)
  verbose && str(verbose, formals, level=-50)
  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up arguments to pass to segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setting up method arguments", level=-50)

  verbose && enter(verbose, "Setting up ", pkgName, " data structure", level=-50)

  sampleName <- "y" # This is going to be the name of the data field

  # Supress all warnings, in order to avoid warnings by DNAcopy::CNA()
  # on "array has repeated maploc positions".  Ideally we should filter
  # just those out. /HB 2013-10-22
  suppressWarnings({
    cnData <- CNA(
      genomdat  = data$y,
      chrom     = data$chrom,
      data.type = "logratio",
      maploc    = data$x,
      sampleid  = sampleName,
      presorted = TRUE
    )
  })
  verbose && str(verbose, cnData, level=-50)
  names(cnData)[3] <- sampleName
  verbose && str(verbose, cnData, level=-50)
  verbose && exit(verbose)

  # Sanity check
  if (R_SANITY_CHECK) {
    # (because all loci with unknown locations have already been dropped)
    .stop_if_not(nrow(cnData) == nrow(data))
  } # if (R_SANITY_CHECK)


  userArgs <- list(...)
  if (length(userArgs) > 0) {
    verbose && cat(verbose, "User arguments:", level=-50)
    verbose && str(verbose, userArgs, level=-50)
  }

  # Check if 'sbdry' can/should be precalculated.  This uses memoization
  # so that next time you segment with same 'nperm', 'alpha' and 'eta'
  # parameters, there will be much less startup overhead.
  if (length(userArgs) > 0 && !is.element("sbdry", names(userArgs))) {
    keys <- c("nperm", "alpha", "eta")
    keep <- is.element(keys, names(userArgs))
    if (any(keep)) {
      verbose && enter(verbose, "Precalculating argument 'sbdry' (with memoization)", level=-50)
      # Precalculate boundaries
      argsT <- formals[keys]
      keys <- keys[keep]
      argsT[keys] <- userArgs[keys]
      argsT$verbose <- less(verbose, 5)
      sbdry <- do.call(getbdry2, args=argsT)
      userArgs$sbdry <- sbdry
      verbose && exit(verbose)
    }
  }

  params <- list()

  if (hasWeights) {
    params$weights <- data$w
  }

  if (undo > 0) {
    params$undo.splits <- "sdundo"
    params$undo.SD <- undo
  }

  verbose && cat(verbose, "Segmentation parameters:", level=-50)
  verbose && str(verbose, params, level=-50)

  # Assign/overwrite by user arguments
  if (length(userArgs) > 0) {
    for (ff in names(userArgs)) {
      params[[ff]] <- userArgs[[ff]]
    }
  }

  verbose && cat(verbose, "Segmentation and user parameters:", level=-50)
  verbose && str(verbose, params, level=-50)

  # Cleaning out unknown parameters
  keep <- (names(params) %in% names(formals))
  params <- params[keep]

  args <- c(list(cnData), params, verbose=as.logical(verbose))
  verbose && cat(verbose, "Final arguments:", level=-50)
  verbose && str(verbose, args, level=-50)

  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calling segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, sprintf("Calling %s() of %s", methodName, pkgName), level=-50)

  # There are a few cases where we can/need to do a dummy segmentation
  # based on a single data points:
  # (a) WORKAROUND for the case when there are no data points.
  # (b) SPEEDUP: When undo=+Inf we don't really have to segment.
  nbrOfNonMissingLoci <- sum(!is.na(cnData$y))
  if (nbrOfNonMissingLoci == 0) {
    args[[1]] <- CNA(genomdat=0, chrom=0, maploc=0)
    if (hasWeights) args$weights <- 1.0
  } else if (undo == +Inf) {
    args[[1]] <- CNA(genomdat=0, chrom=0, maploc=0)
    if (hasWeights) args$weights <- 1.0
    verbose && cat(verbose, "Skipping identification of new change points (undo=+Inf)", level=-50)
  }


  ## Set random seed?
  if (!is.null(seed)) {
    randomSeed("set", seed=seed, kind="L'Ecuyer-CMRG")
    on.exit(randomSeed("reset"), add=TRUE)
    verbose && printf(verbose, "Random seed temporarily set (seed=c(%s), kind=\"L'Ecuyer-CMRG\")\n", paste(seed, collapse=", "))
  }

  # In case the method writes to stdout, we capture it
  # Note: DNAcopy::segment() *does* this.
  stdout <- capture.output({
    # Does not work, because some internal function of the fit function
    # may only be accessible from within the namespace
    # How to do this for DNAcopy::segment()? /HB
##    fit <- do.call(fitFcn, args)
    # This works, but requires that one loads the package and that the
    # function is not masked in the search() path.
    t <- system.time({
      fit <- do.call(methodName, args)
    }, gcFirst = FALSE)
    # Drop the 'call' (because it will be huge due to the do.call() call)
    fit$call <- NULL
  })
  attr(fit, "processingTime") <- t
  attr(fit, "pkgDetails") <- pkgDetails
  attr(fit, "randomSeed") <- seed

  # WORKAROUND for the case when there are no data points.
  if (nbrOfNonMissingLoci == 0) {
    # Drop dummy data point...
    fit$data <- cnData ## fit$data[-1,,drop=FALSE]
    # ...dummy region found
    output <- fit$output
    segRows <- fit$segRows

    # Sanity check
    if (R_SANITY_CHECK) {
      .stop_if_not(nrow(output) == 1)
    } # if (R_SANITY_CHECK)

    # Was a region specified?
    if (nbrOfSegments == 1) {
      seg <- knownSegments[1,]
      output$ID <- sampleName
      output$chrom <- seg$chromosome
      if (is.finite(seg$start)) {
        output$loc.start <- seg$start
      }
      if (is.finite(seg$end)) {
        output$loc.end <- seg$end
      }
      output$num.mark <- 0L
      output$seg.mean <- NA_real_
      segRows[1,] <- NA_integer_
    } else {
      output <- output[-1,,drop=FALSE]
      segRows <- segRows[-1,,drop=FALSE]
    }
    fit$output <- output
    fit$segRows <- segRows
  } else if (undo == +Inf) {
    # Drop dummy data point...
    fit$data <- cnData ## fit$data[-1,,drop=FALSE]
    # ...dummy region found
    output <- fit$output
    segRows <- fit$segRows

    # Sanity check
    if (R_SANITY_CHECK) {
      .stop_if_not(nrow(output) == 1)
    } # if (R_SANITY_CHECK)

    # Was a region specified?
    if (nbrOfSegments == 1) {
      seg <- knownSegments[1,]
      output$ID <- sampleName
      output$chrom <- seg$chromosome
      if (is.finite(seg$start)) {
        output$loc.start <- seg$start
      } else {
        output$loc.start <- min(cnData$maploc, na.rm=TRUE)
      }
      if (is.finite(seg$end)) {
        output$loc.end <- seg$end
      } else {
        output$loc.end <- max(cnData$maploc, na.rm=TRUE)
      }
    }
    output$num.mark <- nrow(fit$data)
    output$seg.mean <- mean(fit$data$y, na.rm=TRUE)
    segRows$endRow <- nrow(fit$data)

    fit$output <- output
    fit$segRows <- segRows
  } # if (undo == +Inf)

  verbose && cat(verbose, "Captured output that was sent to stdout:", level=-50)
  stdout <- paste(stdout, collapse="\n")
  verbose && cat(verbose, stdout, level=-50)

  verbose && cat(verbose, "Fitting time (in seconds):", level=-50)
  verbose && print(verbose, t, level=-50)

  verbose && cat(verbose, "Fitting time per 1000 loci (in seconds):", level=-50)
  verbose && print(verbose, 1000*t/nbrOfLoci, level=-50)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Restructure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Restructuring results", level=-50)

  # Coerce
  fit$output$num.mark <- as.integer(fit$output$num.mark)

  # Coerce 'chrom' to a plain integer
  fit$data$chrom <- unclass(fit$data$chrom)

  # Store genomewide index
  fit$data$index <- data$index

  # Store weights
  fit$data$w <- data$w

  # Not needed anymore
  data <- NULL

  verbose && exit(verbose)


  # Store also interesting parameters to DNAcopy::segment()
  keys <- setdiff(names(formals), c("x", "weights", "sbdry", "verbose"))
  keys <- c(keys, "undo", "seed")
  keep <- is.element(names(params), keys)
  keep <- names(params)[keep]
  params <- params[keep]
  params$undo <- undo
  params$joinSegments <- joinSegments
  params$knownSegments <- knownSegments
  params$seed <- seed
  fit$params <- params

#  class(fit) <- c("CBS", class(fit))
  class(fit) <- c("CBS", "AbstractCBS")

  # Sanity checks
  if (R_SANITY_CHECK) {
    segRows <- fit$segRows
    .stop_if_not(all(segRows[,1] <= segRows[,2], na.rm=TRUE))
    .stop_if_not(all(segRows[-nrow(segRows),2] < segRows[-1,1], na.rm=TRUE))
  } # if (R_SANITY_CHECK)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Renaming column names
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit)
  names <- colnames(data)
  names <- gsub("chrom", "chromosome", names, fixed=TRUE)
  names <- gsub("maploc", "x", names, fixed=TRUE)
  colnames(data) <- names

  # Drop 'CNA' class and DNAcopy attributes
  class(data) <- c("data.frame")
  attr(data, "data.type") <- NULL

  fit$data <- data

  segs <- fit$output

  names <- colnames(segs)
  names <- gsub("ID", "sampleName", names, fixed=TRUE)
  names <- gsub("seg.mean", "mean", names, fixed=TRUE)
  names <- gsub("chrom", "chromosome", names, fixed=TRUE)
  names <- gsub("num.mark", "nbrOfLoci", names, fixed=TRUE)
  names <- gsub("loc.", "", names, fixed=TRUE) # loc.start, loc.end
  colnames(segs) <- names
  fit$output <- segs


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Join segments?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (joinSegments) {
    if (nbrOfSegments == 1) {
      starts <- knownSegments$start
      ends <- knownSegments$end
      if (is.infinite(starts)) starts <- segs$start
      if (is.infinite(ends)) ends <- segs$end
      range <- range(c(starts, ends), na.rm=TRUE)
    } else {
      range <- NULL
    }

    fit <- joinSegments(fit, range=range, verbose=less(verbose, 10))

    # Sanity checks
    if (R_SANITY_CHECK) {
      segRows <- fit$segRows
      .stop_if_not(all(segRows[,1] <= segRows[,2], na.rm=TRUE))
      .stop_if_not(all(segRows[-nrow(segRows),2] < segRows[-1,1], na.rm=TRUE))
    } # if (R_SANITY_CHECK)
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (avg != "mean") {
    verbose && enter(verbose, "Updating mean level using different estimator")
    verbose && cat(verbose, "Estimator: ", avg)
    fit <- updateMeans(fit, avg=avg, verbose=less(verbose, 20))
    verbose && exit(verbose)
  }


  verbose && cat(verbose, "Results object:", level=-10)
  verbose && str(verbose, fit, level=-10)

  verbose && exit(verbose)


  verbose && exit(verbose)

  fit
}) # segmentByCBS()


setMethodS3("segmentByCBS", "data.frame", function(y, ...) {
  # To please R CMD check
  data <- y

  y <- data$y
  if (is.null(y)) {
    y <- data$cn
    if (is.null(y)) {
      y <- data$CT
    }
  }

  segmentByCBS(y=y, chromosome=data$chromosome, x=data$x, index=data$index, w=data$w, ...)
})


setMethodS3("segmentByCBS", "CBS", function(...) {
  resegment(...)
}) # segmentByCBS()
