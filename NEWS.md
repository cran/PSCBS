# Version 0.68.0 [2025-04-18]

## Documentation

 * Fix minor help-page issues.

## Bug Fixes

 * `segmentByCBS()` had an internal bug which triggered an error on
  `could not find function "nbrOfSegments"`, which began with
  **future** 1.40.0.
 

# Version 0.67.0 [2024-02-16]

## Documentation

 * Fix minor help-page issues.

## Miscellaneous

 * Package no longer suggests the **Hmisc** package, because it's
   actually not used. See NEWS entry for v0.66.0 for details.

## Deprecated and Defunct

 * Removed defunct `append()` for `AbtractCBS`.

 * Removed defunct `load()` and `save()` for `AbtractCBS`.
 
 * After removing the above functions, `library(PSCBS)` no longer
   reports on functions being masked.
   

# Version 0.66.0 [2021-10-22]

## Significant Changes

 * Package no longer require the **Hmisc** package, which was used for
   its `wtd.quantile()` function. Instead, we have adopted its GPL
   (>=2) code.  The reason for doing so is that **Hmisc** no longer
   installs out of the box on platforms, e.g. macOS M1, but also that
   **Hmisc** has a large number of package dependencies, which adds an
   unnecessary installation weight to the **PSCBS** package for this
   single function.
   
## Code Refactoring

 * Update package tests to test with 'multisession' futures.
 
## Bug Fixes

 * `segmentByCBS()` and `segmentByPairedPSCBS()` would produce
   warnings on possibly unreliable random numbers due to parallel
   processing in the case they segmented multiple chromosomes or
   segments.

 * Fixed a few partial `getSegments()` argument name matching
   (`splitter` -> `splitters`).

## Deprecated and Defunct

 * `append(x, y)` for `CBS` and `PSCBS` objects, deprecated since
   v0.64.0 is defunct; use `c(x, y)` instead.

 * `save()` for `CBS` and `PSCBS` objects and corresponding
   `CBS$load()` and `PSCBS$load()` methods are defunct.  We recommend
   to use base-R functions `saveRDS()` and `readRDS()` instead.  If
   backward compatibility with the defunct `save()` and `load()`
   methods, the `saveObject()` and `loadObject()` functions from the
   **R.utils** package can be used.


# Version 0.65.0 [2019-05-04]

## New Features

 * The PDF report produced by `report()` for Paired PSCBS results will
   now provide an estimate of normal contamination together with a
   disclaimer.

## Bug Fixes

 * A sanity check in `segmentByPairedPSCBS()` could produce an error
   on `'length(x) = 5 > 1' in coercion to 'logical(1)'` when running
   with `_R_CHECK_LENGTH_1_LOGIC2_=true`. This bug did _not_ affect
   the results of PSCBS.

 * Report templates used by `report()` would produce error `Error in
   unit(x, default.units) : 'x' and 'units' must have length > 0` if
   there are too few loci to plot.  This would typically happen when
   reporting on human chromosome 25 results.

## Deprecated and Defunct

 * Argument `preserveScale` of `segmentByPairedPSCBS()` is now
   defunct.
 

# Version 0.64.0 [2018-08-12]

## New Features

 * Added `c()` for `CBS` and `PSCBS` objects.

## Performance

 * `segmentByCBS()` no longer performs garbage collection, which
   happened indirectly via a `system.time()` call that does GC by
   default.

## Bug Fixes

 * `plotTrack()` for `CBS` objects would produce error `Argument
   'Clim' is not a vector: NULL` when the signal type was unknown.
   Now it will assume (non-logged) copy-number ratios so it can choose
   a default `Clim` range.

## Deprecated and Defunct

 * Removed `bootstrapDHByRegion()`, which was defunct since 0.44.0
   (Feb 2015).

 * `append(x, y)` for `CBS` and `PSCBS` objects is deprecated; use
   `c(x, y)` instead.

 * Use of argument `preserveScale` for `segmentByPairedPSCBS()` is now
   deprecated and ignored.  it's value is now fixed to FALSE, which
   has been the default since **PSCBS** 0.50.0 (Oct 2015).
   

# Version 0.63.0 [2017-06-27]

## Significant Changes

 * Package now depends on R (>= 3.2.0) and Bioconductor (>= 3.1)
   packages.

## New Features

 * `estimateDeltaCN()` for `PairedPSCBS` gained argument `flavor` and
   new estimator `flavor = "delta(mode)"`.

 * Added `isLocallyPhased()` for `PSCBS`.

## Code Refactoring

 * Now package imports **aroma.light** instead of only suggesting it.

 * Package tests no longer test against the deprecated 'lazy' strategy
   of **future**.


# Version 0.62.0 [2016-11-10]

## New Features

 * Added `normalizeTotalCNs()` for `PSCBS` objects.

 * REPORTS: Updated report template for `PairedPSCBS` object such that
   reports are also generated from DH-only data, i.e. when there are
   no BAF signals (which may happen with DNA-Seq data).

 * Added `splitters = TRUE` as the default for `setSegments()`.

## Code Refactoring

 * CLEANUP: Minor internal cleanup and restructuring.



# Version 0.61.0 [2016-02-03]

## Significant Changes

 * Package now requires R (>= 3.1.2) released October 2014, because of
   its dependency on the **listenv** package.

## New Features

 * `segmentByPairedPSCBS()` gained argument `rho` for paired PSCBS
   segmentation on total CNs (TCNs) and decrease-of-heterozygosity
   signals (DHs) as an alternative to for instance TCN and allele B
   fractions (BAFs).

## Bug Fixes

 * Segmentation using futures where not fully reproducible when known
   segments where specified.  Fixed by changing how the random number
   stream is set and used.


# Version 0.60.0 [2015-11-17]

## New Features

 * PARALLEL: Add support for parallel processing via futures by
   utilizing the **future** package.  Parallel segmentation over
   multiple chromosomes (or known segments) done by `segmentByCBS()`
   and `segmentByPairedPSCBS()` is enabled by
   `future::plan("multicore")`.

 * REPRODUCIBILITY: Whenever argument `seed` is given, the
   L'Ecuyer-CMRG stream is now used to generate random numbers.  For
   backward compatibility with other types of random number
   generators, don't specify argument `seed` but instead use
   `set.seed()` to set the seed before calling the method.

## Code Refactoring

 * Bump package dependencies.


# Version 0.50.0 [2015-10-14]

## Significant Changes

 * Now argument `preserveScale` for `segmentByPairedPSCBS()` defaults
   to FALSE.  In the past the default has effectively been TRUE, but
   has given a warning since v0.43.0 (June 2014) that it eventually
   will be changed to FALSE.  Now it will give a warning that it has
   changed, unless the option is explicitly specified.  This new
   warning will eventually be removed in a future version.


# Version 0.45.1 [2015-09-16]

## New Features

 * More informative error messages when `append()` for `CBS` or
   `PSCBS` fail.

## Bug Fixes

 * `segmentByCBS(, ..., w, knownSegments)` would give internal
   assertion errors if one of the priorly known segments would have
   zero data points.  Thanks to Kirill Tsukanov (Moscow) and Eric
   Talevich (UCSF) for reporting on this.


# Version 0.45.0 [2015-09-11]

## Significant Changes

 * Package now requires R (>= 3.1.1) released July 2014. This allows
   us to use Bioconductor (>= 3.0) (October 2014).
  
## New Features

 * `segmentByCBS()` gained argument `avg`.

 * Add `writeWIG()` for `CBS` objects.

 * `pruneByHClust()` no longer gives a message about method `"ward"`
   is now named `"ward.D"`.

 * Added `skip = TRUE` to `report()`.

## Bug Fixes

 * `plotTracks()` for `CBS` ignored arguments `cex`, `col` and
   `meanCol` if two or more chromosomes were plotted.

 * `joinSegments()`, `resetSegments()`, and `pruneBySdUndo()` gave
   errors for multi- chromosome (>= 2) segmentation results.

 * `segmentByCBS()` would ignore argument `w` (weights) if more than one
   chromosome was fitted.

 * `tileChromosomes()` for `CBS` returned incorrect locus data.

 * `gapsToSegments(gaps)` gave an error if `nrow(gaps) == 0`, or
   it contained no `chromosome` column.

 * `findLargeGaps()` could return NULL. Now it always returns a
   data.frame.

 * The `report()` RSP-embedded TeX templates for `CBS` and
   `PairedPSCBS` data did not escape sample and data set names to
   LaTeX in all places needed.

## Code Refactoring

 * ROBUSTNESS: Package test coverage is 62%.

 * ROBUSTNESS: Explicitly importing core R functions.


# Version 0.44.0 [2015-02-22]

 * Package now requires R (>= 3.0.3) and Bioconductor (>= 2.13), which
   are from March 2014 and are in fact old. it's recommended to use a
   more recent R version.

## Deprecated and Defunct

 * CLEANUP: `bootstrapDHByRegion()` is defunct (was deprecated since
   2013).

## Code Refactoring

 * ROBUSTNESS: Package test coverage is 58%.

 * ROBUSTNESS: Forgot to declare some S3 methods in NAMESPACE.

 * SPEEDUP: Now using more functions of **matrixStats**.


# Version 0.43.0 [2014-06-08]

## Significant Changes

 * Now `segmentByPairedPSCBS()` gives a warning about future change of
   the default value of argument `preserveScale` (from current TRUE to
   FALSE).  The warning only appears if the argument is not specified
   explicitly.

 * Package now requires R (>= 3.0.0) and Bioconductor (>= 2.13), which
   were released April 2013 and are in fact old and it's recommended
   to use a more recent version of R.

## Code Refactoring

 * Now using `use()` of **R.utils** where possible.

 * Bumped package dependencies.


# Version 0.42.2 [2014-05-24]

## Code Refactoring

 * Bumped package dependencies.


# Version 0.42.1 [2014-05-05]

## Bug Fixes

 * `pruneByHClust()` for `PairedPSCBS` would give an error on `unable
   to find an inherited method for function 'anyMissing' for signature
   '"PairedPSCNSegments"'`, unless the object contained bootstrap
   statistics.  This is no longer needed.  Thanks to Junsong Zhao, Los
   Angeles, CA for reporting on this.


# Version 0.42.0 [2014-04-25]

## Code Refactoring

 * Minor speedup (a few percents) by now byte compiling the package by
   default.

 * CLEANUP: Dropped unnecessary usage of `::`.

 * Bumped package dependencies.


# Version 0.41.4 [2014-03-30]

## New Features

 * GENERALIZATION: Now `callROH()` also works if paired PSCBS was done
   with only `muN` available (and not `betaN`).  In that case, it
   assumes that all genotype confidence scores are equal.

 * Updated the ordering and the defaults of `testROH()` arguments to
   make it clear that `betaN` is optional and only used if `csN` is
   not given.

 * As an alternative to argument `CT`, `segmentByPairedPSCBS()` now
   accepts arguments `thetaT` and `thetaN`, in case `CT` is calculated
   as `CT = 2 * thetaT / thetaN`.


# Version 0.41.3 [2014-03-29]

## New Features

 * Methods no longer generates warnings on `in max(c(NA_integer_,
   NA_integer_), na.rm = TRUE) : no non-missing arguments to max;
   returning -Inf`.

## Bug Fixes

 * In rare cases, `callROH()` could throw `Error in if (is.na(delta))
   { : argument is of length zero`.


# Version 0.41.2 [2014-03-28]

## New Features

 * Added argument `preserveScale` to `segmentByPairedPSCBS()`, which
   is passed as is to `normalizeTumorBoost()` with the default being
   TRUE corresponding to the previous default behavior.


# Version 0.41.1 [2014-03-28]

## New Features

 * Added `unTumorBoost()` to recalculating the segment means and other
   statistics for a given `PairedPSCBS` profile based on
   non-TumorBoosted tumor BAFs (rather than TumorBoost:ed tumor BAFs).


# Version 0.41.0 [2014-03-26]

## New Features

 * Now `estimateKappaByC1Density()` give more informative error
   messages if it failed to identify modes for estimating the
   parameter.

 * Added argument `from` to `estimateKappaByC1Density()`.

## Bug Fixes

 * `updateMeansC1C2()` for `PairedPSCBS` did not handle missing values
   (= "splitters") in the `c1c2Swap` field.

 * `updateMeans()` for `PairedPSCBS` and `NonPairedPSCBS` returned the
   incorrect DH segment levels for region in AB if `adjustFor = "ab"`
   and likewise for segments in LOH if `adjustFor = "loh"`.  This bug
   does *not* affect any of **PSCBS** methods themselves, because none
   of them utilize those `adjustFor` options.

## Code Refactoring

 * Bumped package dependencies.


# Version 0.40.4 [2014-02-04]

## Bug Fixes

 * `all.equal()` for `CBS` would pass the first/dispatch argument to
   `NextMethod()` as `target = target` and not as `object = target`,
   which would result in it being passed it twice both named and
   non-named where the latter would become argument `tolerance =
   target` in an internal call to `all.equal()` for numerics.  In
   recent R-devel version this would generate `Error in
   all.equal.numeric(target[[i]], current[[i]], check.attributes =
   check.attributes, : 'tolerance' should be numeric Calls: stopifnot
   ...  all.equal.default -> all.equal.list -> all.equal ->
   all.equal.numeric`.


# Version 0.40.3 [2014-01-29]

## New Features

 * ROBUSTNESS: Now `segmentByPairedPSCBS()` asserts that argument
   `muN` is not all NAs.  Similarily, if `muN` is called from `betaN`
   the same assertion is done after calling.


# Version 0.40.2 [2013-12-17]

## New Features

 * Now `estimateDeltaCN()` for `CBS` have the option to estimate the
   size of a copy-number unit based on the change-point magnitutes, in
   addition to the estimator based on the density of segment means.

## Bug Fixes

 * `getChangePoints()` for `CBS` returned empty results.


# Version 0.40.1 [2013-12-09]

## Documentation

 * The CBS vignette referred to C1 and C2 in one of the code examples.

## Code Refactoring

 * Bumped package dependencies.


# Version 0.40.0 [2013-12-07]

## Code Refactoring

 * CLEANUP: No longer a need for an ad-hoc NAMESPACE import.


# Version 0.39.8 [2013-12-04]

## Documentation

 * Now the vignette sections on dropping outliers before segmentation
   explains why outliers are dropped whereas in the original CBS
   publication they were shrunk ("smoothed"). Also, added help for
   `dropSegmentationOutliers()`.


# Version 0.39.7 [2013-11-27]

## New Features

 * Added `callGLAO()` for `CBS`.

 * Added `encodeCalls()` for `data.frame` object returned by
   `getLocusData(..., addCalls = TRUE)`.

## Code Refactoring

 * Bumped package dependencies.


# Version 0.39.6 [2013-11-23]

## New Features

 * Added `clearCalls()` for `AbstractCBS`.

 * Added `extractSegmentDataByLocus()` for `PairedPSCBS`.

## Bug Fixes

 * `estimateDeltaCN()` for `CBS` assumed **aroma.light** was attached.


# Version 0.39.5 [2013-11-15]

## New Features

 * Added `estimateDeltaCN()` for `CBS`.  Added package system test.


# Version 0.39.4 [2013-11-14]

## Bug Fixes

 * `callGainsAndLosses()` for `CBS` would not estimate the median
   median CN level correctly if there were "empty" segments
   (e.g. gaps).  This was/is due to a bug in `segments.summary()` of
   the **DNAcopy** package.  Instead, we are now calculating the
   segment median levels ourselves.  Added a system package test for
   `callGainsAndLosses()`.


# Version 0.39.3 [2013-11-05]

## New Features

 * Added basic implementations of `setLocusData()` and `setSegments()`
   for `AbstractCBS`.


# Version 0.39.2 [2013-10-28]

## New Features

 * Now `plotTracksManyChromosomes()` for `PairedPSCBS` also supports
   tracks `"c1,c2"`, `"c1"`, and `"c2"`.


# Version 0.39.1 [2013-10-25]

## New Features

 * Now `plotTracksManyChromosomes()` uses the locus data field `rho`
   when plotting DH locus-level data.  It only recalculates it from
   the tumor BAFs if the DH signals are not available - if so a
   warning is generated.

## Bug Fixes

 * The `rho` signals returned by `getLocusData(..., fields = "full")`
   for `PairedPSCBS` would have values also for homozygote SNPs.


# Version 0.39.0 [2013-10-23]

## New Features

 * Now all warnings generated by `DNAcopy::CNA()` are suppressed,
   including the common one on `array has repeated maploc positions`.

 * Added `getBootstrapLocusSets()` for `PairedPSCBS`.  Added a package
   system test for it.

 * Added argument `subset` to `applyByRegion()` for `PairedPSCBS`.

 * Added `clearBootstrapSummaries()` for `PairedPSCBS`.

 * SPEEDUP: Added argument `cache` to
   `bootstrapSegmentsAndChangepoints()`, which caches the results to
   file if `cache = TRUE`.


# Version 0.38.6 [2013-10-20]

## Bug Fixes

 * `plotTracks()` for `PairedPSCBS` would use argument `Clim` for
   `Blim` as well, regardless of what argument `Blim` is.  This bug
   was introduced in v0.38.3.

## Code Refactoring

 * Internal restructuring on how bootstrapping of segment means is
   done.


# Version 0.38.5 [2013-10-18]

## Bug Fixes

 * The CBS and Paired PSCBS report templates assumed that the
   **R.utils** package is attached.


# Version 0.38.4 [2013-10-15]

## Code Refactoring

 * CLEANUP: Removed a few unnecessary NAMESPACE imports.

 * Bumped package dependencies.


# Version 0.38.3 [2013-10-14]

## New Features

 * Now `plotTracks()` for `CBS` and `PSCBS` gives a more informative
   error if `Clim` or `Blim` is invalid.  If using `"auto"` (only for
   `CBS`) and the limits could not be inferred due to an unknown or
   unset signal type, an informative error message reports on this as
   well.

## Code Refactoring

 * Now the package vignettes are in `vignettes/`, and not in
   `inst/doc/`, which will not be supported by R (>= 3.1.0).

 * ROBUSTNESS: The overriding of `append()` to become a generic
   function does now call `base::append()` in the default, instead of
   copy the latter.  All this will eventually be removed, when proper
   support for `c`, `[`, `[[` etc. has been added everywhere.

 * CLEANUP: Now explicitly importing only what is needed in NAMESPACE.


# Version 0.38.2 [2013-10-13]

## Bug Fixes

 * While attaching the package, it could cause a cyclic loading of
   namespaces.


# Version 0.38.1 [2013-10-08]

## New Features

 * Now `getSmoothLocusData()` for `CBS` also returns column `count`
   which specifies the number of (finite) loci averaged over in each
   bin.

## Documentation

 * Vignette 'Total copy-number segmentation using CBS' would display
   the same plot as vignette 'Parent-specific copy-number segmentation
   using Paired PSCBS'.

 * Renamed vignette 'Paired PSCBS' to 'Parent-specific copy-number
   segmentation using Paired PSCBS'.

## Bug Fixes

 * `tileChromosomes()` for `CBS` did not set the `tiledChromosomes`
   attribute due to a typo.  This caused `plotTracks()` for `CBS` to
   horizontally misplace the plotted segment levels. Added a system
   tests for this for `CBS` and `PairedPSCBS` objects.  Thanks to
   Ilari Scheinin at VUMC for reporting on this.

## Code Refactoring

 * Bumped package dependencies.


# Version 0.38.0 [2013-09-27]

## Code Refactoring

 * SPEEDUP: `R CMD check` is now significantly faster due to copying of
   pre-generated calculations ("memoization"). For instance, the the same
   segmentation tests are roughly 40% faster compared to version 0.37.2.

 * Now **PSCBS** imports **R.cache** (used to only suggest it).


# Version 0.37.2 [2013-09-27]

## Code Refactoring

 * SPEEDUP: Now utilizing **matrixStats** functions in more places.

 * ROBUSTNESS: Further improved how **aroma.light** is handled for
   backward compatibility.

 * Bumped package dependencies.


# Version 0.37.1 [2013-09-26]

## Code Refactoring

 * CLEANUP: Now package avoids attaching suggested packages such as
   **R.cache**, **aroma.light**, and **Hmisc** by only importing the set
   of functions needed via `::`.  This way those packages are only
   loaded.  Packages that still need to be attached are done so
   "quietly".

 * CLEANUP: Minor adjustments to some of the internal workarounds for
   older versions of **matrixStats** and **aroma.light**.

## Bug Fixes

 * Forgot to import several functions from **matrixStats**. These went
   undetected because **aroma.light** (<= 1.31.5) attaches the
   **matrixStats**.

 * `segmentByPairedPSCBS()` assumed **aroma.light** was attached.

 * One of the system tests assumed **R.utils** was attached.


# Version 0.37.0 [2013-09-21]

## Bug Fixes

 * WORKAROUND: For now, package attaches the **utils** package. This
   is needed due to what appears to be a bug in how **R.oo** finalizes
   Object:s assuming **utils** is attached, which may not be the case
   (unless **R.oo** itself is attached).

 * `callGNL()` for `PairedPSCBS` used non-defined `verbose` object.

## Code Refactoring

 * CLEANUP: Package no longer attaches **R.utils**, only imports it.

 * ROBUSTNESS: Now package imports only what is needed from
   **DNAcopy**.


# Version 0.36.2 [2013-09-18]

## Documentation

 * Added vignette 'Total copy-number segmentation using CBS'.

 * WORKAROUND: For R (< 3.0.0), `hclustCNs()` for `AbstractCBS` would
  generate `Error in rowAlls(ok) : could not find function
  "loadMethod"`. This seems to be a bug in R (< 3.0.0), which we can
  avoid by attaching the **methods** package in `hclustCNs()`.

## Code Refactoring

 * ROBUSTNESS: Now package imports **matrixStats** (previously suggested).

 * ROBUSTNESS: Now package declares S3 methods in the NAMESPACE.

 * ROBUSTNESS: Package vignettes no longer assumes that the **R.rsp**
   package is attached.

 * ROBUSTNESS: Forgot to import `R.methodsS3::appendVarArgs()`.

 * Bumped package dependencies.


# Version 0.36.1 [2013-09-10]

## Code Refactoring

 * CLEANUP: Package no longer utilizes `:::`.


# Version 0.36.0 [2013-08-15]

## New Features

 * Made `extractMinorMajorCNs()` for `PairedPSCBS` acknowledge
   additional fields related to (C1,C2).


# Version 0.35.6 [2013-08-01]

## Code Refactoring

 * Updated the vignettes to utilize the new **R.rsp** features.


# Version 0.35.5 [2013-07-19]

## New Features

 * ROBUSTNESS: Added a sanity check on the estimates of (tauA, tauB)
   when they are estimated from data in `segmentByNonPairedPSCBS()`.


# Version 0.35.4 [2013-07-11]

## Code Refactoring

 * Updated the `Makefile` for the vignettes and added `.Rinstignore`
   such that auxiliary (bib and bst) LaTeX files are not installed but
   part of the build so they are available to `R CMD check`, which is
   recently needed by R devel.

 * Bumped package dependencies.


# Version 0.35.3 [2013-05-25]

## Code Refactoring

 * Minor speedup by replacing all `rm(x)` with `x <- NULL`,
   cf. R-devel thread 'Assigning NULL to large variables is much
   faster than rm() - any reason why I should still use rm()?' on May
   25, 2013.


# Version 0.35.2 [2013-05-20]

## Documentation

 * CRAN POLICY: Now all Rd `\usage{}` lines are at most 90 characters
   long.


# Version 0.35.1 [2013-05-07]

## New Features

 * Now `estimateDeltaCN()` for `PairedPSCBS` adjust for the ploidy, if
   set.

 * Added `ploidy()` and `ploidy()<-` for `AbstractCBS`.

 * Now `tileChromosomes()` no longer gives warnings on `max(i): no
   non-missing arguments to max; returning -Inf`.


# Version 0.35.0 [2013-04-23]

## New Features

 * SPEEDUP: Now `bootstrapTCNandDHByRegion()` for `PairedPSCBS` always
   estimates the default quantiles in addition to any requested ones.

 * SPEEDUP: Made `bootstrapTCNandDHByRegion()` much faster by adding
   `use.names = FALSE` to two internal `unlist()` statements.

## Bug Fixes

 * `updateMeans()` for `PairedPSCBS` and `NonPairedPSCBS` could
   include a signal from a neighboring segment when averaging, iff
   that signal was located at the exact locus of the change
   point. Thanks Ingrid Lonnstedt (WEHI) for reporting on this.


# Version 0.34.9 [2013-04-22]

## Bug Fixes

 * `updateMeans()` would not always preserve the originally specified
   segment-mean level estimator, if different from a (sample) mean
   estimator, e.g. `avgDH = "median"`.  This could result in for
   instance `callAB()` failing on internal sanity checks.

 * Segment levels drawn by `plotTracks()` would have incorrect genomic
   locations for chromosome 2 and beyond.  This bug was introduced in
   v0.34.7.

## Code Refactoring

 * Utilizing new `startupMessage()` of **R.oo**.


# Version 0.34.8 [2013-04-20]

## Deprecated and Defunct

 * Removed previously deprecated methods for `AbstractCBS`.


# Version 0.34.7 [2013-04-18]

## New Features

 * Added more arguments to `plotTracks()`.

 * Now `drawLevels()` and `drawConfidenceBands()` for `CBS` and
   `PairedPSCBS` also works for multiple chromosomes.

## Bug Fixes

 * One of the system tests for `segmentByPairedPSCBS()` failed in the
   case when the data was downsampled (in order to meet the CRAN
   requirements).  The workaround is to use a fix `deltaAB` parameter
   in that case.

 * Internal `calcStatsForCopyNeutralABs()` would give an error if
   there was exactly two AB segments.


# Version 0.34.6 [2013-04-11]

## Bug Fixes

 * `plotTracks(fit, callLoci = TRUE)` would color loci incorrectly if
   more than one chromosome are plotted.


# Version 0.34.5 [2013-04-09]

## New Features

 * Now `callROH()` gives an informative error if called on a
   `NonPairedPSCBS` object.


# Version 0.34.4 [2013-04-05]

## New Features

 * Added more end-user control to `plotTracks()`.


# Version 0.34.3 [2013-04-04]

## Code Refactoring

 * Now package builds with both **R.rsp** (< 0.9.1) and **R.rsp** (>=
   0.9.1).


# Version 0.34.2 [2013-03-28]

## New Features

 * Now `callGainNeutralLoss()`, utilizes `callCopyNeutral()` by
   default.


# Version 0.34.1 [2013-03-21]

## New Features

 * Updated the report generator and its RSP templates.

## Documentation

 * Clarified in the **PSCBS** vignette that the NTCN caller is under
   development, experimental.

## Code Refactoring

 * SPEEDUP: Made `dropChangePoints()` faster by only updating the
   segment statistics/means at the very end.


# Version 0.34.0 [2013-03-19]

## New Features

 * CALLING: Defined a formal hypothesis test for how segments are
   called copy- neutral in TCN (NTCN), with the null hypothesis being
   that a segment is NTCN.  In order for a segment to not be NTCN, its
   confidence interval has to be completely outside the null region.
   This changed how `callCopyNeutralByTCNofAB()` for `PairedPSCBS`
   calls segments; it is now a bit more conservative in rejecting
   NTCN.


# Version 0.33.4 [2013-03-19]

## New Features

 * ROBUSTNESS: Now `calcStatsForCopyNeutralABs()` for `PairedPSCBS`
   does a better job in identifying the TCN mode of the AB segments.

 * Added argument `flavor` to `findNeutralCopyNumberState()`
   specifying how to identify the main mode of the AB segments.

 * VISUALIZATION: Now `plotTracks()` for `PairedPSCBS` displays
   thresholds for calling AB, LOH and and NTCN.


# Version 0.33.3 [2013-03-12]

## Documentation

 * Documented `tauA` and `tauB` in the help for
   `segmentByNonPairedPSCBS()`.


# Version 0.33.2 [2013-03-09]

## New Features

 * Added `getLocusData()` for `PairedPSCBS` and `NonPairedPSCBS`.

## Documentation

 * Updated the vignettes and the report templates to utilize the new
   **ggplot2** themes - **ggplot2** no longer gives a warning on using
   deprecated functions.

 * Now `report()` for `AbstractCBS` also includes files listed in the
   optional file `.install_extras` of the source RSP template
   directory.  The same filename is used by `R CMD build/check` for
   including additional source files needed to build the vignettes.

## Code Refactoring

 * Added an `Authors@R` field to the DESCRIPTION.


# Version 0.33.1 [2013-03-07]

## Documentation

 * Preparing package vignettes for the upcoming R 3.0.0 support for
   non-Sweave vignettes.

## Software Quality

 * Relaxed the internal precision tests of `testROH()`.  This was done
   in response to the CRAN farm lowering its precision on some hosts.


# Version 0.33.0 [2013-03-05]

## New Features

 * Added argument `typeOfWeights` to `estimateKappaByC1Density()` for
   `PairedPSCBS`, and hence indirectly to `estimateKappa()`.  The
   default is `typeOfWeights = "dhNbrOfLoci"`, which may give too much
   overall weight to very long segments causing the estimator to fail
   when there are only a few number of "C1 = 0" segments.  An
   alternative is to use `typeOfWeights = "sqrt(dhNbrOfLoci)"`.


# Version 0.32.6 [2013-03-04]

## Documentation

 * Updated the help usage section for all static methods.


# Version 0.32.5 [2013-02-09]

## Bug Fixes

 * `bootstrapTCNandDHByRegion()` for `PairedPSCBS` did not bootstrap
   from all available loci when calculating total CNs statistics, iff
   the segment had been called run-of-homozygosity (ROH). Internal
   validation tests caught this. Thanks to Oscar Rueda at the Cancer
   Research UK Cambridge Institute for reporting on this.

## Code Refactoring

 * Added a `VignetteBuilder` field to DESCRIPTION.


# Version 0.32.4 [2013-02-07]

## New Features

 * Improved some verbose outputs of `bootstrapTCNandDHByRegion()`.


# Version 0.32.3 [2013-02-05]

## New Features

 * Now `pruneByHClust()` drops any existing segment calls and quantile
   mean-level estimates.


# Version 0.32.2 [2013-02-01]

## New Features

 * Added `resetSegments()` for `AbstractCBS`, which drops extra
   segments columns (e.g. bootstrap statistics and calls) except those
   obtained from the segment algorithm.

## Documentation

 * Added a paragraph on `avgDH = "median"` to the **PSCBS** vignette's
   'Experimental' section.

## Code Refactoring

 * ROBUSTNESS: Now **aroma.light** is explicitly required in cases
   where it is needed.


# Version 0.32.1 [2013-02-01]

## Bug Fixes

 * `segmentByPairedPSCBS(..., avgDH = "median")` only worked for
   single- chromosome data.  Same for `avgTCN = "median"`.  Thanks
   Ritu Roy at UCSF for reporting on this.


# Version 0.32.0 [2013-01-16]

## New Features

 * Added arguments `avgTCN` and `avgDH` to `segmentByPairedPSCBS()`.

 * Now `updateMeans()` and `updateMeansTogether()` methods can
   estimate the mean levels either by the sample mean or the median.


# Version 0.31.0 [2013-01-05]

## Code Refactoring

 * CLEANUP: Now packages **R.methodsS3** and **R.oo** are only
   imported.

 * CLEANUP: Package no longer explicitly imports **digest**.


# Version 0.30.0 [2012-11-05]

## New Features

 * GENERALIZATION: Now `bootstrapTCNandDHByRegion()` works for more
   "flavors", e.g the default (`"tcn"`) used by
   `segmentByNonPairedPSCBS()`.


# Version 0.29.9 [2012-11-05]

## Documentation

 *  FIX: `example(segmentByNonPairedPSCBS)` was for the paired case.

## Code Refactoring

 * CRAN POLICY: Further speed up of examples such that they run faster
   with `R CMD check`.


# Version 0.29.8 [2012-11-04]

## Code Refactoring

 * CLEANUP: Replaced all `whichVector()` with `which()`, because the
   latter is now the fastest again.


# Version 0.29.7 [2012-11-03]

## Code Refactoring

 * Updated deprecated **ggplot2** functions in the RSP reports.


# Version 0.29.6 [2012-11-01]

## Code Refactoring

 * Bumped package dependencies.

 * CRAN POLICY: Made the examples run faster for `R CMD check`.


# Version 0.29.5 [2012-10-16]

## Bug Fixes

 * ROBUSTNESS: No longer passing `...` to `NextMethod()`, cf. R-devel
   thread "Do *not* pass '...' to NextMethod() - it'll do it for you;
   missing documentation, a bug or just me?" on Oct 16, 2012.


# Version 0.29.4 [2012-09-23]

## New Features

 * Now `plotTracks()` [and `plotTracksManyChromosomes()`] draws
   segment levels such that it is easier to see them even when they
   are overlapping.


# Version 0.29.3 [2012-09-21]

## New Features

 * SPEEDUP: By default `bootstrapTCNandDHByRegion()` for `PairedPSCBS`
   no longer do sanity checks within the bootstrap loop.  This
   significantly speed up the method.  To run checks, use argument
   `.debug = TRUE`.  In addition, the `callNnn()` methods that need to
   call this method, does it by decreasing the amount of verbose
   output substantially, which in turn speeds up the process a fair
   bit.

 * Now `getSegments(..., splitters = TRUE)` for `CBS` and `PSCBS`
   inserts NA rows wherever there is a "gap" between segments.  A
   "gap" is when two segments are not connected (zero distance).

 * ROBUSTNESS: Now `append()` for `CBS` and `PSCBS` drops column
   `length` from `knownSegments`, iff it exists.

 * Now `nbrOfChangePoints()` for `AbstractCBS` calculates only change
   points of connected neighboring segments.

## Bug Fixes

 * `seqOfSegmentsByDP()` for `AbstractCBS` would not handle empty
   segments, which could occur if `knownSegments` for instance
   included centromere gaps.

 * `segmentByCBS(... knownSegments)` could return segments for
   chromosome 0 even though it did not exist in the input data.


# Version 0.29.2 [2012-09-18]

## New Features

 * REPORT: Now `report()` for `AbstractCBS` looks for the RSP template
   in `templates/`, and as a backup in `templates,PSCBS/`.  If the
   latter does not exist, it is automatically created as a soft link
   to `templates/` of the **PSCBS** package.  This allows anyone to
   create their own customized copy (in `templates/`) of the default
   **PSCBS** RSP report.

 * REPORT: Now `report(fit, ..., rspTags)` for `AbstractCBS` looks for
   the RSP template named `<className>(,<rspTags>),report.tex.rsp`,
   where `<className>` is `class(fit)[1]` and argument `rspTags` is an
   optional comma-separated character string/vector.  This makes it
   possible to have different types of report for the same class of
   objects.

 * REPORT: Added argument `force` to `report()` for `AbstractCBS`.
   This will copy the RSP template files again, although they are
   already in `reports/` output directory.


# Version 0.29.1 [2012-09-15]

## New Features

 * Added argument `dropMissingCT` to `segmentByPairedPSCBS()`.


# Version 0.29.0 [2012-09-14]

## New Features

 * Added trial version of `pruneByDP()` for `AbstractCBS`.


# Version 0.28.6 [2012-09-13]

## New Features

 * Now `tileChromosomes()` also adjusts `knownSegments`.

 * Added argument `dropGaps` to `gapsToSegments()`.

 * Updated `all.equal()` for `AbstractCBS` to compare locus-level
   data, segments, and other fields.


# Version 0.28.5 [2012-09-13]

## New Features

 * SPEEDUP: Now `segmentByCBS(..., undo = +Inf)` returns much faster,
   which is possible because there is no need to identify new change
   points.


# Version 0.28.4 [2012-09-13]

## New Features

 * CONSISTENCY FIX: Changed the behavior of extreme values of argument
   `undo` to `segmentByCBS()` such that `undo = 0` (was `undo = +Inf`)
   now means that it will not ask `DNAcopy::segment()` to undo the
   segmentation, and such that `undo = +Inf` means that no
   changepoints will be identified. The latter case allows you to
   effectively skip the segmentation but still calculate all the CBS
   statistics across a set of known segments via `segmentByCBS(...,
   undo = +Inf, knownSegments = knownSegments)`.  Arguments `undoTCN`
   and `undoDH` to `segmentByPairedPSCBS()` are adjusted analogously.
   Corresponding system tests were added.


# Version 0.28.3 [2012-08-30]

## Code Refactoring

 * Updated code and Rd cross reference to use the **matrixStats**
   package for `weightedMedian()`, which used to be in
   **aroma.light**.


# Version 0.28.2 [2012-08-20]

## Bug Fixes

 * `segmentByNonPairedPSCBS()` forgot to specify namespace
   **aroma.light** when trying to call `findPeaksAndValleys()`.


# Version 0.28.1 [2012-08-15]

## Documentation

 * Minor grammatical corrections of the Paired PSCBS vignette.


# Version 0.28.0 [2012-07-22]

## New Features

 * Added argument `minLength` to `gapsToSegments()`.  The default is
   no longer to drop zero-length (`minLength == -1L`) segments,
   because if (and only if) such a segment contains a locus, then
   `segmentByNnn()` will currently generate an (internal) error.

## Bug Fixes

 * GENERALIZATION: Now `segmentByPairedPSCBS()` drops loci for which
   CT is missing (regardless of betaT). For instance, in rare cases
   when the reference (e.g. the normal) is missing, then it may be
   that CT is missing while betaT is not.


# Version 0.27.4 [2012-07-22]

## New Features

 * Now verbose output of `segmentByPairedPSCBS()` specifies region
   ranges with greater precision.


# Version 0.27.3 [2012-07-10]

## Documentation

 * Minor updates to the Paired PSCBS vignettes.

## Code Refactoring

 * CLEANUP: One redundancy tests relied on a non-critical function
   that will be removed in **R.utils** 1.16.0 (now in **R.devices**
   2.1.1).


# Version 0.27.2 [2012-07-08]

## Code Refactoring

 * Updated package dependencies.


# Version 0.27.1 [2012-07-02]

## New Features

 * Now we refer to "copy neutral" segments as "neutral TCN" segments
   with acronym `NTCN`.  The corresponding column in the segmentation
   results are labeled correspondingly.  The Paired PSCBS vignette was
   updated accordingly.


# Version 0.27.0 [2012-06-24]

## Documentation

 * Some grammar corrections of the `Paired PSCBS` vignette.

## Code Refactoring

 * (An update that should be ignored)


# Version 0.26.1 [2012-06-05]

## New Features

 * Now `segmentByCBS()` for `data frame`:s does a better job
   identifying the CN signals.


# Version 0.26.0 [2012-06-03]

## New Features

 * Now argument `delta` for `callCopyNeutralByTCNofAB()` of
   `PairedPSCBS` is calculated via `estimateDeltaCN()`, which
   estimates the width of the acceptance regions, used for calling
   copy neutral states, to be a function of the normal contamination.

## Documentation

 * Added details to the Paired PSCBS vignette on how to call segments
   that are copy neutral (typically diploid).


# Version 0.25.3 [2012-06-03]

## Bug Fixes

 * `all.equal(target, current)` for `CBS` objects would give an error
   if either `target` or `current` had zero segments.


# Version 0.25.2 [2012-05-30]

## New Features

 * Added `writeSegments()` for `DNAcopy` objects.

## Bug Fixes

 * `as.CNA()` for `DNAcopy` added incorrect chromosome splitters.

 * `as.CNA()` for `DNAcopy` would ignore argument `sample` and always
   return the first sample.


# Version 0.25.1 [2012-05-30]

## New Features

 * Now `callROH()` records parameter `deltaROH` in the results.

## Documentation

 * Added details to the Paired PSCBS vignette on how to tune the various
   callers.

## Bug Fixes

 * `callLOH(..., force = TRUE)` would append multiple `lohCall`
   columns, if called multiple times.


# Version 0.25.0 [2012-04-20]

## New Features

 * Added a trial (very much true) version of
   `segmentByNonPairedPSCBS()`.


# Version 0.24.0 [2012-04-20]

## New Features

 * Now it is possible to skip the DH segmentation in Paired PSCBS, i.e.
   `segmentByPairedPSCBS(..., flavor = "tcn")`.


# Version 0.23.2 [2012-04-20]

## Bug Fixes

 * `segmentByPairedPSCBS()` would throw `error in $<-.data.frame
   (*tmp*, "rho" ...` if some loci had unknown genomic positions.


# Version 0.23.1 [2012-04-20]

## New Features

 * Added RSP report for `CBS` objects (adopted from ditto for
   `PairedPSCBS`).

## Documentation

 * Updated the `Paired PSCBS` vignette.


# Version 0.23.0 [2012-03-20]

## Documentation

 * Added a package vignette.


# Version 0.22.2 [2012-02-29]

## Bug Fixes

 * `plotTracks(..., add = TRUE)` for `PairedPSCBS` would add TCNs when
   BAFs and DHs were intended.


# Version 0.22.1 [2012-02-28]

 * Updated package dependencies to **R.rsp** (>= 0.7.3) so that
   `report()` for `PairedPSCBS` no longer require non-public packages.

## New Features

 * Now it is possible to turn off usage of the alpha channel in plots
   generated by `report()`, which can be handy on systems where the
   default PNG device does not support the alpha channel.  Example:
   `setOption("PSCBS::report/useAlphaChannel", FALSE)`.


# Version 0.22.0 [2012-02-27]

## New Features

 * Added `report()` for `PairedPSCBS`.


# Version 0.21.0 [2012-02-27]

## New Features

 * Added argument `fields` to `getLocusData()` for `PairedPSCBS`.

 * Added `renameChromosomes()` to `AbstractCBS`.


# Version 0.20.0 [2012-02-26]

## New Features

 * Added alpha version of `callGainNeutralLoss()` for `PairedPSCBS`,
   which certainly will be updated in the future.  This caller is
   tested by the system tests.

 * Added `dropChangePoints()` for `AbstractCBS`.

## Bug Fixes

 * `extractSegments()` for `PairedPSCBS` would return incorrect row
   indices, more precisely, overlapping data chunks.

 * `bootstrapTCNandDHByRegion()` for `PairedPSCBS` would resample from
   a subset of the intended TCNs, iff the DH mean was non-finite while
   there were still heterozygous SNPs.  This introduced a bias in the
   estimates, which was neglectable for large segments, but for very
   small segments (a few loci) it could be relatively large.

## Software Quality

 * ROBUSTNESS: Added more sanity checks validating the correctness of
   what is returned by `extractSegments()` for `CBS` and
   `PairedPSCBS`.

## Code Refactoring

 * Added some internal utility functions for `PairedPSCBS` taken from
   the **aroma.cn** package.  Some of these may become public later,
   but for they should be considered internal.


# Version 0.19.8 [2012-02-23]

## Code Refactoring

 * ROBUSTNESS: Package now explicitly depends on **utils**.


# Version 0.19.7 [2012-02-22]

## Bug Fixes

 * `findLargeGaps()` did not handle missing values for argument
   `chromosome`.

 * `segmentByCBS(..., knownSegments = knownSegments)` would
   incorrectly throw a sanity-check exception if `knownSegments`
   contains a segment with `start` and `stop` positions being equal.

 * Argument `calls` of `plotTracks()` for `PairedPSCBS` was ignored if
   more than one chromosome was plotted.


# Version 0.19.6 [2012-01-24]

## New Features

 * ROBUSTNESS: Now `getCallStatistics()` for `CBS` asserts that calls
   have been made.  If not, an exception is thrown.


# Version 0.19.5 [2012-01-21]

## Documentation

 * Added details to the help of `callLOH()` and `callAB()` on the
   difference between (AB,LOH) = (TRUE,FALSE) and (AB,LOH) =
   (TRUE,NA).

 * Corrected some of verbose messages of
   `estimateDeltaLOHByMinC1ForNonAB()` for `PairedPSCBS` objects.


# Version 0.19.4 [2012-01-10]

## Code Refactoring

 * Now `example(segmentByPairedPSCBS)` and the system tests that are
   run by `R CMD check` are tuned to (by default) run much faster by
   segmenting using fewer data points and bootstrapping using fewer
   samples.  This update was done to meet the new CRAN policy.  By
   setting environment variable `_R_CHECK_FULL_` to `1` the full data
   set is used instead.


# Version 0.19.3 [2012-01-09]

## New Features

 * ROBUSTNESS: Now `extractSegments()` for `PairedPSCBS` gives an
   informative error message that it is not supported if CNs were
   segmented using flavor `"tcn,dh"`.

## Bug Fixes

 * `postsegmentTCN()` for `PairedPSCBS` could generate an invalid
   `tcnSegRows` matrix, where the indices for two consecutive segments
   would overlap, which is invalid.  Thanks to Minya Pu for reporting
   on failed sanity check related to this.


# Version 0.19.2 [2011-12-29]

 * ROBUSTNESS: Explicitly added **digest** to the list of suggested
   packages.


# Version 0.19.1 [2011-12-13]

## New Features

 * Added support for `callGainsAndLosses(..., method = "ucsf-dmad")`
   for `CBS` objects.


# Version 0.19.0 [2011-12-12]

## New Features

 * Added optional argument `indices` to `getLocusData()` to be able to
   retrieve the locus-level data as indexed by input data.

## Bug Fixes

 * `gapsToSegments()` gave invalid segments for chromosomes with more
   than one gap.  Now `gapsToSegments()` validates argument `gaps` and
   asserts that it returns non-overlapping segments.

## Documentation

 * Clarified in `help("segmentByCBS")` how missing values are dealt
   with.


# Version 0.18.2 [2011-12-07]

## New Features

 * Now `plotTracks()` for `CBS` always returns an invisible object.

## Bug Fixes

 * `pruneBySdUndo()` for `CBS` did not work with more than one array.


# Version 0.18.1 [2011-12-03]

## New Features

 * Added `drawChangePoints()` for `AbstractCBS`.

 * Now `pruneByHClust()` for `AbstractCBS` updates the segment means.

 * Added `writeSegments()` for `PSCBS` object.

 * Now `print()` for `AbstractCBS` returns `getSegments(..., simplify
   = TRUE)`.

 * Added argument `simplify` to `getSegments()`.

 * Added arguments `name`, `tags` and `exts` to `writeSegments()` and
   `writeLocusData()` and dropped `filename`.


# Version 0.18.0 [2011-11-28]

## New Features

 * Added `pruneByHClust()` for `AbstractCBS`, with implementation for
   `CBS` and `PairedPSCBS`.

 * `extractCNs()` for `CBS` would not return a matrix but a
   data.frame.

## Bug Fixes

 * `extractTotalCNs()` for `CBS` would give an error.


# Version 0.17.4 [2011-11-26]

## New Features

 * Added argument `updateMeans = TRUE` to `callROH()` for
   `PairedPSCBS`.

 * Now `bootstrapTCNandDHByRegion()` for `PairedPSCBS` preserves NAs
   for DH and (C1,C2) quantiles, if the DH mean level is NA, which can
   happen when a segment is called ROH.  This also makes sure that a
   segment called ROH will not be called AB.

 * An internal sanity check of `bootstrapTCNandDHByRegion()` for
   `PairedPSCBS` would give an error if DH mean levels had been set to
   NA for segments called ROH.


# Version 0.17.3 [2011-11-24]

## New Features

 * Added `callSegmentationOutliers()` and `dropSegmentationOutliers()`
   for data frames.

 * CLEANUP: Renamed field `position` of the example data to `x`. This
   helps us clean up some of the examples.

## Bug Fixes

 * `bootstrapTCNandDHByRegion()` for `PairedPSCBS` would give an
   error, if a segment did not have any TCN signals, which can occur
   when known segments are specified for Paired PSCBS.


# Version 0.17.2 [2011-11-22]

## New Features

 * Added `findLargeGaps()` and `gapsToSegments()`.


# Version 0.17.1 [2011-11-21]

## Bug Fixes

 * The internal sanity check of `testROH()` on weights was slightly
   too conservative (required to high precision) when it came to
   asserting that the sum of the weights equals one.

 * `resegment()` for `PairedPSCBS` called `segmentByCBS()` instead of
   `segmentByPairedPSCBS()`.


# Version 0.17.0 [2011-11-19]

## New Features

 * Now it is possible to run Paired PSCBS (without TumorBoost) when only
   genotypes but not BAFs are available for the matched normal.


# Version 0.16.3 [2011-11-17]

## New Features

 * Added `resegment()` for `CBS` and `PairedPSCBS` for easy
   resegmentation.

 * Adjusted `segmentByCBS()` such that it can handle `knownSegments`
   with chromosome boundaries given as -Inf and +Inf.

 * Now argument `mar` for `plotTracks()` defaults to NULL.

 * ROBUSTNESS: Added redundancy tests for `segmentByCBS()` and
   `segmentByPairedPSCBS()` with argument `knownSegments`.

 * ROBUSTNESS: Now `segmentByCBS()` does more validation of
   `knownSegments`.

 * FIX: `extractRegions()` for `AbstractCBS` would also show verbose
   output.

## Bug Fixes

 * Now argument/parameter `seed` is correctly preserved by
   `segmentByCBS()`. So is `tbn` for `segmentByPairedPSCBS()`.

 * `segmentByPairedPSCBS()` would give an error when trying to segment
   DH if the TCN segment contains no data points, which could happen
   if `knownSegments` specifies an empty segment, e.g. centromere.

 * `extractSegments()` for `CBS` would throw an error when there were
   multiple chromosomes.


# Version 0.16.2 [2011-11-16]

## New Features

 * Now `segmentByCBS(..., w)` stores weights `w`, if given, in the
   locus-level data table of the returned CBS object.

 * Added `pruneBySdUndo()` for `CBS`, which does what `undo.splits =
   "sdundo"` for `DNA::segment()`, but on the already segmented
   results.

 * Now `updateMeans()` uses locus-specific weights, iff available.

 * Added `updateBoundaries()` for `CBS` to update (start,stop) per
   segment.

 * CORRECTNESS: Now `updateMeans()` for `CBS` identifies loci via
   internal `segRows` field and no longer by locations of segment
   boundaries, which gave slightly incorrect estimates for "tied"
   loci.


# Version 0.16.1 [2011-11-15]

## New Features

 * Now more segmentation parameters are stored in the `CBS` object.

 * SPEEDUP: Now `segmentByCBS()` will use memoization to retrieve so
   called "sequential boundaries for early stopping", iff any of the
   `DNAcopy::segment()` arguments `alpha`, `nperm` and `eta` are
   specified.  See also `DNAcopy::getbdry()`.

 * Added `method = "DNAcopy"` to `estimateStandardDeviation()` for
   `CBS`, which estimates the std. dev. using
   `DNAcopy:::trimmed.variance()`.

## Bug Fixes

 * `extractSegments()` for `CBS` would throw an error, because in most
   cases it would created a corrupt internal `segRows` field.


# Version 0.16.0 [2011-11-12]

## New Features

 * Added argument `oma` and `mar` to `plotTracksManyChromosomes()` for
   `PairedPSCBS` for setting graphical parameters when `add = FALSE`.

 * Added `callROH()`.

 * Added arguments `from` and `adjustFor` to `updateMeans()`.


# Version 0.15.5 [2011-11-04]

## Bug Fixes

 * `extractSegment()` for `AbstractCBS` would give an error, because
   it called itself instead of `extractSegments()`.


# Version 0.15.4 [2011-10-30]

## New Features

 * Added `save()` and `load()` methods to `AbstractCBS`, which are
   wrappers for `saveObject()` and `loadObject()` that assert the
   correct class structure.  Also, the `load()` method will
   automatically update the class hierarchy for `CBS` and
   `PairedPSCBS` objects that were saved before adding class
   `AbstractCBS`.


# Version 0.15.3 [2011-10-23]

## Bug Fixes

 * `callAmplifications()` for `CBS` generated an error, if more than
   one chromosome were called.

 * The length of a segment must be defined as 'end-start' and not
   'end-start+1' so that the the total length of all segments adds up
   correctly.

 * `highlightArmCalls()` for `CBS` did not handle empty chromosomes.

 * `getCallStatisticsByArms()` for `CBS` threw a error if argument
   `genomeData` did not contain exactly the same chromosomes as in the
   `CBS` object.


# Version 0.15.2 [2011-10-21]

## New Features

 * Added `mergeThreeSegments()` for `AbstractCBS`.

## Bug Fixes

 * Recent updates caused `segmentByPairedPSCBS(data)` not to work when
   `data` is a data frame.


# Version 0.15.1 [2011-10-21]

## New Features

 * By setting `start` and `end` to NAs in `knownSegments` (chromosome
  must still be specified), it is possible to insert an empty segment
  that disconnects the two flanking segments, e.g. centromere and the
  two arms.


# Version 0.15.0 [2011-10-20]

## New Features

 * Added support for specifying priorly known segments, such as
   chromosome arms and centromeres, in `segmentByCBS()` via argument
   `knownSegments`.

## Bug Fixes

 * CLEANUP: Dropped a stray debug output message in
   `segmentByPairedPSCBS()`.


# Version 0.14.3 [2011-10-17]

## New Features

 * Added argument `asMissing` to `dropRegions()` for `AbstractCBS`.


# Version 0.14.2 [2011-10-16]

## New Features

 * Implemented `extractCNs()` for `CBS` and `PairedPSCBS`.

 * Added `extractTotalCNs()` for `CBS`.


# Version 0.14.1 [2011-10-14]

## New Features

 * Added implementation of `extractRegions()` for `AbstractCBS`, which
   utilizes `extractSegments()`.

 * Added abstract `extractSegments()` and `extractSegment()` for
   `AbstractCBS`.

 * Now `extractTCNAndDHs()` for `PairedPSCBS` passes `...` to
   `getSegments()`.


# Version 0.14.0 [2011-10-10]

## Code Refactoring

 * CLEANUP: Harmonization of several method names.

 * CLEANUP: Internal restructuring of the source code files.


# Version 0.13.5 [2011-10-10]

## New Features

 * Added `dropChangePoint()` for `AbstractCBS`, which is just a "name
   wrapper" for `mergeTwoSegments()`.

 * Added `dropRegion()` and `dropRegions()` for `AbstractPSCBS`, where
   the former is a wrapper for the latter `dropRegions()`.

 * Added `updateMeans()` and `mergeTwoSegments()` for `CBS` in
   addition already available `PairedPSCBS` versions.

 * Relabeled column `id` to `sampleName` returned by `getSegments()`.

## Bug Fixes

 * For so called "splitter" rows, not all columns returned by
   `getSegments()` of `CBS` were missing values.

 * The object returned by `as.CBS()` for `DNAcopy` did not have the
   correct class hierarchy.

## Code Refactoring

 * ROBUSTNESS: Now using `getSegments()` everywhere possible.


# Version 0.13.4 [2011-10-08]

## New Features

 * Added `all.equal()` for `AbstractCBS`, which does not compare
   attributes.

 * Added optional argument `regions` to `getCallStatistics()` for
   `CBS` in order to calculate call statistics on subsets of
   chromosomes, e.g. chromosome arms.

 * Added `drawChromosomes()` for `CBS`.

 * Added `getCallStatisticsByArms()`, `callArms()`, and
   `highlightArmCalls()` for `CBS` objects.

## Code Refactoring

 * Now internal `getChromosomeRanges()` of `CBS` returns a data.frame
   instead of a matrix, and first column is now `chromosome`.


# Version 0.13.3 [2011-10-03]

## New Features

 * GENERALIZATION: Now `segmentByCBS()` and `segmentByPairedPSCBS()`
   also accepts a data.frame of locus-level data with column names
   matching the locus-level arguments accepted by the corresponding
   method.

 * GENERALIZATION: Now all segmentation result classes (`CBS` and
   `PSCBS`) inherits from the `AbstractCBS` class, which provides
   methods such as `getSampleName()`, `getChromosomes()` and
   `getSegments()`.

## Documentation

 * Added lots of more help pages.

## Code Refactoring

 * CLEANUP: Dropped empty `callSegments()` for `PairedPSCBS`.


# Version 0.13.2 [2011-09-30]

## New Features

 * GENERALIZATION: Now `drawLevels()` for `PairedPSCBS` allows for
   drawing segmentation results in 'betaT' space.

## Bug Fixes

 * `plotTracks2(..., panels = "dh")` gave an error due to a forgotten
   assignment.


# Version 0.13.1 [2011-09-06]

## New Features

 * Added formal class `CBS`, which holds the segmentation results
   returned by `segmentByCBS()`.  Several methods are available for
   `CBS` objects, e.g.  `nbrOfLoci()`, `nbrOfSegments()`,
   `nbrOfChromosomes()`, `getChromosomes()`,
   `estimateStandardDeviation()`, etc.

 * Now `segmentByCBS()` always returns a `CBS` object.  To coerce to a
   `DNAcopy` object (as defined in the `DNAcopy` class) use
   `as.DNAcopy()`.

 * Added coerce methods `as.DNAcopy()` for `CBS` objects and
   `as.CBS()` for `DNAcopy` objects.


# Version 0.13.0 [2011-09-01]

## New Features

 * GENERALIZATION: Now `segmentByCBS()` can process multiple
   chromosomes.

 * Added `append()` for `CBS` objects.

## Bug Fixes

 * Internal methods `plotTracksManyChromosomes()` and
   `tileChromosomes()` for `CBS` did not work at all and therefore
   neither `plotTracks()` for `CBS` with more than one chromosome.


# Version 0.12.2 [2011-08-27]

## Code Refactoring

 * CLEANUP: Now `R CMD check` is no longer giving a note that the
   package loads package **DNAcopy** in `.onAttach()`.


# Version 0.12.1 [2011-08-08]

## Bug Fixes

 * If `dropSegmentationOutliers()` would drop an outlier next to a
   change point, such that the total copy-number signal becomes NA,
   then the sanity checks that TCN segments always overlaps DH
   segments would fail.  Now the sanity checks are aware of this
   special case.  These sanity checks were moved from
   `bootstrapTCNandDHByRegion()` to `segmentByPairedPSCBS()`.  Thanks
   Christine To at University of Toronto for reporting on this.


# Version 0.12.0 [2011-07-23]

## Bug Fixes

 * Recently R devel automatically adds a namespace to a package, if
   missing.  This caused some of the **PSCBS** examples to throw an
   exception related to incorrect dispatching of `cat()`.

## Code Refactoring

 * Added a namespace to the package, which will be more or less a
   requirement in the next major release of R.


# Version 0.11.7 [2011-07-15]

## Documentation

 * Added a section to `help("segmentByPairedPSCBS")` on the importance
   of doing a whole-genome PSCBS segmentations if calling AB and LOH
   states afterward.

 * Made it more clear in `help("segmentByPairedPSCBS")` that arguments
   `betaT`, `betaN` and `muN` may contain NAs for non-polymorphic
   loci.


# Version 0.11.6 [2011-07-14]

## Bug Fixes

 * ROBUSTNESS: In some cases, the segmentation table would contain
   column names with incorrect capitalization, e.g. `tcnnbrOfLoci`
   instead of `tcnNbrOfLoci`.  This would cause several downstream
   methods to give an error.  The reason for this is that the
   **Hmisc** package, if loaded after **R.utils**, overrides
   `capitalize()` in **R.utils** with another (buggy?) `capitalize()`
   function.  To avoid this, we now everywhere specify explicitly that
   we want the one in **R.utils**. Thanks Christine To at University
   of Toronto for reporting on this.


# Version 0.11.5 [2011-07-10]

## Bug Fixes

 * `tileChromosomes()` for `PairedPSCBS` was still assuming the old
   naming convention of column names.  This caused `plotTracks()` to
   throw an exception when plotting multiple chromosomes.

## Code Refactoring

 * ROBUSTNESS: Fixed partial argument matchings in `arrowsC1C2()` and
   `arrowsDeltaC1C2()` for `PairedPSCBS`.


# Version 0.11.4 [2011-07-07]

## New Features

 * GENERALIZATION: Now the internal estimator function that
   `estimateDeltaLOH()` uses returns -Inf if all segments are called
   AB, instead of throwing an exception.  This will in turn make
   `callLOH()` call all segments to be non-LOH.

## Documentation

 * Removed obsolete references to the R-forge repository.

## Bug Fixes

 * Consecutive calls to `callAB(..., force = TRUE)` would append
   additional `abCall` columns to the segmentation table instead of
   replacing existing calls.


# Version 0.11.3 [2011-07-06]

## New Features

 * ROBUSTNESS: Added a sanity check to
   `estimateDeltaLOHByMinC1AtNonAB()` for `PairedPSCBS` object. The
   test asserts that there exist segments that are not in allelic
   balance, which are needed in order to estimate DeltaLOH.

## Documentation

 * The description of argument `chromosome` for
   `segmentByPairedPSCBS()` did not describe how to segment multiple
   chromosomes in one call.


# Version 0.11.2 [2011-07-05]

## Bug Fixes

 * Output fields `tcnNbrOfSNPs` and `tcnNbrOfHets` were mistakenly
   labeled as `tcnNbrOr...`.  Thanks Christine Ho at UC Berkeley for
   reporting on this.


# Version 0.11.1 [2011-06-28]

## Documentation

 * Clarified that argument `CT` should be tumor copy number ratios
   relative to the normal.

 * Added Rd help for `as.data.frame()` of `PairedPSCBS`.


# Version 0.11.0 [2011-06-14]

## Significant Changes

 * Renamed all column names of returned data frames such that they
   follow the camelCase naming conventions in addition to be somewhat
   shorter too.

## New Features

 * GENERALIZATION: Added argument `columnNamesFlavor` to
   `segmentByCBS()`.


# Version 0.10.2 [2011-06-07]

## Code Refactoring

 * CLEANUP: Cleaned up the `example()`:s.

 * Added more `biocViews` categories to DESCRIPTION.


# Version 0.10.1 [2011-05-31]

## New Features

 * GENERALIZATION: The package can now be _installed_ without the
   **DNAcopy** package being installed.  If package is loaded without
   **DNAcopy** installed, an informative message will explain how to
   install it.

 * Added `installDNAcopy()`, which will install **DNAcopy** from
   Bioconductor.

## Code Refactoring

 * ROBUSTNESS: Now all **DNAcopy** functions are called as
   `DNAcopy::nnn()`.


# Version 0.10.0 [2011-05-29]

## Significant Changes

 * Renamed package to **PSCBS** (from **psCBS**).

 * Renamed all arguments, variables, and functions referring to `tau`
   to refer to `delta` reflecting the notation of the Paired PSCBS
   paper.

 * Renamed options, example code and help pages to reflect new package
   name.

## Documentation

 * Updated references in help pages.

 * Now the paired PSCBS is formally referred to as 'Paired PSCBS'.


# Version 0.9.54 [2011-04-27]

## New Features

 * Added argument `maxC` to `estimateTauLOHByMinC1ForNonAB()`.


# Version 0.9.53 [2011-04-14]

 * Added argument `max` to `estimateTauAB()` and `estimateTauLOH()`.


# Version 0.9.52 [2011-04-14]

## Bug Fixes

 * Argument `minSize` of `callAB()` and `callLOH()` had no effect.


# Version 0.9.51 [2011-04-12]

## New Features

 * Added argument `minSize` to `callAB()` and `callLOH()` for
   `PairedPSCBS`.

 * Now the a conflicting call in `callLOH()`/`callAB()` with argument
   `xorCalls = TRUE` is set to NA to contrast it from a FALSE call.


# Version 0.9.50 [2011-04-12]

## New Features

 * Added argument `xorCalls` to `callLOH()` and `callAB()` for
   `PairedPSCBS`.  When TRUE (the default), a segment that is already
   called AB will never be called LOH, and vice versa.


# Version 0.9.49 [2011-04-11]

## New Features

 * Updated `estimateTauABBySmallDH()` for `PairedPSCBS` to use a
   "symmetric" quantile estimator.

 * Added argument `midpoint` to `estimateTauLOHByMinC1AtNonAB()`.

## Bug Fixes

 * The recent `callLOH()` would not store the LOH calls.


# Version 0.9.48 [2011-04-10]

## New Features

 * Added `callLOH()` for `PairedPSCBS`, which in turn calls auxiliary
   methods.

 * Added `estimateTauLOH()` for `PairedPSCBS`, which in turn calls
   auxiliary methods.

 * Now `callAB(..., force = FALSE)` skips the caller if
   allelic-balance calls already exist.

## Documentation

 * Update the example for `segmentByPairedPSCBS()` to reflect the
   restructured AB and LOH callers.


# Version 0.9.47 [2011-04-08]

## New Features

 * Added `estimateTauABBySmallDH()`.

 * Added internal `weightedQuantile()`.

## Documentation

 * Added help pages for more methods.

## Code Refactoring

 * CLEANUP: Started to restructure the source code files.


# Version 0.9.46 [2011-04-08]

## Bug Fixes

 * `postsegmentTCN()` for `PairedPSCBS` could generate an invalid
   `tcnSegRows` matrix, where the indices for two consecutive segments
   would overlap, which is invalid.  This was caught with real data,
   but it seems to have required a very rare combination of data in
   order for it to occur.


# Version 0.9.45 [2011-04-05]

## Bug Fixes

 * `estimateHighDHQuantileAtAB()` for `PairedPSCBS` would throw an
   error on an undefined `trim` if verbose output was used.


# Version 0.9.44 [2011-02-18]

## New Features

 * Added `estimateHighDHQuantileAtAB()` for `PairedPSCBS`.


# Version 0.9.43 [2011-02-06]

## Bug Fixes

 * `plotTracks2()` queried non-existing argument `tracks`.


# Version 0.9.42 [2011-02-03]

## New Features

 * Added `estimateKappa()` for estimating the normal contamination.


# Version 0.9.41 [2011-02-02]

## Significant Changes

 * Updated default for `tauAB` of `callABandHighAI()` and
   `callABandLowC1()` to be estimated from data using
   `estimateTauAB()`.

## New Features

 * Added argument `tauTCN` to `estimateTauAB()`.


# Version 0.9.40 [2011-01-27]

## New Features

 * Added argument `flavor` to `estimateTauAB()` for estimating the AB
   threshold using alternative methods.


# Version 0.9.39 [2011-01-19]

## New Features

 * Added trial version of new `plotTracks2()`, which will later
   replace `plotTracks()`.  Currently it only works for single
   chromosomes.

 * Added support functions, e.g. `updateMeans()`.


# Version 0.9.38 [2011-01-18]

## New Features

 * Added arguments `changepoints` and `col` to `plotTracks()` for
   `PairedPSCBS`.

 * Now `plotTracks(..., add = FALSE)` for `PairedPSCBS` only sets up
   subplots if argument `tracks` specifies more than one panel.

## Documentation

 * Documented more `plotTracks()` arguments for `PairedPSCBS`.

## Bug Fixes

 * Now `plotTracks(..., add = TRUE)` for `PairedPSCBS` plots to the
   current figure/panel.


# Version 0.9.37 [2011-01-18]

## Bug Fixes

 * `tcnSegRows` and `dhSegRows` where not updated by
   `extractByRegions()` for `PairedPSCBS`.


# Version 0.9.36 [2011-01-14]

## New Features

 * Added `estimateTauAB()` for estimating the tauAB tuning parameter
   when calling segments in allelic balance.  Updated
   `example(segmentByPairedPSCBS)` to illustrate how to use it.

 * Added `extractByRegions()` for `PairedPSCBS`.


# Version 0.9.35 [2011-01-12]

## New Features

 * Now `postsegmentTCN(..., force = TRUE)` for `PairedPSCBS` also
   updates the TCN estimates even for segments where the DH
   segmentation did not find any additional change points.


# Version 0.9.34 [2010-12-09]

## Bug Fixes

 * When there were multiple chromosomes processed by
   `segmentByPairedPSCBS()`, then the returned data object would
   contain `betaT` identical to `betaTN`.


# Version 0.9.33 [2010-12-07]

## New Features

 * Added `callLowC1ByC1()` and `callABandLowC1()`.


# Version 0.9.32 [2010-12-03]

## Bug Fixes

 * In rare cases the bootstrap sanity checks can indeed produce an
   invalid `range`, more precisely where (`range[,2]` >= `range[,1]`)
   is not true.  This can happen if there is no variation in the
   bootstrap estimates.  Because of this we allow for some tolerance.


# Version 0.9.31 [2010-12-02]

## New Features

 * Added option `psCBS/sanityChecks/tolerance` for specifying the
   tolerance of some internal sanity checks.


# Version 0.9.30 [2010-12-01]

## Code Refactoring

 * Rewrote all code dealing with the identification of loci belong to
   segments.  The code is now utilizing the `segRows` element returned
   by `DNAcopy::segment()`.  Lots of the code was rewritten and
   therefore completely new bugs may have been introduced.


# Version 0.9.25 [2010-11-30]

## Bug Fixes

 * Argument `flavor` of `segmentByPairedPSCBS()` would be ignored if
   multiple chromosomes were segmented.

 * `extractByChromosome()` for `PSCBS` would call it self instead of
   `extractByChromosomes()`.


# Version 0.9.24 [2010-11-28]

## Bug Fixes

 * `postsegmentTCN()` did not handle loci with the same positions and
   that are split in two different segments.  It also did not exclude
   loci with missing values.


# Version 0.9.23 [2010-11-28]

## Bug Fixes

 * The algorithm in `segmentByCBS()` that infers which loci (of the
   ones share the same genomic positions) that should be exclude from
   each segment did not take missing signals into account.

 * Iff argument `chromosome` to `segmentByPairedPSCBS()` was of length
   greater than one and specified exactly one unique chromosome, then
   exception `Number of elements in argument 'chromosome' should be
   exactly 8712 not 86209 value(s)` would be thrown.


# Version 0.9.22 [2010-11-27]

## Bug Fixes

 * `bootstrapTCNandDHByRegion()` would incorrectly include
   non-polymorphic loci in the set of homozygous SNPs during
   resampling.

 * `segmentByPairedPSCBS()` would not accept missing values in
   argument `chromosome`.


# Version 0.9.21 [2010-11-27]

## New Features

 * Now arguments `...` of `segmentByPairedPSCBS()` are passed to the
   two `segmentByCBS()` calls.

 * Added `callSegmentationOutliers()`, which can be used to identify
   single-locus outliers that have a genomic signal that is clearly
   outside the expected range.  The `dropSegmentationOutliers()` sets
   locus outliers detected by this method to missing values. This is
   useful for excluding total copy-number outliers that otherwise can
   have a dramatic impact on the non-robust CBS method.


# Version 0.9.20 [2010-11-26]

## New Features

 * Added optional argument `chromosomes` to `plotTracks()` to plot a
   subset of all chromosomes.

 * Added `extractByChromosomes()` for `PSCBS`.

 * Now the default confidence intervals for `plotTracks()` is
   (0.05,0.95), if existing.

 * Now all call functions estimate symmetric bootstrap quantiles for
   convenience of plotting confidence intervals.

## Bug Fixes

 * `callABandHighAI()` for `PairedPSCBS` used the old DH-only
   bootstrap method.

 * The statistical sanity checks of the bootstrap estimates would give
   an error when only single-sided bootstrap confidence interval was
   calculated.

 * The call functions, for instance `callABandHighAI()`, would throw
   `Error in quantile.default(x, probs = alpha) : missing values and
   NaN's not allowed if 'na.rm' is FALSE`, unless
   `bootstrapTCNandDHByRegion()` was run before.


# Version 0.9.19 [2010-11-23]

## New Features

 * ROBUSTNESS: Added more sanity checks to
   `bootstrapTCNandDHByRegion()`.

 * WORKAROUND: The precision of the mean levels of
   `DNAcopy::segment()` is not great enough to always compare it to
   that of R's estimates.

## Bug Fixes

 * `bootstrapTCNandDHByRegion()` would give an error if there was only
   one segment.

 * `segmentByPairedPSCBS()` and `bootstrapTCNandDHByRegion()` would
   not subset the correct set of DH signals if there were some missing
   values in TCN.


# Version 0.9.18 [2010-11-22]

## New Features

 * Added argument `calls` to `plotTracks()` for highlighting called
   regions.

 * Updated `callAllelicBalanceByDH()` and
   `callExtremeAllelicImbalanceByDH()` to utilize
   `bootstrapTCNandDHByRegion()`.

 * ROBUSTNESS: Now `drawConfidenceBands()` of `PairedPSCBS` silently
   does nothing if the requested bootstrap quantiles are available.

## Bug Fixes

 * `bootstrapTCNandDHByRegion()` for `PairedPSCBS` would not correctly
   detect if bootstrap results are already available.


# Version 0.9.17 [2010-11-21]

## New Features

 * Now `plotTracks()` supports tracks `"tcn,c1"`, `"tcn,c2"`, and
   `"c1,c2"` too.

 * Added support for flavor `"tcn&dh"` in `segmentByPairedPSCBS()`,
   which contrary to `"tcn,dh"` enforces TCN and DH to have the same
   change points.  The default flavor is now `"tcn&dh"`.

 * Added argument `xlim` to `plotTracks()` making it possible to zoom
   in.


# Version 0.9.16 [2010-11-21]

## New Features

 * Now `joinSegments = TRUE` is the default for `segmentByCBS()` and
   `segmentByPairedPSCBS()`.

 * Added argument `quantiles` to `plotTracks()`, which if specified
   draws confidence bands previously estimated from bootstrapping.

 * Added `drawConfidenceBands()` for `PairedPSCBS`.

 * Added `bootstrapTCNandDHByRegion()` for `PairedPSCBS`.

 * Added standalone `joinSegments()` for `CBS` results.

 * Now `segmentByPairedPSCBS()` also returns minor and major copy
   numbers for each segment.


# Version 0.9.15 [2010-11-21]

## New Features

 * Adjusted `postsegmentTCN()` such that the updated TCN segment
   boundaries are the maximum of the DH segment and the support by the
   loci.  This means that `postsegmentTCN()` will work as expected
   both when signals where segmented with `joinSegments` being TRUE or
   FALSE.

 * Updated `plotTracks()` for `PairedPSCBS` such that the TCN
   segmentation is colored 'purple' and the DH segmentation 'orange'
   for TCN and DH only tracks.


# Version 0.9.14 [2010-11-20]

## New Features

 * Now it is possible to specify the boundaries of the regions to be
   segmented as known change points via argument `knownCPs`.

 * Added argument `joinSegments` to `segmentByCBS()` and
   `segmentByPairedPSCBS()` in order to specify if neighboring
   segments should be joined or not.

 * Now `segmentByCBS()` and `segmentByPairedPSCBS()` allow for unknown
   genomic positions as well as missing total CN signals.


# Version 0.9.13 [2010-11-19]

## New Features

 * Added argument `joinSegments` to `segmentByCBS()` in order to
   specify if neighboring segments should be joined or not.


# Version 0.9.12 [2010-11-19]

## New Features

 * Added `plotTracks()` and `drawLevels()` etc. to CBS results.

 * Now `segmentByCBS()` allows for unknown genomic positions.

 * Now `segmentByCBS()` allows for missing signals.

 * Added argument `preserveOrder` to `segmentByCBS()`.  If TRUE, then
   the loci in the returned `data` object are ordered as the input
   data, otherwise it is ordered along the genome.


# Version 0.9.11 [2010-11-16]

## New Features

 * Now the `data` object returned by `segmentByCBS()` contains field
   `index` if and only if the loci had to be reorder along the genome.

## Documentation

 * Added more details, references to papers, and cross links to other
   functions to the help pages.

## Bug Fixes

 * In the rare cases where two loci at the same positions are split up
   into two neighboring segments, then `segmentByPairedPSCBS()` would
   fail to infer which they were if and only if the loci were not
   ordered along the genome.  This could happen with for instance
   Affymetrix `GenomeWideSNP_6` data.


# Version 0.9.10 [2010-11-09]

## New Features

 * Added argument `cex = 1` to `plotTracks()`.

## Bug Fixes

 * It was not possible to plot BAF tracks with `plotTracks()`.


# Version 0.9.9 [2010-11-05]

## Bug Fixes

 * `segmentByCBS()` tried to pass non-existing argument `undo.split`
   to `DNAcopy::segment()`.  It should be `undo.splits`.


# Version 0.9.8 [2010-11-04]

## Bug Fixes

 * There was a stray/debug `stop()` statement left in
   `segmentByPairedPSCBS()` causing an "error" in the rare case when
   loci that have the same physical locations are split into two
   different segments.


# Version 0.9.7 [2010-11-03]

 * ROBUSTNESS: Now `bootstrapDHByRegion()` uses `resample()` of
   **R.utils**.

## Bug Fixes

 * `bootstrapDHByRegion()` did not sample from the correct unit(s)
   when there was only one DH signal.


# Version 0.9.6 [2010-11-02]

## New Features

 * Added arguments `undoTCN` and `undoDH` to `segmentByPairedPSCBS()`.

 * Added argument `undo` to `segmentByCBS()`, which corresponds to
   `undo.splits = "sdundo"` and `undo.SD = undo`, if undo < +Inf.

## Bug Fixes

 * Arguments `alphaTCN` and `alphaDH` of `segmentByPairedPSCBS()` were
   not used when more than one chromosome were segmented.


# Version 0.9.5 [2010-11-01]

## New Features

 * Added arguments `alphaAB` and `alphaHighAI` to `callABandHighAI()`.

## Bug Fixes

 * `bootstrapDHByRegion()` would give an error if only a single
   quantile was requested.

 * `bootstrapDHByRegion()` would give `Error in if (nbrOfUnits >
   segJJ[, "dh.num.mark"]) { : missing value where TRUE/FALSE needed`
   when `dh.num.mark` was NA.


# Version 0.9.4 [2010-10-25]

## New Features

 * Now the default is a 95% confidence interval for calls.

 * Now `segmentByCBS()` also returns element `lociNotPartOfSegment`,
   if there are segments that share end points, which can happen if a
   change point is called in middle of a set of loci that have the
   same genomic positions.  In such cases, `lociNotPartOfSegment`
   specifies which loci are *not* part of which segment.  Then by
   identifying the loci that are within a segment by their positions
   and excluding any of the above, one knows exactly which loci CBS
   included in each segment.

## Bug Fixes

 * Now `bootstrapDHByRegion()` for `PairedPSCBS` handles the rare case
   when markers with the same positions are split in two different
   segments.

 * Now the correct set of loci are extracted from each TCN segment, in
   the rare case that two neighboring TCN segments have the same end
   points.


# Version 0.9.3 [2010-10-25]

## New Features

 * Added argument `ciRange` to `callAllelicBalance()` and
   `callExtremeAllelicImbalance()`.

## Bug Fixes

 * `bootstrapDHByRegion()` for `PairedPSCBS` would bootstrap from the
   incorrect set of loci when the DH region contained only one locus.

 * `bootstrapDHByRegion()` for `PairedPSCBS` would bootstrap from the
   incorrect set of loci if more than one chromosome was available.


# Version 0.9.2 [2010-10-24]

## Bug Fixes

 * `plotTracks()` would give `Error: object 'nbrOfLoci' not found` for
   whole-genome plots.


# Version 0.9.1 [2010-10-20]

## New Features

 * Now `plotTracks()` can plot whole-genome data.


# Version 0.9.0 [2010-10-18]

## New Features

 * Added arguments `alphaTCN` and `alphaDH` to
   `segmentByPairedPSCBS()` with defaults according to the paper.


# Version 0.8.3 [2010-10-18]

## New Features

 * Now `segmentByPairedPSCBS()` can segment multiple chromosomes.


# Version 0.8.2 [2010-10-17]

## New Features

 * Added argument `tbn` to `segmentByPairedPSCBS()` specifying whether
   TumorBoostNormalization should be applied or not.


# Version 0.8.1 [2010-10-10]

## Significant Changes

 * The default for `segmentByPairedPSCBS()` is now to segment TCN on
   the original scale, not the `sqrt()`.

## New Features

 * Added `plotTracks()` for `PairedPSCBS`.


# Version 0.8.0 [2010-10-06]

## Code Refactoring

 * CLEANUP: Removed all old code, native code and help pages.


# Version 0.7.8 [2010-10-03]

## New Features

 * Added optional argument `chromosome` to `segmentByCBS()`. Note that
   at this point it is only used for annotating the results; it can
   not be used to segmented multiple chromosomes at ones.


# Version 0.7.7 [2010-09-26]

## New Features

 * Now `subsetBySegments()` and `postsegmentTCN()` for `PairedPSCBS`
   handles multiple chromosomes.


# Version 0.7.6 [2010-09-24]

## New Features

 * Added support to annotating and subsetting also by chromosomes, as
   well as appending segmentation results from different chromosomes
   together.


# Version 0.7.5 [2010-09-21]

## New Features

 * Added `postsegmentTCN()` for `PairedPSCBS`, which updates the TCN
   segment start and ends, estimates and counts given the DH segments.


# Version 0.7.4 [2010-09-18]

## New Features

 * Added argument `chromosome` to `segmentByPairedPSCBS()`, which, if
   given, adds a chromosome column to the data and segmentation
   results.

## Bug Fixes

 * `plot()` for `PairedPSCBS` used a non-defined variable.


# Version 0.7.3 [2010-09-16]

## New Features

 * Added `callABandHighAI()` for calling paired PSCBS segmentation
   results.

## Code Refactoring

 * Added internal bootstrapping functions.


# Version 0.7.2 [2010-09-15]

## New Features

 * Added more methods for the `PSCBS` class.


# Version 0.7.1 [2010-09-08]

## New Features

 * Added more methods for the `PSCBS` class.

 * Now `segmentByPairedPSCBS()` also returns the TumorBoost normalized
   data.


# Version 0.7.0 [2010-09-04]

## New Features

 * Updated `segmentByPairedPSCBS()` to provide two-step segmentation
   from first segmenting the total copy numbers and then the
   decrease-of-heterozygosity signals.  Added utility functions for
   plotting the results.  The code for calling allelic imbalance and
   LOH is still to be added.


# Version 0.6.3 [2010-09-02]

## New Features

 * Now `segmentByCBS()` also works if there are no data points.


# Version 0.6.2 [2010-07-14]

## New Features

 * Added `callNaiveHeterzygotes()`, which is a cleaned up version of
   `findheterozygous()`. Added Rd example that asserts that the two
   are identical and compares the calls to those of
   `aroma.light::callNaiveGenotypes()`.


# Version 0.6.1 [2010-07-09]

## New Features

 * Added low-level `segmentByPairedPSCBS()`, which runs paired PSCBS
   segmentation on a single sample and a single chromosome. It only
   segments; it does not call segments.  This is only a stub in the
   sense that it still does not adjust p-values etc.

 * Added low-level `segmentByCBS()`, which runs CBS segmentation on a
   single sample and a single chromosome.

 * BACKWARD COMPATIBILITY: Now `psCNA()` returns a list of length 8.

 * Reverted `psSegment()` back to v0.5.6.


# Version 0.6.0 [2010-07-08]

## New Features

 * Now `psSegmentPaired()` returns a data frame (no longer a matrix).

 * CLEANUP: Created `psSegmentPaired()` from `psSegment()`.

## Code Refactoring

 * CLEANUP: Major cleanup, i.e. renaming variables, reordering etc.

 * ROBUSTNESS: Replaced all `1:n` with `seq(length = n)` to deal with
   `n == 0`.

 * ROBUSTNESS: Now all list elements are referenced by name.

 * ROBUSTNESS: Now all iterator variables are written as `ii`, `jj`,
   etc.

 * Using `setMethodS3()` of **R.methodsS3** to define S3 methods.

 * Dropping NAMESPACE while package is finished.  This makes it easier
   to patch methods, etc.


# Version 0.5.6 [2010-07-07]

## Documentation

 * Added `example(psSegment)`.

## Bug Fixes

 * Previous clean up introduced bugs.

 * The dynamic library for `hrmode()` was not loaded.


# Version 0.5.5 [2010-05-05]

## Code Refactoring

 * CLEANUP/ROBUSTNESS: Major code clean up.


# Version 0.5.4 [2010-04-30]

## Code Refactoring

 * Added internal `hrmode()`.

 * CLEANUP: Renamed source files to match function names. Only only
   function per source file.


# Version 0.5.3 [2010-04-22]

## Significant Changes

 * ABO updated the psCBS algorithm.


# Version 0.5.2 [2010-0?-??]

 * ???


# Version 0.5.1 [2010-03-31]

## Significant Changes

 * Now `psSegment(..., matching.reference = TRUE)` does TumorBoost
   normalization on the allele B fractions before segmentation.


# Version 0.5.0 [2010-03-12]

 * Added to R-forge repository.
