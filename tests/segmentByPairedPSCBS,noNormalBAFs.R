library("PSCBS")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# (note to package developers: this example data set may
#  be replaced in a future release of the package)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname <- system.file("data-ex/PairedPSCBS,exData,chr01.Rbin", package="PSCBS")
data <- R.utils::loadObject(pathname)

# Order by chromosome and position
o <- order(data$chromosome, data$position)
data <- data[o,]
str(data)
R.oo::attachLocally(data)
x <- position
J <- length(x)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate that genotypes are known by other means
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library("aroma.light")
muN <- aroma.light::callNaiveGenotypes(betaN, censorAt=c(0,1))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop single-locus outliers
CTs <- dropSegmentationOutliers(CT, chromosome=1, x=x, verbose=-10)

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(CTs, betaT=betaT, muN=muN, tbn=FALSE,
                            chromosome=1, x=x, 
                            seed=0xBEEF, verbose=-10)
print(fit)

# Plot results
plotTracks(fit)

# Sanity check
stopifnot(nbrOfSegments(fit) == 10)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Bootstrap segment level estimates
# (used by the AB caller, which, if skipped here,
#  will do it automatically)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- bootstrapTCNandDHByRegion(fit, verbose=-10)
print(fit)
plotTracks(fit)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calling segments in allelic balance (AB) and
# in loss-of-heterozygosity (LOH)
# NOTE: Ideally, this should be done on whole-genome data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- callAB(fit, verbose=-10)
fit <- callLOH(fit, verbose=-10)
print(fit)
plotTracks(fit)