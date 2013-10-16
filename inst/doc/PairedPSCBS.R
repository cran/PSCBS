###########################################################################
## This 'tangle' R script was created from an RSP document.
## RSP source document: 'PairedPSCBS.tex.rsp'
## Metadata 'title': 'Parent-specific copy-number segmentation using Paired PSCBS'
## Metadata 'author': 'Henrik Bengtsson'
## Metadata 'keywords': 'copy numbers, allele specific, parent specific, genomic aberrations'
## Metadata 'engine': 'R.rsp::rsp'
###########################################################################

t0 <- Sys.time()
library("PSCBS");
library("R.devices");
evalCapture <- R.utils::evalCapture;
PSCBS <- R.oo::Package("PSCBS");
R.rsp <- R.oo::Package("R.rsp");
fixLocations <- function(fit, ...) {
  for (key in grep("(End|Start)$", colnames(fit$output))) {
    fit$output[[key]] <- as.integer(fit$output[[key]]);
  }
  fit;
} # fixLocations()

devOptions("png", width=840);
options(width=85);
options(digits=3);
options(str=strOptions(strict.width="cut"));
R.rsp$version
R.rsp$author
format(as.Date(PSCBS$date), format="%B %d, %Y")
fullname <- "PairedPSCBS,exData,chr01";
evalCapture({
data <- PSCBS::exampleData("paired.chr01")
str(data)
})
evalCapture({
data <- dropSegmentationOutliers(data)
})
evalCapture({
gaps <- findLargeGaps(data, minLength=1e6)
gaps
})
evalCapture({
knownSegments <- gapsToSegments(gaps)
knownSegments
})
evalCapture({
fit <- segmentByPairedPSCBS(data, knownSegments=knownSegments, seed=0xBEEF, verbose=-10)
})
nbrOfSegments(fit)
fit <- fixLocations(fit);
evalCapture({
getSegments(fit, simplify=TRUE)
})
segs <- getSegments(fit, simplify=TRUE)
which(segs$tcnNbrOfLoci == 0)
which(segs$dhNbrOfLoci == 0)
toPNG(fullname, tags=c(class(fit)[1L], "tracks"), aspectRatio=0.6, {
    plotTracks(fit);
  })
evalCapture({
fit <- callROH(fit, verbose=-10)
})
evalCapture({
fit <- callAB(fit, verbose=-10)
})
evalCapture({
fit <- callLOH(fit, verbose=-10)
})
evalCapture({
fit <- callNTCN(fit, verbose=-10)
})
evalCapture({
getSegments(fit, simplify=TRUE)
})
segs <- getSegments(fit, simplify=TRUE)
evalCapture({
fitP <- pruneByHClust(fit, h=0.25, verbose=-10)
})
toPNG(fullname, tags=c(class(fitP)[1L], "pruned", "tracks"), aspectRatio=0.6, {
    plotTracks(fitP);
  })
toLatex(sessionInfo())
dt <- round(Sys.time()-t0, digits=2)
attr(dt, "units")
