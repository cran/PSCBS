###########################################################################
## This 'tangle' R script was created from an RSP document.
## RSP source document: 'CBS.tex.rsp'
## Metadata 'title': 'Total copy-number segmentation using CBS'
## Metadata 'keywords': 'copy numbers, genomic aberrations'
## Metadata 'author': 'Henrik Bengtsson'
## Metadata 'engine': 'R.rsp::rsp'
###########################################################################

t0 <- Sys.time()
library("PSCBS");
library("R.devices");
evalCapture <- R.utils::evalCapture;
PSCBS <- R.oo::Package("PSCBS");
R.rsp <- R.oo::Package("R.rsp");
fixLocations <- function(fit, ...) {
  for (key in grep("(end|start)$", colnames(fit$output))) {
    fit$output[[key]] <- as.integer(fit$output[[key]]);
  }
  fit;
} # fixLocations()

signalType <- "TCN";

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
data <- data[,c("chromosome", "x", "CT")]
colnames(data)[3] <- "y"
str(data)
})
signalType
signalType
evalCapture({
data <- dropSegmentationOutliers(data)
})
signalType
evalCapture({
gaps <- findLargeGaps(data, minLength=1e6)
gaps
})
evalCapture({
knownSegments <- gapsToSegments(gaps)
knownSegments
})
signalType
signalType
evalCapture({
fit <- segmentByCBS(data, knownSegments=knownSegments, seed=0xBEEF, verbose=-10)
})
signalType
nbrOfSegments(fit)
signalType
fit <- fixLocations(fit);
evalCapture({
getSegments(fit, simplify=TRUE)
})
segs <- getSegments(fit, simplify=TRUE)
which(segs$nbrOfLoci == 0)
signalType
signalType
signalType
toPNG(fullname, tags=c("tracks"), aspectRatio=0.35, {
    plotTracks(fit);
  })
signalType
signalType
signalType
evalCapture({
getSegments(fit, simplify=TRUE)
})
segs <- getSegments(fit, simplify=TRUE)
signalType
signalType
evalCapture({
fitP <- pruneByHClust(fit, h=0.25, verbose=-10)
})
toPNG(fullname, tags=c("pruned", "tracks"), aspectRatio=0.35, {
    plotTracks(fitP);
  })
signalType
toLatex(sessionInfo())
dt <- round(Sys.time()-t0, digits=2)
attr(dt, "units")
