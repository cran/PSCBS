###########################################################################
### This 'tangle' R script was created from an RSP document.
### RSP source document: 'CBS.tex.rsp'
### Metadata 'title': 'Total copy-number segmentation using CBS'
### Metadata 'author': 'Henrik Bengtsson'
### Metadata 'engine': 'R.rsp::rsp'
### Metadata 'keywords': 'copy numbers, genomic aberrations'
###########################################################################

t0 <- Sys.time()
R.utils::use("R.utils")

# RSP specific
R.rsp <- R.oo::Package("R.rsp")
withCapture <- R.utils::withCapture
options("withCapture/newline"=FALSE)
options(str=strOptions(strict.width="cut"))
options(width=85)
options(digits=3)

# Graphics
use("R.devices")
options("devEval/args/field"="fullname") # Preferred for LaTeX
devOptions("png", width=840)

# Analysis
use("PSCBS")
PSCBS <- R.oo::Package("PSCBS")
fixLocations <- function(fit, ...) {
  for (key in grep("(end|start)$", colnames(fit$output))) {
    fit$output[[key]] <- as.integer(fit$output[[key]])
  }
  fit
} # fixLocations()

signalType <- "TCN"
R.rsp$version
R.rsp$author
format(as.Date(PSCBS$date), format="%B %d, %Y")
fullname <- "PairedPSCBS,exData,chr01"
withCapture({
data <- PSCBS::exampleData("paired.chr01")
data <- data[,c("chromosome", "x", "CT")]
colnames(data)[3] <- "y"
str(data)
})
signalType
signalType
withCapture({
data <- dropSegmentationOutliers(data)
})
signalType
withCapture({
gaps <- findLargeGaps(data, minLength=1e6)
gaps
})
withCapture({
knownSegments <- gapsToSegments(gaps)
knownSegments
})
signalType
signalType
withCapture({
fit <- segmentByCBS(data, knownSegments=knownSegments, seed=0xBEEF, verbose=-10)
})
signalType
nbrOfSegments(fit)
signalType
fit <- fixLocations(fit)
withCapture({
getSegments(fit, simplify=TRUE)
})
segs <- getSegments(fit, simplify=TRUE)
which(segs$nbrOfLoci == 0)
signalType
signalType
signalType
toPNG(fullname, tags=c(class(fit)[1L], "tracks"), aspectRatio=0.35, {
    plotTracks(fit)
  })
signalType
signalType
signalType
signalType
withCapture({
fitP <- pruneByHClust(fit, h=0.25, verbose=-10)
})
toPNG(fullname, tags=c(class(fitP)[1L], "pruned", "tracks"), aspectRatio=0.35, {
    plotTracks(fitP)
  })
signalType
toLatex(sessionInfo())
dt <- round(Sys.time()-t0, digits=2)
attr(dt, "units")
