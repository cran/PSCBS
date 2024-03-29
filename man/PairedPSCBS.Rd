%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{PairedPSCBS}
\docType{class}
\alias{PairedPSCBS}


\title{The PairedPSCBS class}

\description{
 Package:  PSCBS \cr
\bold{Class PairedPSCBS}\cr

\code{list}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[PSCBS]{AbstractCBS}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[PSCBS]{PSCBS}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{PairedPSCBS}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract class \bold{PairedPSCBS}\cr
extends \emph{\link[PSCBS]{PSCBS}}\cr



 A PairedPSCBS is an object containing the results from the
 Paired PSCBS method.
}

\usage{PairedPSCBS(fit=list(), ...)}

\arguments{
  \item{fit}{A \code{\link[base]{list}} structure containing the Paired PSCBS results.}
  \item{...}{Not used.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{callAB} \tab  -\cr
 \tab \code{callCopyNeutral} \tab  -\cr
 \tab \code{callGNL} \tab  -\cr
 \tab \code{callGNLByTCNofAB} \tab  -\cr
 \tab \code{callGainNeutralLoss} \tab  -\cr
 \tab \code{callLOH} \tab  -\cr
 \tab \code{callNTCN} \tab  -\cr
 \tab \code{callROH} \tab  -\cr
 \tab \code{estimateDeltaAB} \tab  -\cr
 \tab \code{estimateDeltaLOH} \tab  -\cr
 \tab \code{estimateKappa} \tab  -\cr
 \tab \code{extractCNs} \tab  -\cr
 \tab \code{hasBootstrapSummaries} \tab  -\cr
 \tab \code{plotTracks} \tab  -\cr
 \tab \code{segmentByNonPairedPSCBS} \tab  -\cr
 \tab \code{segmentByPairedPSCBS} \tab  -\cr
 \tab \code{seqOfSegmentsByDP} \tab  -\cr
}


 \bold{Methods inherited from PSCBS}:\cr
as.data.frame, c, drawChangePoints, extractChromosomes, extractWIG, getLocusData, getLocusSignalNames, getSegmentTrackPrefixes, isLocallyPhased, isSegmentSplitter, normalizeTotalCNs, writeSegments

 \bold{Methods inherited from AbstractCBS}:\cr
adjustPloidyScale, all.equal, as.data.frame, clearCalls, drawChangePoints, drawKnownSegments, dropChangePoint, dropChangePoints, dropRegion, dropRegions, extractCNs, extractChromosome, extractChromosomes, extractRegions, extractSegments, extractWIG, getChangePoints, getChromosomeOffsets, getChromosomeRanges, getChromosomes, getLocusData, getLocusSignalNames, getMeanEstimators, getSampleName, getSegmentSizes, getSegmentTrackPrefixes, getSegments, mergeThreeSegments, mergeTwoSegments, nbrOfChangePoints, nbrOfChromosomes, nbrOfLoci, nbrOfSegments, normalizeTotalCNs, ploidy, ploidy<-, plotTracks, print, pruneByDP, pruneByHClust, renameChromosomes, report, resegment, resetSegments, sampleCNs, sampleName, sampleName<-, seqOfSegmentsByDP, setLocusData, setMeanEstimators, setPloidy, setSampleName, setSegments, shiftTCN, tileChromosomes, updateMeans, writeWIG

 \bold{Methods inherited from list}:\cr
Ops,nonStructure,vector-method, Ops,structure,vector-method, Ops,vector,nonStructure-method, Ops,vector,structure-method, all.equal, as.data.frame, attachLocally, callHooks, coerce,ANY,list-method, relist, type.convert, within


}

\author{Henrik Bengtsson}

\seealso{
  The \code{\link{segmentByPairedPSCBS}}() method returns an object of this class.
}
\keyword{classes}
