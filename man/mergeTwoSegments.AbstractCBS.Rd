%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.RESTRUCT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{mergeTwoSegments.AbstractCBS}
\alias{mergeTwoSegments.AbstractCBS}
\alias{AbstractCBS.mergeTwoSegments}
\alias{mergeTwoSegments,AbstractCBS-method}
\alias{AbstractCBS.dropChangePoint}
\alias{dropChangePoint.AbstractCBS}
\alias{dropChangePoint,AbstractCBS-method}


\title{Merge two neighboring segments}

\description{
  Merge two neighboring segments into one segment, which is done by dropping their
  common change point and recalculating the segment statistics.
}

\usage{
\method{mergeTwoSegments}{AbstractCBS}(...)
}

\arguments{
 \item{...}{Not used.}
}

\value{
  Returns an \code{\link{AbstractCBS}} of the same class with one less segment.
}

\author{Henrik Bengtsson}

\seealso{
  To merge a segment and its two flanking segments, see
  \code{\link[PSCBS:mergeThreeSegments.AbstractCBS]{*mergeThreeSegments}()}.
  To drop regions (a connected set of segments)
  see \code{\link[PSCBS:dropRegions.AbstractCBS]{*dropRegions}()}.
  For more information see \code{\link{AbstractCBS}}.
}
\keyword{internal}
\keyword{methods}
