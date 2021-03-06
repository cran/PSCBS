%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CBS.IO.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{writeSegments.CBS}
\alias{writeSegments.CBS}
\alias{CBS.writeSegments}
\alias{writeSegments,CBS-method}
\alias{writeWIG}
\alias{writeWIG.AbstractCBS}

\title{Writes the table of segments to file}

\description{
 Writes the table of segments to file.
}

\usage{
\method{writeSegments}{CBS}(fit, name=getSampleName(fit), tags=NULL, ext="tsv", path=NULL,
  addHeader=TRUE, createdBy=NULL, sep="\t", nbrOfDecimals=4L, splitters=FALSE,
  overwrite=FALSE, skip=FALSE, ...)
}

\arguments{
  \item{name, tags}{Name and optional tags part of the filename}.
  \item{path}{The directory where the file will be written.}
  \item{addHeader}{If \code{\link[base:logical]{TRUE}}, header comments are written.}
  \item{createdBy}{A header comment of whom created the file.}
  \item{splitters}{If \code{\link[base:logical]{TRUE}}, each chromosome is separated by a row
    of missing values.}
  \item{overwrite, skip}{If an output file already exists, these
    arguments specifies what should happen.}
  \item{...}{Additional arguments pass to \code{getSegments()}.}
}

\value{
  Returns the pathname of the the file written.
}

\author{Henrik Bengtsson}

\seealso{
  Utilizes \code{\link[PSCBS:getSegments.CBS]{*getSegments}()}.
  For more information see \code{\link{CBS}}..
}


\keyword{internal}
\keyword{methods}
