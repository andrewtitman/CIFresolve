\name{align_CIF}
\alias{align_CIF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Align a list of separate cumulative incidence functions into a single matrix.
}
\description{
Utility function for converting a set of separate digitized cumulative incidence functions into a single object, for use in \link{CIF_resolve}.
}
\usage{
align_CIF(..., res_digits=4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{
 One or more lists containing named vectors 'time' and 'cif' representing the times and corresponding cumulative incidence function values for one competing risk.
}
\item{res_digits}{
Number of digits to round the times. This defines which increment points from different CIFs are assumed to be at the same time. Defaults to 4 corresponding to rounding to nearest 0.0001.}
}
\details{
The \link{CIF_resolve} function requires the cumulative incidence functions to be supplied on a common set of times. This function allows independently digitized CIFs to be combined into the appropriate matrix.
}
\value{
Returns a list containing the common set of times and the matrix of CIF values at those times.
}

\author{
Andrew Titman \email{a.titman@lancaster.ac.uk}
}

\seealso{
\code{\link{CIF_resolve}}}
