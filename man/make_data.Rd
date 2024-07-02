\name{make_data}
\alias{make_data}
\title{
Generate a dataset with individual patient data from a reconstructed dataset.
}
\description{
Generates a dataset with individual patient data, that can be directly used in e.g. the survival package, from a CIFresolve or KMresolve object.
}
\usage{
make_data(cif, cen_method=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cif}{
An object of class CIFresolve or KMresolve, generated using \code{CIF_resolve} or \code{KM_resolve}.
}
\item{cen_method}{
Method for assigning timing of censoring events: "start" - Censoring time taken at the start of the interval that censoring could occur (e.g. at an event time); "mid" - Censoring time taken at the middle of the potential interval. Can be omitted, in which case the method depends on whether marked times were given in the reconstruction; if ticks were included in the reconstruction default is "start", otherwise is "mid".
}
}
\value{
Returns a data frame of pseudo-IPD giving the times (in column 'time') and event indicators (in column 'event').
}
\author{
Andrew Titman \email{a.titman@lancaster.ac.uk}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{KM_resolve}}, \code{\link{CIF_resolve}}}
