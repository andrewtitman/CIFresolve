\name{plot.KMresolve}
\alias{plot.KMresolve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot results from a reconstruction of individual patient data from a Kaplan-Meier curve.
}
\description{
Plotting method for KMresolve objects made using \code{KM_resolve}.
}
\usage{
\method{plot}{KMresolve}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  A KMresolve model object generated from \code{KM_resolve}.
}

\item{...}{further arguments passed to or from other methods.}
}
\value{
Produces a two or three panel plot comparing the original and generated Kaplan-Meier curve, and also comparing the estimated patients at risk with those supplied. If supplied, the third panel compares the estimated cumulative events with those given.
}

\author{
Andrew Titman \email{a.titman@lancaster.ac.uk}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{KM_resolve}}}
