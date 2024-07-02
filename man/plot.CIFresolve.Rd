\name{plot.CIFresolve}
\alias{plot.CIFresolve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot results from a reconstruction of individual patient data from a set of cumulative incidence functions.
}
\description{
Plotting method for CIFresolve objects made using \code{CIF_resolve}.
}
\usage{
\method{plot}{CIFresolve}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A CIFresolve model object generated from \code{CIF_resolve}.
}
\item{...}{further arguments passed to or from other methods.}
}
\value{
Produces a multi-panel plot comparing the original and generated cumulative incidence function curves, and also comparing the estimated patients at risk with those supplied.
}

\author{
Andrew Titman \email{a.titman@lancaster.ac.uk}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CIF_resolve}}}
