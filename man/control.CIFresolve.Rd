\name{control.CIFresolve}
\alias{control.CIFresolve}
\title{
Control parameters for CIFresolve fitting functions.
}
\description{
Function for defining the control parameters for the CIF_resolve and KM_resolve functions.
}
\usage{
control.CIFresolve(strict_tick=TRUE,strict_dec=TRUE,cen_penalty=1e-3,
constr_tol=1e-8,nprobe=0,epagap=1e-6,epgap=1e-4,tilim=1000,trace=1,
ttol=0.001,cen_max=NULL,ceventinc=FALSE,timeunit=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{strict_tick}{
Logical to indicate whether censoring should be constrained to only occur at the tick marks. Default is FALSE
}
\item{strict_dec}{Logical to indicate whether all identified decrements (of S(t)) should be associated with at least 1 event. Default is TRUE}
\item{cen_penalty}{Penalty term for square of the number of censoring events to avoid non-convex optimization problems. Should be a small positive value. Larger values lead to fewer ties in censoring times. Defaults to 1e-3 }
\item{constr_tol}{Numerical value for tolerance value on constraints in `approx` method to avoid solve.QP erroneously indicating no solutions. Defaults to 1e-8. You can increase this value e.g. to 1e-3 or 1e-2 if this error occurs. Switching to opt_method='miqp' may also resolve the issue.}
\item{nprobe}{Control parameter for miqp method. Defaults to 0. See IBM CPLEX manual for more details}
\item{epagap}{Control parameter for miqp method. Defaults to 1e-6. See IBM CPLEX manual for more details}
\item{epgap}{Control parameter for miqp method. Defaults to 1e-4. See IBM CPLEX manual for more details}
\item{tilim}{Time limit for miqp method optimization (in seconds). Defaults to 1000.}
\item{trace}{Trace control parameter for miqp method. 1=print progress, 0=do not print any progress of optimization.}
\item{ttol}{Tolerance by which the totaltime value is to be recovered. Optimization will constrain total time to (1-ttol) < TT < (1+ttol).}
\item{cen_max}{Maximum number of censoring times at any particular point in the solution (another way of preventing solutions with over clumpy censoring distributions than varying cen_penalty).}
\item{ceventinc}{Logical to indicate whether the cumulative event counts (if supplied) include events occurring exactly at the time (TRUE), or only include events occurring strictly before time t (FALSE). Defaults to FALSE.}
\item{timeunit}{Value to add to maximum observed time if n.risk does not include a zero value. Defaults to 1.}
}
\details{
The resulting list can be used as the control argument in \link{KM_resolve} or \link{CIF_resolve}.
}
\value{
Returns a list with the specified control arguments.
}

\author{
Andrew Titman \email{a.titman@lancaster.ac.uk}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CIF_resolve}}, \code{\link{plot.KMresolve}}}
