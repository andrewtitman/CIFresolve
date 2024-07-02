\name{CIF_resolve}
\alias{CIF_resolve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reconstruct patient-level data from cumulative incidence function estimates.
}
\description{
Outputs estimates of the number of competing risks events and number of censored observation from supplied digitized data points.
}
\usage{
CIF_resolve(S, t.risk, n.risk, nevent= NULL, ndeath = NULL, ticks = NULL,
strict_tick = FALSE, strict_dec=TRUE, totaltime=NULL,
totaltime_power=1, c.event=NULL, t.event=NULL, optmethod = "approx",
cen_penalty=1e-3, constr_tol = 1e-8, nprobe = 0,
epagap = 1e-06, epgap = 1e-04, tilim = 1000, trace=1,
ttol=0.01, cen_max=NULL, ceventinc=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
 List containing a vector `time` corresponding to the times at which the cumulative incidence curve has been obtained, and a matrix `CIF` giving the corresponding cumulative incidence function values where the number of columns should correspond to the number of risks.
}
  \item{t.risk}{
 Vector containing the timepoints at which the number of patients at risk/under observation is known.
}
\item{n.risk}{Vector containing the number of patients at risk/under observation at the times given in the vector `t.risk`}
\item{nevent}{Vector of number of events of each type (should be of same length as number of risks/number of columns in S$cif). Can be omitted. See also c.event if interim counts are available.}
\item{ndeath}{Vector of total number of events across all types. Can be omitted.}
\item{ticks}{Optional vector of time points where the CIF plot had tick marks representing censoring times}
\item{strict_tick}{Logical to indicate whether censoring should be constrained to only occur at the tick marks. Default is FALSE}
\item{strict_dec}{Logical to indicate whether all identified decrements (of S(t)) or increments (of F_j(t)) should be associated with at least 1 event. Default is TRUE}
\item{totaltime}{Total time at risk in the dataset (sum of all patient's follow-up times) or sum of transformed follow-up times (sum t_{i}^pow). Can be omitted.}
\item{totaltime_power}{The power by which each of the individual times has been raised in the calculation of totaltime, defaults to 1, potentially useful if the MLEs of the Weibull rate and shape parameters are known.}
\item{c.event}{Vector or matrix of cumulative number of events; should either be a vector (in which case is interpreted as the cumulative total number of events) or a matrix where each column represents the cumulative number of events of each time. Can be omitted}
\item{t.event}{Vector of times at which number of cumulative events is supplied. If omitted it is assumed t.event=t.risk}
\item{optmethod}{Character string to indicate method of finding solution: `approx` - solve the continuous QP and then perform integer rounding. `miqp` - use full mixed integer QP optimization using Rcplex. Default is `approx`}
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
}
\details{
The \code{CIF_resolve} function is used to reconstruct pseudo patient-level data from a set of digitized cumulative incidence function curves (estimated via the Aalen-Johansen estimator).
}
\value{
Returns a data frame containing the table of estimated numbers at risk, numbers of events and number of censorings at time points.
}
\author{
Andrew Titman \email{a.titman@lancaster.ac.uk}
}

\seealso{
\code{\link{KM_resolve}}, \code{\link{plot.CIFresolve}}, \code{\link{make_data}}}
