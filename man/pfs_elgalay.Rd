\name{pfs_elgalay}
\alias{pfs_elgalay}
\docType{data}
\title{
Extracted points from progression-free survival curve in El Galay et al (2015).
}
\description{
Manually extracted values of the Kaplan-Meier estimate of progression-free survival and the corresponding marked censoring times (tick points). To be used as an example input for the KM_resolve function.
}
\usage{data("pfs_elgalay")}
\format{
  A list with the following components
\describe{
   \item{times}{Time in months of extracted decrement points}
   \item{surv}{Corresponding Kaplan-Meier estimate values}
   \item{ticks}{Extracted times of the marked censoring times from the plot.}
   \item{nrisk}{Vector of periodically reported numbers at risk}
   \item{trisk}{Times at which number of patients at risk is reported}
 }
}
\keyword{datasets}
