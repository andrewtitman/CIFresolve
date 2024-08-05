\name{cif_elgalay}
\alias{cif_elgalay}
\docType{data}
\title{
Extracted points from cumulative incidence function curves in El Galay et al (2015)
}
\description{
Manually extracted values of the estimated cumulative incidence function curves of treatment or FL related death and death from unrelated cause before progression, and the corresponding marked censoring times (tick points). To be used as an example input for the CIF_resolve function.
}
\usage{data("cif_elgalay")}
\format{
  A list with the following components
\describe{
   \item{times}{Time in months of extracted increment points}
   \item{cif}{Matrix of the corresponding values of the CIFs; 'Treatment or FL related death' and 'Death from unrelated cause before progression' in columns one and two, respectively.}
   \item{ticks}{Extracted times of the marked censoring times from the plot.}
   \item{nrisk}{Vector of periodically reported numbers at risk}
   \item{trisk}{Times at which number of patients at risk is reported}
 }
}
\keyword{datasets}
