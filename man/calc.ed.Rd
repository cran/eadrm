% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ea_funcs.R
\name{calc.ed}
\alias{calc.ed}
\title{Finds the dose that corresponds to a particular level of the response}
\usage{
calc.ed(eadrm.obj, response, ci.obj = NULL)
}
\arguments{
\item{eadrm.obj}{Fitted eadrm model object.}

\item{response}{The desired response level.}

\item{ci.obj}{Output from eadrm.ci that is used to compute
confidence intervals. Defaults to NULL, in which case no confidence
intervals are computed.}
}
\value{
If no confidence intervals are requested, it returns the
concentration corresponding to the specified response.
Otherwise returns a list of three values ed.med, ed.l95, and
ed.u95, which correspond to the median concentration and the
corresponding lower/upper 95\% confidence bounds.
}
\description{
Calculates the concentration that induces a particular level of the
response.
}
\examples{
ea.fit <- eadrm(CarboA$y, CarboA$x)
calc.ec(ea.fit)
}
\seealso{
\code{\link{eadrm}}, \code{\link{eadrm.ci}}
}
