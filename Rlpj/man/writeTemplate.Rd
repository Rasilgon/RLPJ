% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/writeTemplate.R
\name{writeTemplate}
\alias{writeTemplate}
\title{A writing template function for LPJ}
\usage{
writeTemplate(template1 = NULL, parameterList = NULL, runDir = NULL)
}
\arguments{
\item{template1}{a character string providing the general model template,
e.g, global.ins. Provide only the file name, not the path}

\item{parameterList}{a named list containing the parameters to be calibrated}

\item{runDir}{a character string indicating path to the run directory}
}
\value{
none
}
\description{
This function reads a template, and replaces parameters with the
provides parameters list.If any parameters values is not provided, the function will
set it to the default values. The function assumes that a copy of the template
 is already placed in the run folder.
}
\note{
Based an older code of Istem Fer, Uni Potsdam
}
\examples{
\dontrun{
writeTemplate("global.ins", list(run_lamda_max = 0.5, run_emax= 5),
              "/home/lpjRun/runDirectory1")
}
}
\author{
Florian Hartig, Ramiro Silveyra Gonzalez, Maurizio Bagnara
}
\keyword{Rlpj}

