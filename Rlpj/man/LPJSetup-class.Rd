% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LPJSetupClass.R
\docType{class}
\name{LPJSetup-class}
\alias{LPJSetup-class}
\title{A LPJSetup class object}
\value{
a S4 class object with five slots
\itemize{
\item  clusterType a character string indicating the type of cluster to be created (MPI o SOCK)
\item  numCores a integer specifying number of cores of the cluster
\item  mainDir  a character string indicating the path to the directory where
 the model link and the templates are located
 the directory structure for the outputs
\item  runDir a character string indicating the path to the directory where
 the needed files for model run are located
\item  outDir a character string indicating the path to the directory where
 the model outputs will be saved
}
}
\description{
This is a S4 class object to store the configuration for parallel
runs to be used by runLPJ.
}
\author{
Ramiro Silveyra Gonzalez
}
\keyword{Rlpj}

