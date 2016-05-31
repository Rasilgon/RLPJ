# Rlpj

A R package that wraps the LPJ-GUESS model

#### Installing Rlpj

Rlpj depends on three packages:

- zoo
- snow (for SOCK cluster)
- Rmpi (for MPI clusters)

Installing Rmpi might be complicated. If you are thinking of using Rlpj on a laptop or workstation, you will be dealing with SOCK clusters and you do not need Rmpi.


#### Recommended installation

Install the latest stable [release](link) with the devtools package as follows  

```{r}
library(devtools)
install_url("xxxxx")
library(Rlpj)
?Rlpj
```
You may want to check if the link above is really the latest stable realease. 

#### Build the package yourself 

Clone the package to your computer, and build with hand or Rstudio. If you need help see here http://biometry.github.io/APES/R/R70-PackageDevelopment.html


In Rstudio, the vignette may not be built per default. You will turn this on in your project options, or run 

```{r}
devtools::build_vignettes("BayesianTools")
```
