% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{hexafake.geno.dist}
\alias{hexafake.geno.dist}
\title{Simulated autohexaploid dataset with genotype probabilities.}
\format{
An object of class \code{mappoly.data} which contains a
    list with the following components:
\describe{
    \item{ploidy}{ploidy level = 6}
    \item{n.ind}{number individuals = 300}
    \item{n.mrk}{total number of markers = 1500}
    \item{ind.names}{the names of the individuals}
    \item{mrk.names}{the names of the markers}
    \item{dosage.p1}{a vector containing the dosage in
      parent P for all \code{n.mrk} markers}
    \item{dosage.p2}{a vector containing the dosage in
      parent Q for all \code{n.mrk} markers}
    \item{chrom}{a vector indicating which sequence each marker
      belongs. Zero indicates that the marker was not assigned to any
      sequence}
    \item{genome.pos}{Physical position of the markers into the
      sequence}
    \item{prob.thres = 0.95}{probability threshold to associate a marker 
                             call to a dosage. Markers with maximum genotype 
                             probability smaller than 'prob.thres' are considered 
                             as missing data for the dosage calling purposes}
    \item{geno}{a data.frame 
      containing the probability distribution for each combination of
      marker and offspring. The first two columns represent the marker
      and the offspring, respectively. The remaining elements represent
      the probability associated to each one of the possible
      dosages}
      \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
      for each individual (columns). Missing data are represented by 
      \code{ploidy_level + 1 = 7}}
    \item{n.phen}{There are no phenotypes in this simulation}
    \item{phen}{There are no phenotypes in this simulation}
}
}
\usage{
hexafake.geno.dist
}
\description{
A dataset of a hypothetical autohexaploid full-sib population 
containing three homology groups. This dataset contains the
probability distribution of the genotypes and 2\% of missing data, 
but is essentially the same dataset found in \code{\link[mappoly]{hexafake}}
}
\keyword{datasets}
