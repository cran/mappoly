% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_genome_vs_map.R
\name{plot_genome_vs_map}
\alias{plot_genome_vs_map}
\title{Physical versus genetic distance}
\usage{
plot_genome_vs_map(
  map.list,
  phase.config = "best",
  same.ch.lg = FALSE,
  alpha = 1/5,
  size = 3
)
}
\arguments{
\item{map.list}{A list or a single object of class \code{mappoly.map}}

\item{phase.config}{A vector containing which phase configuration should be
plotted. If \code{'best'} (default), plots the configuration
with the highest likelihood for all elements in \code{'map.list'}}

\item{same.ch.lg}{Logical. If \code{TRUE} displays only the scatterplots between the 
chromosomes and linkage groups with the same number. Default is \code{FALSE}.}

\item{alpha}{transparency factor for SNPs points}

\item{size}{size of the SNP points}
}
\description{
This function plots scatterplot(s) of physical distance (in Mbp) versus the genetic 
distance (in cM). Map(s) should be passed as a single object or a list of objects 
of class \code{mappoly.map}.
}
\examples{
  plot_genome_vs_map(solcap.mds.map, same.ch.lg = TRUE)
  plot_genome_vs_map(solcap.mds.map, same.ch.lg = FALSE, 
                     alpha = 1, size = 1/2)
 
}
\references{
Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
    analysis and haplotype phasing in experimental autopolyploid
    populations with high ploidy level using hidden Markov
    models, _G3: Genes, Genomes, Genetics_. 
    \doi{10.1534/g3.119.400378}
}
\author{
Marcelo Mollinari, \email{mmollin@ncsu.edu}
}
