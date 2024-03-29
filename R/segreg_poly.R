#' Polysomic segregation frequency
#'
#' Computes the polysomic segregation frequency given a ploidy level
#' and the dosage of the locus in both parents. It does not consider
#' double reduction.
#'
#' @param ploidy the ploidy level
#'
#' @param dP the dosage in parent P
#'
#' @param dQ the dosage in parent Q
#'
#' @return a vector containing the expected segregation frequency for
#'     all possible genotypic classes.
#'
#' @examples
#' # autohexaploid with two and three doses in parents P and Q,
#' # respectively
#' seg <- segreg_poly(ploidy = 6, dP = 2, dQ = 3)
#' barplot(seg, las = 2)
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400378}
#'     
#'     Serang O, Mollinari M, Garcia AAF (2012) Efficient Exact 
#'     Maximum a Posteriori Computation for Bayesian SNP 
#'     Genotyping in Polyploids. _PLoS ONE_ 7(2): 
#'     e30906.
#'     
#'
#' @importFrom stats dhyper
#' @export segreg_poly
#'
segreg_poly <- function(ploidy, dP, dQ) {
    if (ploidy%%2 != 0)
        stop("m must be an even number")
    p.dose <- numeric((ploidy + 1))
    p.names <- character((ploidy + 1))
    seg.p1 <- dhyper(x = c(0:(ploidy + 1)), m = dP, n = (ploidy - dP), k = ploidy/2)
    seg.p2 <- dhyper(x = c(0:(ploidy + 1)), m = dQ, n = (ploidy - dQ), k = ploidy/2)
    M <- tcrossprod(seg.p1, seg.p2)
    for (i in 1:nrow(M)) {
        for (j in 1:ncol(M)) {
            p.dose[i + j - 1] <- p.dose[i + j - 1] + M[i, j]
        }
    }
    p.dose <- p.dose[!is.na(p.dose)]
    for (i in 0:ploidy) p.names[i + 1] <- paste(paste(rep("A", i), collapse = ""), paste(rep("a", (ploidy - i)), collapse = ""), sep = "")
    names(p.dose) <- p.names
    return(p.dose)
}
