#' Import from updog
#'
#' Read objects with information related to genotype calling in polyploids.
#' Currently this function supports output objects created with the
#' \code{updog} (output of \code{multidog} function) package.
#' This function creates an object of class \code{mappoly.data}
#'
#' @param object the name of the object of class \code{multidog}
#'     
#' @param prob.thres probability threshold to associate a marker call to a 
#'     dosage. Markers with maximum genotype probability smaller than 'prob.thres' 
#'     are considered as missing data for the dosage calling purposes
#'     
#' @param filter.non.conforming if \code{TRUE} (default) exclude samples with non 
#'     expected genotypes under random chromosome pairing and no double reduction
#'     
#' @param chrom a vector indicating which sequence each marker
#'       belongs. Zero indicates that the marker was not assigned to any
#'       sequence
#'       
#' @param genome.pos vector with physical position of the markers into the
#'       sequence
#'
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#' 
#' @return An object of class \code{mappoly.data} which contains a
#'     list with the following components:
#'     \item{ploidy}{ploidy level}
#'     \item{n.ind}{number individuals}
#'     \item{n.mrk}{total number of markers}
#'     \item{ind.names}{the names of the individuals}
#'     \item{mrk.names}{the names of the markers}
#'     \item{dosage.p1}{a vector containing the dosage in
#'       parent P for all \code{n.mrk} markers}
#'     \item{dosage.p2}{a vector containing the dosage in
#'       parent Q for all \code{n.mrk} markers}
#'     \item{chrom}{a vector indicating which sequence each marker
#'       belongs. Zero indicates that the marker was not assigned to any
#'       sequence}
#'     \item{genome.pos}{physical position of the markers into the
#'       sequence}
#'     \item{prob.thres}{probability threshold to associate a marker call to a 
#'     dosage. Markers with maximum genotype probability smaller than 'prob.thres' 
#'     were considered as missing data in the 'geno.dose' matrix}
#'     \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
#'       for each individual (columns). Missing data are represented by 
#'       \code{ploidy_level + 1}}
#'     \item{geno}{a data.frame 
#'       containing the probability distribution for each combination of
#'       marker and offspring. The first two columns represent the marker
#'       and the offspring, respectively. The remaining elements represent
#'       the probability associated to each one of the possible
#'       dosages. Missing data are converted from \code{NA} to the expected
#'       segregation ratio using function \code{\link[mappoly]{segreg_poly}}}
#'     \item{n.phen}{number of phenotypic traits}
#'     \item{phen}{a matrix containing the phenotypic data. The rows
#'                 correspond to the traits and the columns correspond
#'                 to the individuals}
#'     \item{chisq.pval}{a vector containing p-values related to the chi-squared 
#'     test of Mendelian segregation performed for all markers}
#'     
#' @examples
#' \donttest{
#' if(requireNamespace("updog", quietly = TRUE)){
#' library("updog")
#' data("uitdewilligen")
#' mout = multidog(refmat = t(uitdewilligen$refmat), 
#'                 sizemat = t(uitdewilligen$sizemat), 
#'                 ploidy = uitdewilligen$ploidy, 
#'                 model = "f1",
#'                 p1_id = colnames(t(uitdewilligen$sizemat))[1],
#'                 p2_id = colnames(t(uitdewilligen$sizemat))[2],
#'                 nc = 2)
#' mydata = import_from_updog(mout)
#' mydata
#' plot(mydata)
#' }
#' }
#'
#' @author Gabriel Gesteira, \email{gdesiqu@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400378}
#'     
#' @export import_from_updog
#' 
import_from_updog = function(object, prob.thres = 0.95, 
                             filter.non.conforming = TRUE,
                             chrom = NULL,
                             genome.pos = NULL,
                             verbose = TRUE){
  # Case 1: updog
  if (inherits(object, "multidog")){
    ploidy = object$snpdf$ploidy[1]
    dosage.p1 = as.integer(object$snpdf$p1geno)
    names(dosage.p1) <- object$snpdf$snp
    dosage.p2 = as.integer(object$snpdf$p2geno)
    names(dosage.p2) <- object$snpdf$snp
    if (is.null(prob.thres)) prob.thres = 0.95
    # Selecting conforming markers
    dp = abs(abs(dosage.p1-(ploidy/2))-(ploidy/2))
    dq = abs(abs(dosage.p2-(ploidy/2))-(ploidy/2))
    id = dp+dq != 0
    dosage.p1 = dosage.p1[id]
    dosage.p2 = dosage.p2[id]
    mrk.names = unique(as.character(object$snpdf$snp))[id]
    mrk.sel = which(object$inddf$snp %in% mrk.names)
    geno <- object$inddf[mrk.sel,]
    geno <- geno[c(1,2,grep(pattern = "Pr_", colnames(geno)))]
    colnames(geno) <- c("mrk", "ind", 0:ploidy)
    geno <- geno[order(geno$ind, geno$mrk, decreasing = TRUE),]
    ind.names = as.character(unique(geno$ind))
    mrk.names = as.character(unique(geno$mrk))
    n.ind = length(ind.names)
    n.mrk = length(mrk.names)
    if(n.ind * n.mrk != nrow(geno))
      stop("Check your dataset.")
    if(!is.null(chrom) & length(chrom) != length(object$snpdf$snp)) {
      stop("Check 'chrom' input. The vector should have length equal to the number of markers in the updog output.")
    } else if(!is.null(chrom)){ 
      names(chrom) <- object$snpdf$snp
    }
    if(!is.null(genome.pos) & length(genome.pos) != length(object$snpdf$snp)){ 
      stop("Check 'genome.pos' input. The vector should have length equal to the number of markers in the updog output.")
    } else if(!is.null(genome.pos)){
      names(genome.pos) <- object$snpdf$snp
    }
    ## dosage info
    if(filter.non.conforming){
      geno.dose = matrix(NA,1,1)      
    } else {
      geno.dose <- dist_prob_to_class(geno = geno, prob.thres = prob.thres)
      if(geno.dose$flag)
      {
        geno <- geno.dose$geno
        geno.dose <- geno.dose$geno.dose
      } else {
        geno.dose <- geno.dose$geno.dose
      }
      geno.dose[is.na(geno.dose)] <- ploidy + 1
    }
    mrk.names <- unique(geno$mrk)
    ind.names <- unique(geno$ind)
    nphen = 0
    phen = NULL
    res <- structure(list(ploidy = ploidy,
                          n.ind = length(ind.names),
                          n.mrk = length(mrk.names),
                          ind.names = ind.names,
                          mrk.names = mrk.names,
                          dosage.p1 = dosage.p1[mrk.names],
                          dosage.p2 = dosage.p2[mrk.names],
                          chrom = chrom[mrk.names],
                          genome.pos = genome.pos[mrk.names],
                          prob.thres = prob.thres,
                          geno = geno,
                          geno.dose = geno.dose,
                          nphen = nphen,
                          phen = phen,
                          chisq.pval = NULL),
                     class = "mappoly.data")
    if(filter.non.conforming){
      if (verbose) cat("    Filtering non-conforming markers.\n    ...")
      res <- filter_non_conforming_classes(res)
      if (verbose) cat("\n    Performing chi-square test.\n    ...")
      ##Computing chi-square p.values
      Ds <- array(NA, dim = c(ploidy+1, ploidy+1, ploidy+1))
      for(i in 0:ploidy)
        for(j in 0:ploidy)
          Ds[i+1,j+1,] <- segreg_poly(ploidy = ploidy, dP = i, dQ = j)
      Dpop <- cbind(res$dosage.p1, res$dosage.p2)
      M <- t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
      dimnames(M) <- list(res$mrk.names, c(0:ploidy))
      M <- cbind(M, res$geno.dose)
      res$chisq.pval <- apply(M, 1, mrk_chisq_test, ploidy = ploidy)
      if (verbose) cat("\n    Done.\n")
    }
    return(res)
  }
  # # Case 2: polyRAD
  # else if (inherits(object, 'RADdata')){
  #   outfile = paste0(getwd(), '/import_temp')
  #   polyRAD::Export_MAPpoly(object, file = outfile)
  #   res = read_geno_prob(outfile)
  #   return(res)
  # }
  else stop("You must provide an object of class 'multidog' (from package 'updog') in order to continue importing data.")
}
