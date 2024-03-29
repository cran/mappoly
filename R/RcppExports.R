# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

calc_genoprob_cpp <- function(m, geno, ph1, ph2, rf, probs, verbose) {
    .Call('_mappoly_calc_genoprob_cpp', PACKAGE = 'mappoly', m, geno, ph1, ph2, rf, probs, verbose)
}

calc_genprob_haplo_cpp <- function(m, n_mrk, n_ind, haplo, emit, rf, probs, verbose) {
    .Call('_mappoly_calc_genprob_haplo_cpp', PACKAGE = 'mappoly', m, n_mrk, n_ind, haplo, emit, rf, probs, verbose)
}

calc_genprob_haplo_highprec_cpp <- function(m, n_mrk, n_ind, haplo, emit, rf, probs, verbose) {
    .Call('_mappoly_calc_genprob_haplo_highprec_cpp', PACKAGE = 'mappoly', m, n_mrk, n_ind, haplo, emit, rf, probs, verbose)
}

loglike_hmm_cpp <- function(m, geno, ph1, ph2, rf, verbose) {
    .Call('_mappoly_loglike_hmm_cpp', PACKAGE = 'mappoly', m, geno, ph1, ph2, rf, verbose)
}

.vcf_get_probabilities <- function(mat, pl_pos) {
    .Call('_mappoly_vcf_get_probabilities', PACKAGE = 'mappoly', mat, pl_pos)
}

.vcf_transform_dosage <- function(mat, gt_pos) {
    .Call('_mappoly_vcf_transform_dosage', PACKAGE = 'mappoly', mat, gt_pos)
}

.vcf_get_ploidy <- function(mat, gt_pos) {
    .Call('_mappoly_vcf_get_ploidy', PACKAGE = 'mappoly', mat, gt_pos)
}

.vcf_get_depth <- function(mat, dp_pos) {
    .Call('_mappoly_vcf_get_depth', PACKAGE = 'mappoly', mat, dp_pos)
}

