#' A benchmark CRISPRn pooled screen data from Evers et al.
#'
#' @format The data object is a list and contains below information:
#' \describe{
#'   \item{count}{The count matrix from Evers et al.'s paper and contains the CRISPRn screening result using RT112 cell-line.
#'     It contains three different replicates for T0 (before) and contains different three replicates for T1 (after).}
#'   \item{egenes}{The list of 46 essential genes used in Evers et al.'s study.}
#'   \item{ngenes}{The list of 47 non-essential genes used in Evers et al.'s study.}
#'   \item{design}{The data.frame contains study design.}
#'   \item{sg_stat}{The data.frame contains the sgRNA-level statistics.}
#'   \item{gene_stat}{The data.frame contains the gene-level statistics.}
#' }
#' 
#' @docType data
#'
#' @usage data(Evers_CRISPRn_RT112)
#'
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/27111720/}
"Evers_CRISPRn_RT112"


#' A benchmark CRISPRn pooled screen data from Sanson et al. 
#'
#' @format The data object is a list and contains below information:
#' \describe{
#'   \item{count}{The count matrix from Sanson et al.'s paper and contains the CRISPRn screening result using A375 cell-line. 
#'     It contains a sample of plasimd, and three biological replicates after three weeks.}
#'   \item{egenes}{The list of 1,580 essential genes used in Sanson et al.'s study.}
#'   \item{ngenes}{The list of 927 non-essential genes used in Sanson et al.'s study.}
#'   \item{design}{The data.frame contains study design.}
#' }
#' 
#' @docType data
#'
#' @usage data(Sanson_CRISPRn_A375)
#'
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/30575746/}
"Sanson_CRISPRn_A375"
