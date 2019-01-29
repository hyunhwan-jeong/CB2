#' @param sgcount 
#'
#' @export
get_CPM <- function(sgcount) {
  nmat <- rep.row(colSums(sgcount), nrow(sgcount))
  rowMeans(sgcount / nmat * 10^6)
}