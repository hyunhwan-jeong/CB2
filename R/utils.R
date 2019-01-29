#' @param sgcount The input data.frame contains counts of sgRNAs for each sample
#'
#' A function to calculate the CPM (Counts Per Million) (required)
#' 
#' @return a normalized CPM data.frame will be returned
#' 
#' @examples
#' 
#' data(Evers$count)
#' get_CPM(Evers$count)
#' 
#' @export
get_CPM <- function(sgcount) {
  nmat <- rep.row(colSums(sgcount), nrow(sgcount))
  sgcount / nmat * 10^6
}

#' @param sgcount The input matrix contains counts of sgRNAs for each sample (required)
#'
#' @param df_design The data.frame contains study design(required)
#'
#' @export
plot_PCA <- function(sgcount, df_design) {
  pca_obj <- sgcount %>% t %>% prcomp 
  pca_obj$x %>% as.data.frame() %>% 
    tibble::rownames_to_column("sample_name") %>% 
    dplyr::left_join(df_design, by = "sample_name") %>% 
    ggplot2::ggplot(ggplot2::aes(x=PC1, y=PC2)) + 
    ggplot2::geom_point(ggplot2::aes(color=group), size=2) +
    ggplot2::geom_text(ggplot2::aes(label=sample_name))
}

#' @param count_obj 
#'
#' @param df_design 
#'
#' @export
calc_mappability <- function(count_obj, df_design) {
  csum <- count_obj$count %>% colSums()
  mp <- csum / count_obj$total * 100.0
  df_design %>% 
    dplyr::mutate(mappability =  mp) %>% 
    dplyr::select(-fastq_path)
}

#' @param sgcount 
#'
#' @param df_design 
#'
#' @export
plot_count_distribution <- function(sgcount, df_design) {
  sgcount %>% tibble::rownames_to_column("sgRNA") %>% 
    tidyr::gather(sample_name, count, -sgRNA) %>% 
    dplyr::left_join(df_design, by = "sample_name") %>% 
    ggplot2::ggplot(ggplot2::aes(log2(1+count))) + 
    ggplot2::geom_density(aes(fill=group)) +
    ggplot2::facet_wrap(~sample_name, ncol=1)
}

#' @param sgcount 
#'
#' @param df_design 
#' @param gene 
#'
#' @export
plot_dotplot <- function(sgcount, df_design, gene) {
  sgcount %>% tibble::rownames_to_column("sgRNA") %>% 
    tidyr::gather(sample_name, count, -sgRNA) %>% 
    dplyr::left_join(df_design, by = "sample_name") %>% 
    dplyr::filter(stringr::str_detect(sgRNA, glue::glue("^{gene}")))%>% 
  ggplot2::ggplot(ggplot2::aes(x=group, y=count)) +
    ggplot2::geom_dotplot(ggplot2::aes(fill=group, color=group),
                 binaxis='y',
                 stackdir='center',
                 stackratio=1.5,
                 dotsize=1.2) + 
    ggplot2::facet_wrap(~sgRNA, scales="free_y")
}
  