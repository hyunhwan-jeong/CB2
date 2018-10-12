rep.row <- function(x, n) matrix(rep(x,each=n),nrow=n)

#' @export
run_sgrna_quant <- function(lib_path, design) {
  # `design` has to be data.frame
  # `design` has to have four columns: sample_name, group, fastq_path
  lib_path <- normalizePath(lib_path)
  design$fastq_path <- normalizePath(design$fastq_path)
  quant_ret <- quant(lib_path, design$fastq_path)
  df_count <- as.data.frame(quant_ret$count)
  rownames(df_count) <- quant_ret$sgRNA
  colnames(df_count) <- design$sample_name
  df_count
}

read_sgrna_mat <- function(mat_file) {
  
}

#' @importFrom magrittr %>%
#' @export %>%
#' @importFrom tibble tibble
#' @export tibble
#' @importFrom dplyr arrange
#' @export arrange
#' @export
run_estimation <- function(sgcount, design, group_a, group_b) {
  group_a <- design$sample_name[design$group == group_a]
  group_b <- design$sample_name[design$group == group_b]
  sgcount_a <- as.matrix(sgcount[,group_a])
  sgcount_b <- as.matrix(sgcount[,group_b])
  nmat_a <- rep.row(colSums(sgcount_a), nrow(sgcount_a))
  nmat_b <- rep.row(colSums(sgcount_b), nrow(sgcount_b))
  est_a <- fit_ab(sgcount_a, nmat_a)
  est_b <- fit_ab(sgcount_b, nmat_b)
  
  # if you have gene colum then get read of this part
  est <- tibble(sgRNA = rownames(sgcount))
  est$gene <- stringr::str_split(est$sgRNA, "_", simplify=T)[,1]
  est$n_a <- length(group_a)
  est$n_b <- length(group_b)
  est$phat_a <- est_a$phat
  est$vhat_a <- est_a$vhat
  est$phat_b <- est_b$phat
  est$vhat_b <- est_b$vhat
  est$cpm_a <- rowMeans(sgcount_a / nmat_a * 10^6)
  est$cpm_b <- rowMeans(sgcount_b / nmat_b * 10^6)
  est$logFC <- log2(est$cpm_a+1) 
  - log2(est$cpm_b+1)
  zero_var <- 1 * ((est$vhat_a == 0) & (est$vhat_b == 0))
  eps <- .Machine$double.eps
  est$t_value <- (est$phat_b - est$phat_a) / sqrt(est$vhat_a + est$vhat_b + eps * zero_var)
  
  est$df <- ((est$vhat_a+est$vhat_b)^2+zero_var)/ ((est$vhat_a^2) / (est$n_a-1) + (est$vhat_b^2) / (est$n_b-1) + zero_var)
  est$df[est$df==0] <- 1
  est$df[is.nan(est$df)] <- 1
  
  est$p_ts <- 2 * pt(-abs(est$t_value), df=est$df)
  est$p_pa <- pt(est$t_value, df=est$df)
  est$p_pb <- pt(-est$t_value, df=est$df)
  
  est$fdr_ts <- p.adjust(est$p_ts, method="fdr")
  est$fdr_pa <- p.adjust(est$p_pa, method="fdr")
  est$fdr_pb <- p.adjust(est$p_pb, method="fdr")
  est %>% arrange(sgRNA)
}

#' @importFrom magrittr %>%
#' @export %>%
#' @importFrom tibble tibble
#' @export tibble
#' @importFrom dplyr group_by
#' @export group_by
#' @export
measure_gene_stats <-
  function(sgrna_stat) {
    sgrna_stat %>% 
      dplyr::group_by(gene) %>%
      dplyr::summarise(
        n_sgrna = n(),
        cpm_a = mean(cpm_a),
        cpm_b = mean(cpm_b),
        logFC = mean(logFC),
        p_ts = metap::sumlog(p_ts)$p,
        p_pa = metap::sumlog(p_pa)$p,
        p_pb = metap::sumlog(p_pb)$p
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        fdr_ts = p.adjust(p_ts, method="fdr"),
        fdr_pa = p.adjust(p_pa, method="fdr"),
        fdr_pb = p.adjust(p_pb, method="fdr")
      )
  }

