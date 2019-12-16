#' A function to run a sgRNA quantification algorithm from NGS sample
#' 
#' @param lib_path The path of the FASTA file.
#' @param design A table contains the study design. It must contain `fastq_path` and `sample_name.`
#' @param map_path The path of file contains gene-sgRNA mapping.
#' @param ncores The number that indicates how many processors will be used with a parallelization. 
#'   The parallelization will be enabled if users do not set the parameter as `-1`` 
#'   (it means the full physical cores will be used) or greater than `1`.
#' @param verbose Display some logs during the quantification if it is set to `TRUE`
#' 
#' @importFrom tools file_ext
#' @importFrom readr read_csv read_tsv
#' @importFrom dplyr left_join
#' @importFrom parallel detectCores makeCluster clusterExport clusterApply stopCluster
#' @importFrom R.utils gunzip
#' @return It will return a list, and the list contains three elements. 
#'   The first element (`count') is a data frame contains the result of the quantification for each sample. 
#'   The second element (`total') is a numeric vector contains the total number of reads of each sample.
#'   The last element (`sequence') a data frame contains the sequence of each sgRNA in the library.
#'   
#' @examples
#' library(CB2)
#' library(magrittr)
#' library(tibble)
#' library(dplyr)
#' library(glue)
#' FASTA <- system.file("extdata", "toydata", "small_sample.fasta", package = "CB2")
#' ex_path <- system.file("extdata", "toydata", package = "CB2")
#' 
#' df_design <- tribble(
#'   ~group, ~sample_name,
#'   "Base", "Base1",  
#'   "Base", "Base2", 
#'   "High", "High1",
#'   "High", "High2") %>% 
#'     mutate(fastq_path = glue("{ex_path}/{sample_name}.fastq"))
#' 
#' cb2_count <- run_sgrna_quant(FASTA, df_design)
#'
#' @export
run_sgrna_quant <- function(lib_path, 
                            design, 
                            map_path = NULL,
                            ncores = 1, 
                            verbose = FALSE) {
    # `design` has to be table `design` has to have four columns: sample_name, group, fastq_path
    if(!all(c("sample_name", "fastq_path") %in% colnames(design))) {
        stop("The design data frame should have both `sample_name` and `fastq_path` columns.")    
    }
    
    if(!file.exists(lib_path)) {
        stop("The library annotation file (FASTA) does not exist.")
    }
    
    if(!is.null(map_path) && !file.exists(map_path)) {
        stop("The mapping file does not exist.")
    }
    
    if(!is.null(map_path)&&!tolower(file_ext(map_path)) %in% c("csv", "tsv")) {
        stop("The mapping file should be either CSV or TSV file.")
    }
    if(!all(file.exists(design$fastq_path))) {
        stop("Some of sample FASTQ files does not exist.")
    }
    
    lib_path <- normalizePath(lib_path)
    design$fastq_path <- normalizePath(design$fastq_path)
    
    is_gzipped <- endsWith(tolower(design$fastq_path), ".gz")
    
    quant_ret <- NULL
    fastq_path <- design$fastq_path
    
    for(i in 1:length(is_gzipped)) {
        if(is_gzipped[i]) {
            tmp_path <- tempfile(fileext = '.fastq')
            R.utils::gunzip(fastq_path[i], tmp_path, remove=FALSE)
            fastq_path[i] <- tmp_path
        }
    }
    if(ncores == -1 || ncores >= 2) {
        max_cores <- detectCores(logical = F)
        if(ncores == -1) {
            ncores <- max_cores
        }
        cl <- makeCluster( min(ncores, max_cores) )
        
        clusterExport(cl, "lib_path", envir = environment() )
        clusterExport(cl, "fastq_path", envir = environment() )
        clusterExport(cl, "is_gzipped", envir = environment() )
        clusterExport(cl, "verbose", envir = environment() )
        
        tmp <- clusterApply(cl, x = 1:length(fastq_path), function(i) {
            CB2::quant(lib_path, fastq_path[i], verbose)
        })
        
        stopCluster(cl)
        quant_ret$sgRNA <- tmp[[1]]$sgRNA
        quant_ret$sequence <- tmp[[1]]$sequence
        quant_ret$count <- sapply(tmp, function(x) x$count)
        quant_ret$total <- sapply(tmp, function(x) x$total)
    } else {
        quant_ret <- quant(lib_path, fastq_path)
    }
    
    if(is.null(map_path)) {
        df_count <- as.data.frame(quant_ret$count)
        rownames(df_count) <- quant_ret$sgRNA
        colnames(df_count) <- design$sample_name
    } else {
        df_count <- as.data.frame(quant_ret$count)
        colnames(df_count) <- design$sample_name
        df_count <- cbind(
            data.frame(id = quant_ret$sgRNA),
            df_count
            )
        if(tolower(file_ext(map_path)) == "csv") {
            df_map <- read_csv(map_path)
        } else {
            df_map <- read_tsv(map_path)
        }
        df_count <-
            left_join(df_map,
                      df_count, by = "id")
    }
    
    total <- quant_ret$total
    names(total) <- design$sample_name
    sequence <- data.frame(sgRNA = quant_ret$sgRNA, sequence=quant_ret$sequence)
    list(count = df_count, total = total, sequence = sequence)
}

#' A function to perform a statistical test at a sgRNA-level, deprecated.
#' @param sgcount This data frame contains read counts of sgRNAs for the samples.
#' @param design This table contains study design. It has to contain `group.`
#' @param group_a The first group to be tested.
#' @param group_b The second group to be tested.
#' @param delim The delimiter between a gene name and a sgRNA ID. It will be used if only rownames contains sgRNA ID.
#' @param ge_id The column name of the gene column.
#' @param sg_id The column/columns of sgRNA identifiers.
#' @return A table contains the sgRNA-level test result, and the table contains these columns: 
#' \itemize{
#'  \item `sgRNA': The sgRNA identifier.
#'  \item `gene': The gene is the target of the sgRNA 
#'  \item `n_a': The number of replicates of the first group.
#'  \item `n_b': The number of replicates of the second group.
#'  \item `phat_a': The proportion value of the sgRNA for the first group.
#'  \item `phat_b': The proportion value of the sgRNA for the second group.
#'  \item `vhat_a': The variance of the sgRNA for the first group.
#'  \item `vhat_b': The variance of the sgRNA for the second group.
#'  \item `cpm_a': The mean CPM of the sgRNA within the first group.
#'  \item `cpm_b': The mean CPM of the sgRNA within the second group.
#'  \item `logFC': The log fold change of sgRNA between two groups.
#'  \item `t_value': The value for the t-statistics.
#'  \item `df': The value of the degree of freedom, and will be used to calculate the p-value of the sgRNA.
#'  \item `p_ts': The p-value indicates a difference between the two groups.
#'  \item `p_pa': The p-value indicates enrichment of the first group.
#'  \item `p_pb': The p-value indicates enrichment of the second group.
#'  \item `fdr_ts': The adjusted P-value of `p_ts'.
#'  \item `fdr_pa': The adjusted P-value of `p_pa'.
#'  \item `fdr_pb': The adjusted P-value of `p_pb'.
#' }
#' @export
run_estimation <- function( sgcount, design, 
                            group_a, group_b, 
                            delim = "_",
                            ge_id = NULL,
                            sg_id = NULL) {
    .Deprecated('measure_sgrna_stats')
    measure_sgrna_stats(
        sgcount, design, 
        group_a, group_b, 
        delim,
        ge_id,
        sg_id
    )
}


#' A function to perform a statistical test at a sgRNA-level
#' @param sgcount This data frame contains read counts of sgRNAs for the samples.
#'
#' @param design This table contains study design. It has to contain `group.`
#' @param group_a The first group to be tested.
#' @param group_b The second group to be tested.
#' @param delim The delimiter between a gene name and a sgRNA ID. It will be used if only rownames contains sgRNA ID.
#' @param ge_id The column name of the gene column.
#' @param sg_id The column/columns of sgRNA identifiers.
#'
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @importFrom dplyr arrange_
#' @importFrom stats p.adjust pt
#' @return A table contains the sgRNA-level test result, and the table contains these columns: 
#' \itemize{
#'  \item `sgRNA': The sgRNA identifier.
#'  \item `gene': The gene is the target of the sgRNA 
#'  \item `n_a': The number of replicates of the first group.
#'  \item `n_b': The number of replicates of the second group.
#'  \item `phat_a': The proportion value of the sgRNA for the first group.
#'  \item `phat_b': The proportion value of the sgRNA for the second group.
#'  \item `vhat_a': The variance of the sgRNA for the first group.
#'  \item `vhat_b': The variance of the sgRNA for the second group.
#'  \item `cpm_a': The mean CPM of the sgRNA within the first group.
#'  \item `cpm_b': The mean CPM of the sgRNA within the second group.
#'  \item `logFC': The log fold change of sgRNA between two groups.
#'  \item `t_value': The value for the t-statistics.
#'  \item `df': The value of the degree of freedom, and will be used to calculate the p-value of the sgRNA.
#'  \item `p_ts': The p-value indicates a difference between the two groups.
#'  \item `p_pa': The p-value indicates enrichment of the first group.
#'  \item `p_pb': The p-value indicates enrichment of the second group.
#'  \item `fdr_ts': The adjusted P-value of `p_ts'.
#'  \item `fdr_pa': The adjusted P-value of `p_pa'.
#'  \item `fdr_pb': The adjusted P-value of `p_pb'.
#' }
#' @examples 
#' library(CB2)
#' data(Evers_CRISPRn_RT112)
#' measure_sgrna_stats(Evers_CRISPRn_RT112$count, Evers_CRISPRn_RT112$design, "before", "after")
#' 
#' @export
measure_sgrna_stats <- function(sgcount, design, 
                                group_a, group_b, 
                                delim = "_",
                                ge_id = NULL,
                                sg_id = NULL) {

    if(!is.data.frame(sgcount) && !is.matrix(sgcount)) {
        stop("sgcount has to be either a data.frame or a matrix.")    
    }
    
    if(!all(c("sample_name", "group") %in% colnames(design))) {
        stop("The design table should have both `sample_name` and `group` columns.")    
    }
    
    if(!all(design$sample_name %in% colnames(sgcount))) {
        stop("Some samples are missing in sgcount.")
    }
    
    if(!all(group_a %in% design$group)) {
        stop("group_a must be one of the groups in the design data frame.")
    }
    if(!all(group_b %in% design$group)) {
        stop("group_b must be one of the groups in the design data frame.")
    }

    if(is.matrix(sgcount)) {
        if(!is.null(ge_id)||!is.null(sg_id)) {
            stop("ge_id and sg_id should be NULL if sgcount is a matrix.")
        }
    }
    
    if(xor(is.null(ge_id), is.null(sg_id))) {
        stop("Both of ge_id and sg_id should be null or non-null.")    
    }
    
    cname <- colnames(sgcount)
    if(is.null(ge_id)) {
        if(!all(sapply(sgcount, class) == "numeric")) {
            stop(
                paste0("sgcount contains some character columns. ", 
                 "It may need to specify both ge_id and sg_id."))
        }
        sgcount <- as.matrix(sgcount)
        cnt_delim <- stringr::str_count(rownames(sgcount), delim)
        if(!all(cnt_delim == 1)) {
          stop("Every rownames should contains exact one delimiter.")
        }
    }
    else {
        if(length(ge_id) != 1) {
            stop("ge_id should be a character variables.")
        }
        if(!(ge_id %in% cname)) {
            stop("Column of ge_id was not found in sgcount.")
        }
        if(!all(sg_id %in% cname)) {
            stop("Column/columns of sg_id was not found in sgcount.")
        }
    }

    group_a <- which(cname %in% design$sample_name[design$group == group_a])
    group_b <- which(cname %in% design$sample_name[design$group == group_b])
    sgcount_a <- as.matrix(sgcount[, group_a])
    sgcount_b <- as.matrix(sgcount[, group_b])
    
    nmat_a <- rep.row(colSums(sgcount_a), nrow(sgcount_a))
    nmat_b <- rep.row(colSums(sgcount_b), nrow(sgcount_b))
    est_a <- fit_ab(sgcount_a, nmat_a)
    est_b <- fit_ab(sgcount_b, nmat_b)
    
    # if you have gene colum then get read of this part
    
    if(is.null(ge_id)) {
        est <- data.frame(sgRNA = rownames(sgcount), stringsAsFactors=F) %>% as_tibble()
        est$gene <- stringr::str_split(est$sgRNA, "_", simplify = T)[, 1]
    } else {
        est <- data.frame(gene=sgcount[,ge_id])
        colnames(est)[1] <- "gene"
        est <- cbind(est, as.data.frame(sgcount[,sg_id]))
        est <- data.frame(gene = sgcount[,ge_id])
        est <- cbind(est, sgcount[,sg_id])
    }
    est$n_a <- length(group_a)
    est$n_b <- length(group_b)
    est$phat_a <- est_a$phat
    est$vhat_a <- est_a$vhat
    est$phat_b <- est_b$phat
    est$vhat_b <- est_b$vhat
    est$cpm_a <- rowMeans(sgcount_a/nmat_a * 10^6)
    est$cpm_b <- rowMeans(sgcount_b/nmat_b * 10^6)
    
    # if `ncol(sgcount_a)' or `nocl(sgcount_b)' are 1, `vhat_a' or `vhat_b' will only contain `na'. 
    # In this case we need to treat all the `na' values as 0
    # for further procedure.
    est$vhat_a[is.na(est$vhat_a)] <- 0
    est$vhat_b[is.na(est$vhat_b)] <- 0
    
    est$logFC <- log2(est$cpm_b + 1) - log2(est$cpm_a + 1)
    zero_var <- 1 * ((est$vhat_a == 0) & (est$vhat_b == 0))
    eps <- .Machine$double.eps
    est$t_value <- (est$phat_b - est$phat_a)/sqrt(est$vhat_a + est$vhat_b + eps * zero_var)
    
    
    est$df <- ((est$vhat_a + est$vhat_b)^2 + zero_var)/((est$vhat_a^2)/max(1,est$n_a - 1) + (est$vhat_b^2)/max(1,est$n_b - 1) + zero_var)
    est$df[est$df == 0] <- 1
    est$df[is.nan(est$df)] <- 1
    
    est$p_ts <- 2 * pt(-abs(est$t_value), df = est$df)
    est$p_pa <- pt(est$t_value, df = est$df)
    est$p_pb <- pt(-est$t_value, df = est$df)
    
    est$fdr_ts <- p.adjust(est$p_ts, method = "fdr")
    est$fdr_pa <- p.adjust(est$p_pa, method = "fdr")
    est$fdr_pb <- p.adjust(est$p_pb, method = "fdr")
    est
}

#' A function to perform gene-level test using a sgRNA-level statistics.
#' 
#' @param sgrna_stat A data frame created by `measure_sgrna_stats'
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom dplyr group_by
#' @importFrom stats p.adjust
#'
#' @return A table contains the gene-level test result, and the table contains these columns: 
#' \itemize{
#'   \item `gene': Theg gene name to be tested.
#'   \item `n_sgrna': The number of sgRNA targets the gene in the library.
#'   \item `cpm_a': The mean of CPM of sgRNAs within the first group.
#'   \item `cpm_b': The mean of CPM of sgRNAs within the second group.
#'   \item `logFC': The log fold change of the gene between two groups.
#'   \item `p_ts': The p-value indicates a difference between the two groups at the gene-level.
#'   \item `p_pa': The p-value indicates enrichment of the first group at the gene-level.
#'   \item `p_pb': The p-value indicates enrichment of the second group at the gene-level.
#'   \item `fdr_ts': The adjusted P-value of `p_ts'.
#'   \item `fdr_pa': The adjusted P-value of `p_pa'.
#'   \item `fdr_pb': The adjusted P-value of `p_pb'.
#' }
#' 
#' @examples
#' data(Evers_CRISPRn_RT112)
#' measure_gene_stats(Evers_CRISPRn_RT112$sg_stat)
#' 
#' @export
measure_gene_stats <- function(sgrna_stat) {
    if(!all(c("gene", "cpm_a", "cpm_b", "logFC", "p_ts", "p_pa", "p_pb") %in% colnames(sgrna_stat))) {
        stop("It looks like `sgrna_stat` does not contain any result of a statistical test.")
    }
    sgrna_stat %>% dplyr::group_by_(~gene) %>%
        dplyr::summarise_(
            n_sgrna = ~dplyr::n(),
            cpm_a = ~mean(cpm_a),
            cpm_b = ~mean(cpm_b),
            logFC = ~mean(logFC),
            p_ts = ~ifelse(dplyr::n() > 1, metap::sumlog(p_ts)$p, sum(p_ts)),
            p_pa = ~ifelse(dplyr::n() > 1, metap::sumlog(p_pa)$p, sum(p_pa)),
            p_pb = ~ifelse(dplyr::n() > 1, metap::sumlog(p_pb)$p, sum(p_pb))
        ) %>%
        dplyr::ungroup() %>% dplyr::mutate_(
            fdr_ts = ~p.adjust(p_ts, method = "fdr"),
            fdr_pa = ~p.adjust(p_pa, method = "fdr"),
            fdr_pb = ~p.adjust(p_pb, method = "fdr")
        )
}

