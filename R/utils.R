rep.row <- function(x, n) matrix(rep(x, each = n), nrow = n)

#' A function to normalize sgRNA read counts.
#' 
#' @param sgcount The input table contains read counts of sgRNAs for each sample
#'
#' A function to calculate the CPM (Counts Per Million) (required)
#' 
#' @return a normalized CPM table will be returned
#' 
#' @examples
#' library(CB2)
#' data(Evers_CRISPRn_RT112)
#' get_CPM(Evers_CRISPRn_RT112$count)
#' 
#' @export
get_CPM <- function(sgcount) {
  cols <- which(sapply(sgcount, class) == "numeric")
  nmat <- rep.row(colSums(sgcount[,cols]), nrow(sgcount[,cols]))
  sgcount[,cols] <- sgcount[,cols]/nmat * 10^6
  sgcount
}

#' A function to calculate the mappabilities of each NGS sample.
#' 
#' @param count_obj A list object is created by `run_sgrna_quant`.
#' @param df_design The table contains a study design.
#' @importFrom magrittr %>%
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
#' calc_mappability(cb2_count, df_design)
#'
#' @export
calc_mappability <- function(count_obj, df_design) {
  csum <- count_obj$count %>% colSums()
  mp <- csum/count_obj$total * 100
  df_design %>% 
    dplyr::mutate_(total_reads = ~count_obj$total) %>% 
    dplyr::mutate_(mapped_reads = ~csum) %>% 
    dplyr::mutate_(mappability = ~mp) %>% 
    dplyr::select_(.dots = c("-fastq_path"))
}

#' A function to join a count table and a design table.
#'
#' @param sgcount The input matrix contains read counts of sgRNAs for each sample.
#' @param df_design The table contains a study design.
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @return A tall-thin and combined table of the sgRNA read counts and study design will be returned.
#'
#' @examples 
#' library(CB2)
#' data(Evers_CRISPRn_RT112) 
#' head(join_count_and_design(Evers_CRISPRn_RT112$count, Evers_CRISPRn_RT112$design))
#' 
#' @export
join_count_and_design <- function(sgcount, df_design) {
  cols <- colnames(sgcount)
  
  if(all(sapply(sgcount, class) == "numeric")) {
    sgcount %>% as.data.frame(stringsAsFactors=F) %>% 
      tibble::rownames_to_column("sgRNA") %>% 
      tidyr::gather_(key_col = "sample_name", value_col = "count", 
                     gather_cols = cols) %>% 
      dplyr::left_join(df_design, by = "sample_name") 
  } else {
    cols <- colnames(sgcount)[sapply(sgcount, class) == "numeric"]
    sgcount %>% as.data.frame(stringsAsFactors=F) %>% 
      tidyr::gather_(key_col = "sample_name", value_col = "count", 
                     gather_cols = cols) %>% 
      dplyr::left_join(df_design, by = "sample_name") 
  }
}

library(magrittr)
library(readr)
library(dplyr)
library(stringr)
input_path <- "~/Downloads/broadgpp-dolcetto-targets-seta.txt"
out_path <- "~/Downloads/broadgpp-dolcetto-targets-seta.fa"


#' Convert a CSV/TSV format file to a FASTA file. 
#'
#' @param input_path A path to a CSV/TSV (comma/tab-separated values) format file. 
#' @param output_path A path to the output file.
#' @param seq_col A column index or column name whose column contains guide RNA sequences.
#' @param id_col A column index or column name whose column contains gene/guide name.  
#' @param id_type Type of `id_col`, it can be specified either `gene` or `sgRNA`.
#'
#' @return It does not return anything. Only create a FASTA file at `output_path`.
#' 
#' @import dplyr
#' @import readr
#' @importFrom magrittr %>%
#' @importFrom stringr str_c
#' @export
#'
#' @examples
#' library(CB2)
#' 
#' fa_path <- system.file(
#'  "extdata",
#'  "toydata",
#'  "small_sample.fasta",
#'  package = "CB2"
#' ) 
#' csv_path <- system.file(
#'   "extdata",
#'   "toydata",
#'   "small_sample.csv",
#'   package = "CB2"
#' )
#' 
#' tmp_fasta_path <- tempfile(fileext = ".fa")
#' table2fa(
#'   input_path = csv_path,
#'   output_path = tmp_fasta_path,
#'   id_col = 1,
#'   seq_col = 2,
#'   id_type = "sgRNA"
#' )
#' cat(readLines(tmp_fasta_path)[1:5], sep="\n") # To see the head of the file.
table2fa <- function(input_path, 
                     output_path,
                     seq_col = 1,
                     id_col = 2,
                     id_type = "gene") {
  if(endsWith(tolower(input_path), ".csv")) {
    df_lib <- read_csv(input_path)
  } else {
    df_lib <- read_tsv(input_path)
  }
  
  if(typeof(seq_col) != typeof(id_col)) {
    stop("Both `seq_col` and `id_col` have to be the same type.")
  }
  
  message(typeof(seq_col), " type is given for column indexing.")
  
  if(typeof(seq_col) == "character") {
    vars <- c(id_col, seq_col)
    if(!all(vars %in% colnames(df_lib))) {
      stop("Either `seq_col` or `id_col` is not existed.")
    }
  } else {
    vars <- c(id_col, seq_col)
    if(!all(vars >= 1) || !all(vars <= ncol(df_lib))) {
      stop("Either `seq_col` or `id_col` is out of the range.")
    }
    vars <- colnames(df_lib)[vars]
  }
  
  if(id_type == "gene") {
    df_sel <- df_lib %>% 
      select(all_of(vars)) %>%
      select(id = 1, seq = 2) %>% 
      mutate(line = str_c(">",
                          .data$id,
                          "_",
                          .data$seq,
                          "\n",
                          .data$seq))
  } else {
    df_sel <- df_lib %>% 
      select(all_of(vars)) %>%
      select(id = 1, seq = 2) %>% 
      mutate(line = str_c(">",
                          .data$id,
                          "\n",
                          .data$seq))
  }
  
  cat(file=output_path, df_sel$line, sep="\n")
  message("A FASTA file has been created at ", out_path)
  NULL
}

#' Convert a FASTA file to a data.frame
#'
#' @param fasta_path a path to the fasta file.
#'
#' @return It returns a data.frame that contains information of guide RNA library. 
#' @export
#'
#' @examples
#' library(CB2)
#' 
#' fa_path <- system.file(
#'  "extdata",
#'  "toydata",
#'  "small_sample.fasta",
#'  package = "CB2"
#' ) 
#' cat(readLines(fa_path)[1:5], sep="\n") # To see the head of the file.
#' 
#' fasta2df(fa_path)
fasta2df <- function(fasta_path) {
  ids <- readLines(fasta_path)[c(TRUE, FALSE)]
  seqs <- readLines(fasta_path)[c(FALSE, TRUE)]
  
  message("Performing sanity checks.")
  
  valid_ids <- startsWith(ids, ">")
  if(!all(valid_ids)) {
    stop("Something is not correct with the guide IDs. Please check the line ",
         which(!valid_ids)[1] * 2 -1)
    
  }
  
  valid_seqs <- sapply(seqs, function(x) nchar(x) == nchar(seqs[1]))
  
  if(!all(valid_seqs)) {
    stop("Something is not correct with the guide RNAs. Please check the line ",
         which(!valid_seqs)[1] * 2)
  }
  
  ids <- stringr::str_remove(ids, ">")  
  
  message("Generating a data.frame")
  data.frame(
    id = ids,
    sgRNA = seqs
  ) 
}
