context("Check whether it works")

library(CB2)
library(dplyr)
library(tidyr)
library(tibble)
FASTA <- system.file("extdata",
                     "toydata", "small_sample_dup.fasta",
                     package = "CB2")
df_design <- data.frame()
for(g in c("Low", "High", "Base")) {
  for(i in 1:2) {
    FASTQ <- system.file("extdata", "toydata",
                         sprintf("%s%d.fastq", g, i), 
                         package = "CB2")
    df_design <- rbind(df_design, 
                       data.frame(
                         group = g, 
                         sample_name = sprintf("%s%d", g, i),
                         fastq_path = FASTQ, 
                         stringsAsFactors = F)
    )
  }
}

rep.row <- function(x,n) matrix(rep(x,each=n),nrow=n)


test_that("A very simple sanity check", {
  expect_length((cb2_count <- run_sgrna_quant(FASTA, df_design)), 3)
  expect_error(run_sgrna_quant(FASTA, df_design[,c("group", "sample_name")]))
  expect_length((sgstat <- measure_sgrna_stats(cb2_count$count, df_design, "Base", "Low")), 19)
  expect_error(measure_sgrna_stats(cb2_count$count, df_design, "Base", "Lows"),  "group_b must be one of the groups in the design data frame.")
  expect_error(measure_sgrna_stats(cb2_count$count, df_design, "Bae", "Low"),  "group_a must be one of the groups in the design data frame.")  
  expect_error(measure_sgrna_stats(cb2_count$count[,-1], df_design, "Bae", "Low"),  "Some samples are missing in sgcount.")  
  expect_length(measure_gene_stats(sgstat), 11)
  expect_error(measure_gene_stats(cb2_count$count), "It looks like `sgrna_stat` does not contain any result of a statistical test.")
  
  expect_error(measure_sgrna_stats(cb2_count, df_design, "Base", "Low"),  "sgcount has to be either a data.frame or a matrix.")
  
})


# test_that("The result has to be consistent to the publication", {
#   data("Sanson_CRISPRn_A375")
#   published <- read.csv("https://raw.githubusercontent.com/hyunhwaj/CB2-Experiments/master/01_gene-level-analysis/results/Sanson/CRISPRn-A375/FDR/CB2.csv")
#   sgrna <- measure_sgrna_stats(Sanson_CRISPRn_A375$count, Sanson_CRISPRn_A375$design, "ctl", "trt")
#   gene <- measure_gene_stats(sgrna)
# 
#   test <- dplyr::left_join(
#     published,
#     gene %>% select(gene, fdr_test = fdr_pa)
#   )
# 
#   testthat::expect(sum(abs(test$fdr-test$fdr_test)) < 1e-10)
# })

test_that("It should allow a table with two columns contains sgRNA information", {
  data(Evers_CRISPRn_RT112)
  df_cb2 <- separate(rownames_to_column(Evers_CRISPRn_RT112$count, "gRNA"), "gRNA", c("gene", "id"))
  testthat::expect_type(
    measure_sgrna_stats(df_cb2, Evers_CRISPRn_RT112$design, "before", "after", ge_id = "gene", sg_id = "id"),
    "list"
  )
})


test_that("Testing error handling of the measure_sgrna_stats.", {
  data(Evers_CRISPRn_RT112)
  df_cb2 <- separate(rownames_to_column(Evers_CRISPRn_RT112$count, "gRNA"), "gRNA", c("gene", "id"))
  
  
  expect_error(measure_sgrna_stats(df_cb2, Evers_CRISPRn_RT112$design, "before", "after"),  
               paste0("sgcount contains some character columns. ", 
                      "It may need to specify both ge_id and sg_id."))

  expect_error(measure_sgrna_stats(Evers_CRISPRn_RT112$count, 
                                   Evers_CRISPRn_RT112$design, 
                                   "before", "after", "-"),  
               "Every rownames should contains exact one delimiter.")
  
  expect_error(measure_sgrna_stats(Evers_CRISPRn_RT112$count, 
                                   Evers_CRISPRn_RT112$design, 
                                   "before", "after", 
                                   ge_id = "A"),  
               "Both of ge_id and sg_id should be null or non-null.")
  
  expect_error(measure_sgrna_stats(df_cb2, Evers_CRISPRn_RT112$design, 
                                   "before", "after", ge_id = c("A", "B"),
                                   sg_id = c("A", "B")),  
               "ge_id should be a character variables")
  
  
})

