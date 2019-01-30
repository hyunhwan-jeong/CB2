context("Check whether it works")

library(CB2)
library(tidyverse)
FASTA <- system.file("extdata",
                     "toydata", "small_sample.fasta",
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
  expect_length((cb2_count <- run_sgrna_quant(FASTA, df_design)), 2)
  expect_error(run_sgrna_quant(FASTA, df_design[,c("group", "sample_name")]))
  expect_length((sgstat <- run_estimation(cb2_count$count, df_design, "Base", "Low")), 19)
  expect_error(run_estimation(cb2_count$count, df_design, "Base", "Lows"),  "group_b must be one of groups in the design data.frame.")
  expect_error(run_estimation(cb2_count$count, df_design, "Bae", "Low"),  "group_a must be one of groups in the design data.frame.")  
  expect_error(run_estimation(cb2_count$count[,-1], df_design, "Bae", "Low"),  "Some samples are missing in sgcount.")  
  expect_length(measure_gene_stats(sgstat), 11)
  expect_error(measure_gene_stats(cb2_count$count), "It seems that the `sgrna_stat` does not contain any result of statistical test.")
})
