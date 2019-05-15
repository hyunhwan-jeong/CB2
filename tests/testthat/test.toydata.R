context("Check whether it works")

library(CB2)
library(dplyr)
library(readr)
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
  expect_length((cb2_count <- run_sgrna_quant(FASTA, df_design)), 3)
  expect_error(run_sgrna_quant(FASTA, df_design[,c("group", "sample_name")]))
  expect_length((sgstat <- run_estimation(cb2_count$count, df_design, "Base", "Low")), 19)
  expect_error(run_estimation(cb2_count$count, df_design, "Base", "Lows"),  "group_b must be one of the groups in the design data frame.")
  expect_error(run_estimation(cb2_count$count, df_design, "Bae", "Low"),  "group_a must be one of the groups in the design data frame.")  
  expect_error(run_estimation(cb2_count$count[,-1], df_design, "Bae", "Low"),  "Some samples are missing in sgcount.")  
  expect_length(measure_gene_stats(sgstat), 11)
  expect_error(measure_gene_stats(cb2_count$count), "It looks like `sgrna_stat` does not contain any result of a statistical test.")
})


test_that("Testing whether subsamplings with the same seed are identical.", {
  set.seed(123)
  cb2_count <- run_sgrna_quant(FASTA, df_design, sampling_ratio = 0.5)
  
  set.seed(123)
  cb2_count2 <- run_sgrna_quant(FASTA, df_design, sampling_ratio = 0.5)
  
  expect_identical(cb2_count, cb2_count2)
})

test_that("Testing whether the subsampling code works as intended.", {
  cb2_count_nosubs <- run_sgrna_quant(FASTA, df_design)

  for( i in seq(0.1, 0.9, 0.1)) {
    set.seed(123)
    cb2_count <- run_sgrna_quant(FASTA, df_design, sampling_ratio = i)
    expect_equal(sum(dplyr::between(cb2_count$total / cb2_count_nosubs$total, i-0.1, i+0.1)), length(cb2_count$total))
  }
})

test_that("Testing the result are consistent to the publication", {
  data("Sanson_CRISPRn_A375")
  published <- read.csv("https://raw.githubusercontent.com/hyunhwaj/CB2-Experiments/master/01_gene-level-analysis/results/Sanson/CRISPRn-A375/FDR/CB2.csv")
  sgrna <- run_estimation(Sanson_CRISPRn_A375$count, Sanson_CRISPRn_A375$design, "ctl", "trt")
  gene <- measure_gene_stats(sgrna)
  
  test <- dplyr::left_join(
    published,
    gene %>% select(gene, fdr_test = fdr_pa)
  )
    
  testthat::expect(sum(abs(test$fdr-test$fdr_test)) < 1e-10)
})
