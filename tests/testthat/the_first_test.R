context("Check whether it works")


library(tidyverse)
FASTA <- system.file("extdata", "small_sample.fasta", package = "CRISPRExpress")
FASTQ <- c()
for(g in c("Base", "High", "Low")) {
  for(i in 1:2) {
    FASTQ <- c(FASTQ, system.file("extdata", sprintf("%s%d.fastq", g, i), package="CRISPRExpress"))
  }
}

rep.row <- function(x,n) matrix(rep(x,each=n),nrow=n)


test_that("The first test", {
  lst.quant <- quant(FASTA, FASTQ) %>%
    expect_success()
  
  fit_ab(lst.quant$count, 
         lst.quant$count %>% 
           colSums %>% 
           rep.row( nrow(lst.quant$count))) %>%
    expect_success()
  
})