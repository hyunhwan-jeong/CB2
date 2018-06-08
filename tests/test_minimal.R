FASTA <- system.file("extdata", "small_sample.fasta", package = "CRISPRExpress")
FASTQ <- c()
for(g in c("Low", "High", "Base")) {
  for(i in 1:2) {
    FASTQ <- c(FASTQ, system.file("extdata", 
                                  sprintf("%s%d.fastq", g, i), 
                                  package = "CRISPRExpress"))
  }
}
library(tidyverse)
sample_quant_output <- quant(FASTA, FASTQ)
rep.row<-function(x,n) matrix(rep(x,each=n),nrow=n)
xvec <- sample_quant_output$count 
xvec <- matrix(sample(1:10000, 1e5*10, replace=T), ncol=10)
xvec[1,] <- 0
xvec
nvec <- xvec %>% colSums() %>% rep.row(nrow(xvec))
as.data.frame(fit_ab(xvec, nvec))
