library(tidyverse)
Sanson_CRISPRn_A375 <- list()
Sanson_CRISPRn_A375$count <- 
  read_tsv("https://raw.githubusercontent.com/hyunhwaj/CB2-Experiments/master/01_gene-level-analysis/data/Sanson/CRISPRn-A375.tsv") %>% 
  unite(sgRNA, c("gene", "sgRNA"), sep="_") %>% 
  column_to_rownames("sgRNA")

Sanson_CRISPRn_A375$egenes <- scan("https://raw.githubusercontent.com/hyunhwaj/CB2-Experiments/master/01_gene-level-analysis/data/Sanson/essential-genes.txt", what="character")
Sanson_CRISPRn_A375$ngenes <- scan("https://raw.githubusercontent.com/hyunhwaj/CB2-Experiments/master/01_gene-level-analysis/data/Sanson/non-essential-genes.txt", what="character")

Sanson_CRISPRn_A375$design <- 
  tribble(
    ~group, ~sample_name,
    "ctl", "pDNA",
    "trt", "RepA",
    "trt", "RepB",
    "trt", "RepC"
  )

usethis::use_data(Sanson_CRISPRn_A375)

