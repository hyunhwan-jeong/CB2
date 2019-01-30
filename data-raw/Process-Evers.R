library(tidyverse)
library(CB2)
Evers_CRISPRn_RT112 <- list()
Evers_CRISPRn_RT112$count <- 
  read_csv("https://raw.githubusercontent.com/hyunhwaj/CB2-Experiments/master/01_gene-level-analysis/data/Evers/CRISPRn-RT112.csv") %>% 
  select(-gene) %>% 
  column_to_rownames("sgRNA")

Evers_CRISPRn_RT112$egenes <- scan("https://raw.githubusercontent.com/hyunhwaj/CB2-Experiments/master/01_gene-level-analysis/data/Evers/essential-genes.txt", what="character")
Evers_CRISPRn_RT112$ngenes <- scan("https://raw.githubusercontent.com/hyunhwaj/CB2-Experiments/master/01_gene-level-analysis/data/Evers/non-essential-genes.txt", what="character")

Evers_CRISPRn_RT112$design <- 
  tribble(
    ~group, ~sample_name,
    "before", "B1",
    "before", "B2",
    "before", "B3",
    "after", "A1",
    "after", "A2",
    "after", "A3"
  )

Evers_CRISPRn_RT112$sg_stat <- run_estimation(Evers_CRISPRn_RT112$count, Evers_CRISPRn_RT112$design, "before", "after")
Evers_CRISPRn_RT112$gene_stat <- measure_gene_stats(Evers_CRISPRn_RT112$sg_stat)
usethis::use_data(Evers_CRISPRn_RT112, overwrite = T)
