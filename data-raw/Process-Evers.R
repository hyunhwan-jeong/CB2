library(tidyverse)
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

usethis::use_data(Evers_CRISPRn_RT112)
