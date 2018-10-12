rm(list=ls())

library(CRISPRExpress)

FASTA <- "~/Downloads/Human_CRISPR_Evers.fasta"
df_design <- data.frame(
  group = c("T0", "T0", "T0", "T1", "T1", "T1"),
  sample_name = c("T0_1", "T0_2", "T0_3", "T1_1", "T1_2", "T1_3"),
  fastq_path = sort(Sys.glob("~/Downloads/CE/*.fastq")),
  stringsAsFactors = F) 
  
sgrna_count <- run_sgrna_quant(FASTA, df_design)
sgrna_stat <- run_estimation(sgrna_count, df_design, "T1", "T0")
gene_stat <- measure_gene_stats(sgrna_stat)

