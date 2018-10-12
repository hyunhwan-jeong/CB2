rm(list=ls())
FASTA <- system.file("extdata",
                     "small_sample.fasta",
                     package = "CB2")
df_design <- data.frame()
for(g in c("Low", "High", "Base")) {
  for(i in 1:2) {
    FASTQ <- system.file("extdata", 
                         sprintf("%s%d.fastq", g, i), 
                         package = "CB2")
    df_design <- rbind(df_design, data.frame(group = g, sample_name = sprintf("%s%d", g, i),fastq_path = FASTQ, stringsAsFactors = F))
  }
}

sgrna_count <- run_sgrna_quant(FASTA, df_design)
sgrna_stat <- run_estimation(sgrna_count, df_design, "Base", "Low")
gene_stat <- measure_gene_stats(sgrna_stat)

