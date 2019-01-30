library(CB2)

sgstat <- run_estimation(Evers_CRISPRn_RT112$count, Evers_CRISPRn_RT112$design, "before", "after")
gestat <- measure_gene_stats(sgstat)

plot(gestat$fdr_pa)



sgstat <- run_estimation(Sanson_CRISPRn_A375$count, Sanson_CRISPRn_A375$design, "ctl", "trt")
gestat <- measure_gene_stats(sgstat)

plot(gestat$fdr_pa)

df_i <- gestat %>% 
  mutate(essential = gene %in% Sanson_CRISPRn_A375$egenes) %>% 
  mutate(non_essential = gene %in% Sanson_CRISPRn_A375$ngenes) %>% 
  filter(essential | non_essential)
df_i %>% ggplot(aes(-log10(fdr_pa))) +
  geom_density(aes(fill=essential, color=essential), alpha=0.5) 

sum(df_i$fdr_pa[df_i$essential]<0.01) / sum(df_i$essential)
sum(df_i$fdr_pa[df_i$non_essential]<0.01) / sum(df_i$non_essential)

library(precrec)

eval <- mmdata(join_scores(1-df_i$fdr_pa), join_labels(df_i$essential)) %>% evalmod() 

eval

eval %>% autoplot()


