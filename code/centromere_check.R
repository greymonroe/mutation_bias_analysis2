library(polymorphology)
#load mutations
mutations<-fread("data/Monroe_all_mutations.csv")

pdf("figures/used_mutations.pdf", width=6, height=2)
som<-mutations[src %in% c("MA_somatic")]
ggplot(som, aes(x=POS/1000000, fill=modelused))+
  geom_histogram()+
  facet_grid(modelused~CHROM, scales="free_x")+
  theme_classic(base_size = 6)+
  scale_x_continuous(name="Position (Mb)")
dev.off()
