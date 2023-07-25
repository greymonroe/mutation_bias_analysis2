library(polymorphology)
#load mutations
mutations<-fread("data/Monroe_all_mutations.csv")
mutations$start<-mutations$POS
mutations$end<-mutations$POS
mutations$ID<-1:nrow(mutations)
mutations$CHROM<-as.character(mutations$CHROM)

setkey(mutations, CHROM, start, end)

#load genes
genes<-fread("data/A_thal_genes.csv")
setkey(genes, CHROM, start, end)
genes$muts<-mutations_in_features(genes, mutations)
genes$homopolymer_mutations<-mutations_in_features(genes, mutations[potential_homopolymer==T])


out_mut<-polymorphology::plot_peaks(encode_data = genes, "Genes","Gene body",deciles = 2,var_data = mutations, gene=T)
out_mut_hp<-polymorphology::plot_peaks(encode_data = genes, "Genes","Gene body",deciles = 2,var_data = mutations[potential_homopolymer==T], gene=T)
out_mut[[2]]$hp<-out_mut_hp[[2]]$mut

pdf("figures/HP_gene_body.pdf", width=1, height=1.7)
sums<-out_mut[[2]][,.(pcthp=sum(hp)/sum(mut)),by=.(pos,region)]
ggplot(sums, aes(x=pos, y=pcthp))+
  geom_vline(xintercept = c(2.5,4.5), linetype="dashed", size=0.25)+
  geom_bar(stat="identity", position="dodge", col="black",fill="brown1")+
  scale_y_continuous(name="proportion of variants that are potential\nhomopolymer bleed-through errors")+
  theme_classic(base_size = 6)+
  scale_x_continuous(breaks=c(1,2,3.5,5,6), labels=c("-2kb","-1kb","Gene\nbodies","+1kb","+2kb"), name="Region relative to Gene\nbodies")

sums<-out_mut[[2]][,.(pcthp=sum(hp)/sum(mut), mutpct=sum(mut)/sum(length)),by=.(pos,region)]
sums$propmut<-prop.table(sums$mutpct)
ggplot(sums, aes(x=pos, y=propmut))+
  geom_vline(xintercept = c(2.5,4.5), linetype="dashed", size=0.25)+
  geom_bar(stat="identity", position="dodge", col="black",fill="brown1")+
  scale_y_continuous(name="proportion of variants that are potential\nhomopolymer bleed-through errors")+
  theme_classic(base_size = 6)+
  scale_x_continuous(breaks=c(1,2,3.5,5,6), labels=c("-2kb","-1kb","Gene\nbodies","+1kb","+2kb"), name="Region relative to Gene\nbodies")

dev.off()




pdf("figures/HP_genes_essentiality.pdf", width=1.2, height=1.7)

sum<-genes[essentiality!="",.( muts=sum(all_mutations), hp_muts=sum(homopolymer_mutations)), by=.(grp=essentiality)]

sum$hpmutpct<-sum$hp_muts/sum$muts
ggplot(sum, aes(x=grp, y=hpmutpct))+
  geom_bar(stat="identity", position="dodge", fill="brown1",col="black")+
  scale_y_continuous(name="proportion of variants that are potential\nhomopolymer bleed-through errors")+
  theme_classic(base_size = 6)+
  scale_x_discrete(name="Gene Function")

sum<-genes[essentiality!="",.( muts=sum(muts), hp_muts=sum(homopolymer_mutations), length=sum(stop-start),homopolymers=sum(homopolymers)/sum(stop-start)), by=.(grp=essentiality)]

sum$mutpct<-sum$muts/sum$length
ggplot(sum, aes(x=grp, y=mutpct))+
  geom_bar(stat="identity", position="dodge", fill="brown1",col="black")+
  scale_y_continuous(name="proportion of variants that are potential\nhomopolymer bleed-through errors")+
  theme_classic(base_size = 6)+
  scale_x_discrete(name="Gene Function")

ggplot(sum, aes(x=grp, y=homopolymers))+
  geom_bar(stat="identity", position="dodge", fill="brown1",col="black")+
  scale_y_continuous(name="Homopolymers/bp")+
  theme_classic(base_size = 6)+
  scale_x_discrete(name="Gene Function")

dev.off()



