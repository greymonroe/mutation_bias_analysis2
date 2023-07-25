#load variants and format for analysis
library(polymorphology)
variants<-fread("~/repos/mutation_bias_analysis2/data/leaves_variants.csv")
variants$CHROM<-as.numeric(variants$CHROM)
variants$start<-as.numeric(variants$POS)
variants$end<-as.numeric(variants$POS)
setkey(variants, CHROM, start, end)
## load genes
genes<-fread("~/repos/mutation_bias_analysis2/data/A_thal_genes.csv")
genes$CHROM<-as.numeric(genes$CHROM)
setkey(genes, CHROM, start, end)
windows <- make_feature_windows(genes, 
                                2, gene = T)
windows$CHROM<-as.numeric(windows$chr)
setkey(windows, CHROM, start, stop)



##### Comparing singletons to variants with multiple read support
t.test(variants$QD~variants$Alt_depth==1)
t.test(variants$BaseQRankSum<0~variants$Alt_depth==1)
t.test(variants$MQ<60~variants$Alt_depth==1)
t.test(variants$MQRankSum<0~variants$Alt_depth==1)
t.test(variants$potential_homopolymer~variants$Alt_depth==1)

# compare genes
overlap<-foverlaps(genes, variants,type="any")

sum<-overlap[essentiality!="" , 
             .(N=sum(!is.na(POS)), 
               MQ=mean(MQ, na.rm=T),
               QD=mean(QD, na.rm=T)), 
             by=.(gene, essentiality)]


sum$length<-genes$length[match(sum$gene, genes$gene)]

sum2<-sum[essentiality!="", 
             .(N=sum(N), 
               MQ=mean(MQ, na.rm=T),
               QD=mean(QD, na.rm=T), 
               length=sum(length)),
             by=.(grp=essentiality)]


pdf("~/repos/mutation_bias_analysis2/figures/genes_essentiality_leaves.pdf", width=1, height=1.7)

ggplot(sum2, aes(x=grp, y=N/length))+
  geom_bar(stat="identity", position="dodge", fill="dodgerblue",col="black",lwd=0.5)+
  scale_y_continuous(name="Total variants/b.p.")+
  theme_classic(base_size = 6)+
  scale_x_discrete(name="Gene function")+
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggplot(sum2, aes(x=grp, y=MQ))+
  geom_bar(stat="identity", position="dodge", fill="brown1",col="black",lwd=0.5)+
  scale_y_continuous(name="Mean MQ")+
  theme_classic(base_size = 6)+
  scale_x_discrete(name="Gene function")+
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggplot(sum2, aes(x=grp, y=QD))+
  geom_bar(stat="identity", position="dodge", fill="brown1",col="black",lwd=0.5)+
  scale_y_continuous(name="Mean QD")+
  theme_classic(base_size = 6)+
  scale_x_discrete(name="Gene function")+
  theme(axis.text.x = element_text(angle=45, hjust=1))

dev.off()


windows_vars <- foverlaps(windows, variants, type = "any")
windows_means <- windows_vars[, 
                              .(N = sum(!is.na(SRR)), 
                                length = mean(length),
                                QD=mean(QD, na.rm=T),
                                MQ=mean(MQ, na.rm=T)), 
                                by = .(pos,  region)]


pdf("~/repos/mutation_bias_analysis2/figures/genes_bodies_leaves.pdf", width=2.2, height=1.7)

ggplot(windows_means, aes(x=pos, y=N/(length*nrow(genes))))+
  geom_bar(stat="identity", position="dodge", fill="dodgerblue",col="black",lwd=1)+
  geom_vline(xintercept = c(2.5,4.5))+
  scale_y_continuous(name="Total variants/b.p.")+
  theme_classic(base_size = 6)+
  scale_x_continuous(breaks=c(1,2,3.5,5,6), labels=c("-2kb","-1kb","Gene bodies","+1kb","+2kb"), name="Region relative to gene bodies")

ggplot(windows_means, aes(x=pos, y=MQ))+
  geom_bar(stat="identity", position="dodge", fill="brown1",col="black",lwd=1)+
  #geom_vline(xintercept = c(21,40))+
  scale_y_continuous(name="Average MQ")+
  theme_classic(base_size = 6)+
  scale_x_continuous(breaks=c(1,2,3.5,5,6), labels=c("-2kb","-1kb","Gene bodies","+1kb","+2kb"), name="Region relative to gene bodies")

ggplot(windows_means, aes(x=pos, y=QD))+
  geom_bar(stat="identity", position="dodge", fill="brown1",col="black",lwd=1)+
  #geom_vline(xintercept = c(21,40))+
  scale_y_continuous(name="Average QD")+
  theme_classic(base_size = 6)+
  scale_x_continuous(breaks=c(1,2,3.5,5,6), labels=c("-2kb","-1kb","Gene bodies","+1kb","+2kb"), name="Region relative to gene bodies")

dev.off()



