library(polymorphology)
mutations<-fread("data/germline_mutations_meta.csv")
features_data<-fread("data/features_data.csv")

out<-lapply(c("g1001_pct","PnPs","DnDs","NI","Pn_length","Ps_length","Ds_length","Dn_length","exp_cv","alpha","LOF","H3K4me1","log_exp_env_variance","log_exp_gen_variance"), function(var){
  pctile(unique(features_data), c("CDS","intron","five_prime_UTR", "three_prime_UTR"), var, x_name=var, title="Gene body")
})
out2<-list()
out2<-append(out2, list(pctile(features_data, types=c("CDS","intron","five_prime_UTR", "three_prime_UTR"), char=T, variable="alltissues", x_name="Expressed all tissues", title="Gene body")))
out2<-append(out2, list(pctile(features_data, c("CDS","intron","five_prime_UTR", "three_prime_UTR"), "phen", char=T, x_name="Lethal", title="Gene body")))
out2<-append(out2, list(pctile(features_data, c("CDS","intron","five_prime_UTR", "three_prime_UTR"), "isESN", char=T, x_name="Essential", title="Gene body")))
out2<-append(out2, list(pctile(features_data, c("CDS","intron","five_prime_UTR", "three_prime_UTR"), "essentiality", char=T, x_name="Function", title="Gene body")))


genes<-fread("data/A_thal_genes.csv")
genes$ID<-1:nrow(genes)
genes$CHROM<-as.character(genes$CHROM)
setkey(genes, CHROM, start, end)

mutations$start<-mutations$POS
mutations$stop<-mutations$POS
mutations$ID<-1:nrow(mutations)
mutations$CHROM<-as.character(mutations$CHROM)

setkey(mutations, CHROM, start, stop)
out_mut<-plot_peaks(genes, "Genes","Gene body", 30,mutations, gene=T)
plot(out_mut[[1]])

