library(polymorphology)
mutations<-fread("data/germline_mutations_meta.csv")
features_data<-fread("data/features_data.csv")

pdf("figures/germline_only_constraint.pdf", width=1.5, height=2)
lapply(c("g1001_pct","PnPs","DnDs","NI","Pn_length","Ps_length","Ds_length","Dn_length","exp_cv","alpha","LOF","H3K4me1","log_exp_env_variance","log_exp_gen_variance"), function(var){
  pctile(features_data, c("CDS","intron","five_prime_UTR", "three_prime_UTR"), var, x_name=var, title="Gene body")
})
pctile(features_data, types=c("CDS","intron","five_prime_UTR", "three_prime_UTR","gene"), variable="tissue_breadth", x_name="# tissues expressed", title="Gene body")
pctile(features_data, c("CDS","intron","five_prime_UTR", "three_prime_UTR","gene"), "phen", char=T, x_name="Lethal", title="Gene body")
pctile(features_data, c("CDS","intron","five_prime_UTR", "three_prime_UTR","gene"), "isESN", char=T, x_name="Essential", title="Gene body")
dev.off()

pdf("figures/germline_only_constraint_noncoding.pdf", width=1.5, height=2)
list(lapply(c("g1001_pct","PnPs","DnDs","NI","Pn_length","Ps_length","Ds_length","Dn_length","exp_cv","alpha","LOF","H3K4me1","log_exp_env_variance","log_exp_gen_variance"), function(var){
  pctile(features_data, c("intron","five_prime_UTR", "three_prime_UTR"), var, x_name=var, title="non-coding gene body")
}))
pctile(features_data, types=c("intron","five_prime_UTR", "three_prime_UTR"), variable="tissue_breadth", x_name="# tissues expressed", title="non-coding gene body")
pctile(features_data, c("intron","five_prime_UTR", "three_prime_UTR"), "phen", char=T, x_name="Lethal", title="non-coding gene body")
pctile(features_data, c("intron","five_prime_UTR", "three_prime_UTR"), "isESN", char=T, x_name="Essential", title="non-coding gene body")
dev.off()

pdf("figures/germline_only_constraint_noncoding_nonzero.pdf", width=1.5, height=2)
list(lapply(c("g1001_pct","PnPs","DnDs","NI","Pn_length","Ps_length","Ds_length","Dn_length","exp_cv","alpha","LOF","H3K4me1","log_exp_env_variance","log_exp_gen_variance"), function(var){
  pctile(features_data[mutations>0], c("intron","five_prime_UTR", "three_prime_UTR"), var, x_name=var, title="non-coding gene body")
}))
pctile(features_data[mutations>0], types=c("intron","five_prime_UTR", "three_prime_UTR"), variable="tissue_breadth", x_name="# tissues expressed", title="non-coding gene body")
pctile(features_data[mutations>0], c("intron","five_prime_UTR", "three_prime_UTR"), "phen", char=T, x_name="Lethal", title="non-coding gene body")
pctile(features_data[mutations>0], c("intron","five_prime_UTR", "three_prime_UTR"), "isESN", char=T, x_name="Essential", title="non-coding gene body")
dev.off()


