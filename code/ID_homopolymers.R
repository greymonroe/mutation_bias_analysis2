library(polymorphology)
## create database of homopolymer sequences in the TAIR10 genome
genome<-read.fasta("data/TAIR10_chr_all.fas.gz")
names(genome)<-gsub("Chr","",names(genome))
TAIR10_homopolymers<-make_homopolymer(genome, 3)
TAIR10_homopolymers[]<-lapply(TAIR10_homopolymers, function(x) x[length>=3])
TAIR10_homopolymers<-rbindlist(lapply(names(TAIR10_homopolymers), function(x) {
  out<-data.table(TAIR10_homopolymers[[x]])
  out$CHROM<-gsub("Chr","",x)
  return(out)}))
#to identify variants at the end or beginning of homopolymers:
TAIR10_homopolymers$ID<-1:nrow(TAIR10_homopolymers)
TAIR10_homopolymers<-split(TAIR10_homopolymers, by="CHROM")

mutations<-fread("data/Monroe_all_mutations.csv")
mutations$chr<-as.character(mutations$CHROM)
mutations$start<-mutations$POS
mutations$stop<-mutations$POS
mutations$unique<-paste(mutations$chr, mutations$POS, mutations$ALT)

mutations_hp<-homopolymer_var_annotate(vars = mutations, homopolymers = TAIR10_homopolymers, size=3, dist=1)
mutations$potential_homopolymer<-mutations_hp

