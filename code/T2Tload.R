library(polymorphology)
library(Hmisc)
source("code/funcs.R")

#reference genome
genome<-read.fasta("data/Col-CEN_v1.2.fasta.gz")

#gene annotations

gff<-fread("data/T2T-col.handerson_team.gff3")
colnames(gff)[1:7]<-c("CHROM","SRC","TYPE","START","STOP","X","direction")
genes<-gff[TYPE=="gene"]
genes$length<-genes$STOP-genes$START
genes$stop<-genes$STOP
genes$start<-genes$START
genes$chr<-genes$CHROM
genes$ID<-1:nrow(genes)
genes$gene<-gsub("ID=(.+);Note.+","\\1", genes$V9)
genes$note<-gsub(".+Note=(.+);Name.+","\\1", genes$V9)
genes<-genes[note=="protein_coding_gene"]
genes<-genes[CHROM %in% paste0("Chr",1:5)]
contents<-genecontents(genes)
genes<-merge(genes, contents)

genes_3k<-genes
genes_3k$start=genes_3k$start-3000
genes_3k$stop=genes_3k$stop+3000

exons<-gff[CHROM %in% paste0("Chr",1:5) & TYPE=="exon"]
exons$gene<-gsub("Parent=(.+)\\..+;extra.+","\\1", exons$V9)
exons<-exons[gene %in% genes$gene]


windows<-make_feature_windows(encode_data = genes, 
                              2, gene = T)
setkey(windows, chr, start, stop)

#Homopolymers
load("data/T2T_homopolymers.Rda")
T2T_homopolymers<-rbindlist(lapply(names(homopolymers), function(x) {
  out<-data.table(homopolymers[[x]])
  out$CHROM<-gsub("Chr","",x)
  return(out)}))
T2T_homopolymers$teststart<-T2T_homopolymers$start-1
T2T_homopolymers$testend<-T2T_homopolymers$end+1
T2T_homopolymers$POSITION<-T2T_homopolymers$start
T2T_homopolymers$ID<-1:nrow(T2T_homopolymers)
T2T_homopolymers$chr<-paste("Chr",T2T_homopolymers$CHROM, sep="")
setkey(T2T_homopolymers, "chr","teststart","testend")


#VCF data
vcf<-read.vcfR("data/MA/Col-0.merge.vcf.gz")

vcffix<-data.table(vcf@fix)
vcfgt<-data.table(vcf@gt)

vcfmerge<-cbind(vcffix, vcfgt)
vcfmerge$POS<-as.numeric(vcfmerge$POS)

vcfmerge$AC<-as.numeric(gsub("AC=(.+);AF.+","\\1", vcfmerge$INFO))
vcfmerge$AN<-as.numeric(gsub(".+AN=(.+);DP.+","\\1", vcfmerge$INFO))
vcfmerge$homozygous<-apply(vcfmerge, 1, function(x) sum(grepl("1/1", x)))
vcfmerge$phase<-apply(vcfmerge, 1, function(x) sum(grepl("1\\|1", x)))
vcfmerge$heterozygous<-apply(vcfmerge, 1, function(x) sum(grepl("0/1", x)))
vcfmerge$whichhomozygous<-apply(vcfmerge, 1, function(x) paste(grep("1/1|1\\|1", x), collapse = "_"))

vcfmerge$unique<-paste(vcfmerge$CHROM, vcfmerge$POS, vcfmerge$ALT, sep="_")
germline<-vcfmerge[AC==2 & (phase==1 | homozygous==1) & AN>107]
germline$call<-sapply(1:nrow(germline), function(x)  unlist(germline[x,as.numeric(unlist(germline$whichhomozygous[x])), with=F]))
germline$sample<-apply(germline, 1, function(x) colnames(germline)[as.numeric(x["whichhomozygous"])])



#Reads data
reads<-rbindlist(lapply(list.files("data/MA/", full.names = T, pattern="readcount"), function(f){
  data<-fread(cmd=paste("cat ~/repos/mutation_bias_analysis2/data/MA/header.txt",f), fill=T)
  data<-data[-1,]
  data$src<-f
  return(data)
}))

reads<-reads[depth!=0 & alt1!=""]

refs<-data.table(t(matrix(unlist(strsplit(reads$`base:reads:strands:avg_qual:map_qual:plus_reads:minus_reads`,":")), nrow=8)))
colnames(refs)<-c("base","reads","strands","avg_qual","map_qual","plus_reads","minus_reads","extra")
refs[]<-lapply(refs, type.convert)

alts<-data.table(t(matrix(unlist(strsplit(reads$alt1,":")), nrow=7)))
colnames(alts)<-c("base.a","reads.a","strands.a","avg_qual.a","map_qual.a","plus_reads.a","minus_reads.a")
alts[]<-lapply(alts, type.convert)

all_alts<-reads[, c(1,2,grep("alt", colnames(reads))), with=F]

all_alts<-data.table(apply(all_alts, 2, function(x) {
  
  gsub(":.+","",x)
  
}))

all_alts$alts<-apply(all_alts, 1, function(x) sum(!is.na(x) & x!="", na.rm=T))

uniques<-unlist(apply(all_alts, 1, function(x) {
  
  sapply(3:as.numeric(x["alts"]), function(i){
    paste(x["chrom"],x["position"],x[i], sep="_")
  })
}))
uniques_counts<-data.table(table(uniques))


all<-cbind(reads, refs)
all<-cbind(all, alts)

all$unique<-paste(all$chrom, all$position, all$base.a, sep="_")
all$counts<-uniques_counts$N[match(all$unique, uniques_counts$uniques)]

all$unique2<-paste(all$chrom, all$position, all$base.a, all$strands.a, sep="_")
counts2<-data.table(table(all$unique2))
all$counts2<-counts2$N[match(all$unique2, counts2$V1)]

all$unique3<-paste(all$chrom, all$position, sep="_")
counts3<-data.table(table(all$unique3))
all$counts3<-counts3$N[match(all$unique3, counts3$V1)]

all$REF=all$ref_base
all$ALT=all$base.a

all$depth<-as.numeric(all$depth)
all$q20_depth<-as.numeric(all$q20_depth)
all$AD_DP<-all$reads.a/(as.numeric(all$reads.a+all$reads))
all$ID<-1:nrow(all)
all$POSITION<-as.numeric(all$position)
all$POS<-as.numeric(all$position)

all$POSITION2<-as.numeric(all$position)
all$CHROM<-all$chrom
all<-all[CHROM %in% paste0("Chr",1:5)]

all$TYPE<-ifelse(nchar(all$base)==1 & nchar(all$base.a)==1, "SNP","INDEL")
all$sample<-gsub(".+//(.+).sorted.markdup.bwa.bam.readcount", "\\1", all$src)
all_simp<-all[,c("CHROM","POSITION","POSITION2","base.a","ID")]
setkey(all_simp, "CHROM", "POSITION", "POSITION2")
overlap<-foverlaps(T2T_homopolymers, all_simp,type="any")
marked<-overlap[base.a==var]
all$potential_homopolymer<-all$ID %in% marked$ID & all$TYPE=="SNP"

germline_merge<-merge(germline, all, by=c("unique","sample","CHROM","POS"))



