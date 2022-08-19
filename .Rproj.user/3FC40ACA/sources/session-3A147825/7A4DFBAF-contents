library(data.table)
library(vcfR)
library(seqinr)
library(polymorphology)


genome<-read.fasta("~/Downloads/Col-CEN_v1.2.fasta")
annotation<-fread("~/Dropbox/Research/MSH6_tudor/data/T2T-col.handerson_team.gff3")
colnames(annotation)<-c("chr", "source", "type", "start", 
                        "stop", "V1", "direction", "V2", "info")
MA_vcf<-read.vcfR("~/Dropbox/Research/MSH6_tudor/data/T2T_ngm.SNP.vcf.recode.vcf")
#MA_vcf<-read.vcfR("~/Downloads/T2T_bwa.SNP.vcf.recode.vcf")

strandbias<-fread("~/Downloads/T2T_ngm.SNP.one.strand.txt")

fix<-data.table(MA_vcf@fix)
fix$ID<-1:nrow(fix)
gt<-data.table(MA_vcf@gt)
data<-cbind(fix, gt)
data$TYPE<-ifelse(nchar(data$REF)==1 & nchar(data$ALT)==1, "SNP","InDel")

genos<-colnames(gt)[-1]
parsed<-rbindlist(apply(data[!grepl(",", ALT) & TYPE=="SNP"], 1, function(row){
  alts<-unlist(strsplit(row["ALT"], split=","))
  CHROM=row["CHROM"]
  POS=as.numeric(row["POS"])
  REF=row["REF"]
  ID=row["ID"]
  TYPE=row["TYPE"]
  NUM_ALTS<-length(alts)
  alt_annotation<-rbindlist(lapply(1:NUM_ALTS, function(a){
    ALT=alts[a]
    RDs<-as.numeric(sapply(row[genos], function(x) unlist(strsplit(unlist(strsplit(x, split=":"))[2], split=","))[1]))
    ADs<-as.numeric(sapply(row[genos], function(x) unlist(strsplit(unlist(strsplit(x, split=":"))[2], split=","))[a+1]))
    OBS<-sum(ADs>0)
    GENOS<- paste(genos[ADs>0], collapse=" ")
    CALLS<-paste(row[genos[ADs>0]], collapse=" ")
    AD<-sum(ADs)
    AD1<-sum(ADs[ADs>0])
    ALT_RD<-sum(RDs[which((ADs>0))])
    RD<-sum(RDs)
    RD0<-sum(RDs==0)
    INFO<-row["INFO"]
    BaseQRankSum<-as.numeric(unlist(strsplit(gsub(".+BaseQRankSum=(.+);DP.+","\\1",INFO), split=","))[a])
    SOR<-as.numeric(unlist(strsplit(gsub(".+SOR=(.+)","\\1",INFO), split=","))[a])
    FS<-as.numeric(unlist(strsplit(gsub(".+FS=(.+);InbreedingCoeff.+","\\1",INFO), split=","))[a])
    MQM<-as.numeric(unlist(strsplit(gsub(".+MQ=(.+);MQRank.+","\\1",INFO), split=","))[a])
    MQMR<-as.numeric(gsub(".+MQRankSum=(.+);QD.+","\\1",INFO))
    QD<-as.numeric(gsub(".+QD=(.+);ReadPosRankSum.+","\\1",INFO))
    
    return(data.table(CHROM, POS, ID, REF, ALT, INFO, TYPE, CALLS, OBS, NUM_ALTS, AD, AD1, RD, ALT_RD,GENOS, BaseQRankSum, QD, FS, SOR, RD0, MQM, MQMR))
  }))
  return(alt_annotation)
}))

singles<-parsed[OBS==1]
singles$unique<-paste(singles$CHROM, singles$POS, sep="-")

genes<-annotation[type=="gene"]

singles$genic<-apply(singles, 1, function(x){
  CHROM<-x["CHROM"]
  POS<-as.numeric(x["POS"])
  sum(nrow(genes[chr==CHROM & POS>=start-3000 & POS<=stop+3000]))
})

snp_split<-lapply(genos, function(geno) split(singles[grepl(geno, GENOS)],  by=c("CHROM")))
names(snp_split)<-genos

singles$dist<-unlist(lapply(1:nrow(singles), function(i){
  x<-unlist(singles[i])
  chr<-x["CHROM"]
  POS<-as.numeric(x["POS"])
  SRR<-x["GENOS"]
  snp_split
  dist<-snp_split[[SRR]][[chr]]$POS-POS
  dist<-min(abs(dist[dist!=0]))
  return(dist)
}))

#load("data/T2T_homopolymers.Rda")

singles_homopolymer<-homopolymer_context(vars = singles, genome, homopolymers, size = 3)
save(singles_homopolymer, file="data/singles_homopolymer.Rda")
#load("data/singles_homopolymer.Rda")

singles$nexttohomopolymer<-singles_homopolymer$nexttohomopolymer

final<-singles[dist>10 & genic>0 & CHROM %in% paste0("Chr",1:5) &  !unique %in% strandbias$SNP & nexttohomopolymer==0 & MQM==60 & AD/ALT_RD<1 & AD>1]

pdf("figures/Weng_somatic_tss_reanalsysis_T2T_ngm.pdf", width=3, height=1.5)
tss<-tss_tts.variants(gff=annotation, vcf=final) 
tss_tts.variants.plot(tss, window=200)
dev.off()




