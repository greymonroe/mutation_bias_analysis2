source("code/T2Tload.R")

pdf("figures/T2T_analysis.pdf", width=1.75, height=2.5)
  mutations<-germline_merge[AD_DP==1 & q20_depth>10 & counts<10 & avg_qual.a>30 & q20_depth<100]
  QCmuts(mutations,
  "Called homozygous in 1 sample
  #ALT=1
  ALT=100%
  ALT Qual>30
  10<Q20 Depth<100
  ALT Obs<10
  dist>10", windows, fill=T,file="data/germline.csv")

  mutations<-all[alt2=="" & counts<10 & counts2==1 & q20_depth>10 & strands.a==2 & plus_reads.a>1 & minus_reads.a>1 & AD_DP<.9 & avg_qual.a>30& q20_depth<100]
  QCmuts(mutations,
  "Read support F&R in 1 sample
  #ALT=1
  F > 1
  R > 1
  ALT<90%
  ALT Qual>30
  10<Q20 Depth<100
  ALT Obs<10
  dist>10", windows, fill=T,file="data/somatic2.csv")
  
  mutations<-all[alt2=="" & counts<10 & counts2==1 & q20_depth>10 & strands.a==2 & plus_reads.a>0 & minus_reads.a>0 & AD_DP<.9 & avg_qual.a>30& q20_depth<100 ]
  QCmuts(mutations,
         "Read support on F&R in 1 sample
  #ALT=1
  F > 0
  R > 0
  ALT<90%
  ALT Qual>30
  10<Q20 Depth<100
  Alt Obs<10
  dist>10", windows, fill=T, file="data/somatic.csv")
 
dev.off()

pdf("figures/T2T_analysis_ALTpct.pdf", width=1.75, height=1.5)
  mutations<-germline_merge[AD_DP==1 & q20_depth>10 & counts<10& avg_qual.a>30 & q20_depth<100]
  ALTfreqplot(mutations)
  
  mutations<-all[alt2=="" & counts<10 & counts2==1 & q20_depth>10 & strands.a==2 & plus_reads.a>1 & minus_reads.a>1 & AD_DP<.9 & avg_qual.a>30 & q20_depth<100]
  ALTfreqplot(mutations)
  
  mutations<-all[alt2=="" & counts<10 & counts2==1 & q20_depth>10 & strands.a==2 & plus_reads.a>0 & minus_reads.a>0 & AD_DP<.9 & avg_qual.a>30 & q20_depth<100]
  ALTfreqplot(mutations)
  
dev.off()


pdf("figures/T2T_analysis_noHP.pdf", width=1.75, height=2.5)
mutations<-germline_merge[AD_DP==1 & q20_depth>10 & counts<10& avg_qual.a>30 & potential_homopolymer==F& q20_depth<100]
QCmuts(mutations,
       "Called homozygous in 1 sample
  #ALT=1
  ALT=100%
  ALT Qual>30
  10<Q20 Depth<100
  ALT Obs<10", windows)

mutations<-all[alt2=="" & counts<10 & counts2==1 & q20_depth>10 & strands.a==2 & plus_reads.a>1 & minus_reads.a>1 & AD_DP<.9 & avg_qual.a>30 & potential_homopolymer==F& q20_depth<100]
QCmuts(mutations,
       "Read support F&R in 1 sample
  #ALT=1
  F > 1
  R > 1
  ALT<90%
  ALT Qual>30
  10<Q20 Depth<100
  ALT Obs<10", windows)

mutations<-all[alt2=="" & counts<10 & counts2==1 & q20_depth>10 & strands.a==2 & plus_reads.a>0 & minus_reads.a>0 & AD_DP<.9 & avg_qual.a>30 & potential_homopolymer==F& q20_depth<100]
QCmuts(mutations,
       "Read support on F&R in 1 sample
  #ALT=1
  F > 0
  R > 0
  ALT<90%
  ALT Qual>30
  10<Q20 Depth<100
  Alt Obs<10", windows)

dev.off()















