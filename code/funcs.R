

calcdists<-function(mutations){
  dist<-apply(mutations, 1, function(x){
    S<-x["sample"]
    C<-x["CHROM"]
    P<-as.numeric(x["POS"])
    check<-mutations[CHROM==C & sample==S]
    D=abs(check$POS-P)
    if(length(D)==1){
      return(Inf)
    } else return(min(D[D!=0]))
    
  })
}

plotgenebodies<-function(variants, title="Mutations relative to gene bodies", windows){
  variants<-variants[,c("CHROM","POS")]
  variants$POS2<-variants$POS
  variants$start=variants$POS
  setkey(variants, "CHROM", "POS","POS2")
  variants<-unique(variants)
  variants$ID<-1:nrow(variants)
  n=(nrow(variants))
  
  p<-plot_peaks(genes, "Genes, SBS","Gene body", 2,variants, gene=T, window=windows)[[2]]
  sum<-p[pos %in% 1:6,.(N=sum(mut), notmut=sum(length)-sum(mut), length=sum(length)), by=.(region=="gene body")]
  print(sum[,2:3])
  chi<-chisq.test(sum[,2:3])
  
  print(chi)
  plot<-ggplot(p, aes(x=pos, y=pct, fill=region=="gene body"))+
    geom_vline(xintercept = c(2.5,4.5), linetype="dashed", size=0.25)+
    geom_bar(stat="identity", col="black")+
    #geom_line()+
    scale_fill_manual(values=c("dodgerblue4", "dodgerblue"), guide="none")+
    scale_y_continuous(name="Mutations/b.p.")+
    theme_classic(base_size = 6)+
    scale_x_continuous(breaks=c(1,2,3.5,5,6), labels=c("-2kb","-1kb","Gene bodies","+1kb","+2kb"), name="Region relative to gene bodies")+
    ggtitle(paste0("Gene bodies vs. Intergenic Chi.sq \nn=",n, " X2=", round(chi$statistic,digits = 1), " df=",chi$parameter, " p=",round(chi$p.value, digits = 8)))+
    labs(caption = title)+
    theme(title=element_text(size=5))
  
  plot(plot)
}


plot_peaks<-function (encode_data, ggtitle, xtitle, deciles, var_data, gene = F, 
                      marky = F, windows) 
{
  
  windows_vars <- foverlaps(windows, var_data, type = "any")
  windows_means <- windows_vars[, .(mut = sum(!is.na(ID)), 
                                    N = .N, length = mean(length, na.rm = T)), by = .(pos, 
                                                                                      region, i.ID)]
  windows_vars[is.na(i.ID)]
  windows_means <- windows_means[, .(pct = sum(mut)/sum(length), 
                                     mut = sum(mut), length = sum(length)), by = .(pos, region)]
  
  if (marky == T) {
    windows_means$pct <- windows_means$mut
  }
  plot <- ggplot(windows_means, aes(x = pos, y = pct, col = region == 
                                      "gene body", group = region)) + geom_vline(xintercept = c(deciles, 
                                                                                                deciles * 2), linetype = "dashed", size = 0.25) + geom_line() + 
    scale_color_manual(values = c("gray75", "green3"), guide = "none") + 
    theme_classic(base_size = 6) + scale_y_continuous(name = "Mutations/bp") + 
    ggtitle(ggtitle) + scale_x_continuous(breaks = c(0, max(windows_means$pos)/3, 
                                                     max(windows_means$pos)/3 * 2, max(windows_means$pos)), 
                                          labels = c("-2kb", "0%", "100%", "+2kb"), name = xtitle)
  return(list(plot, windows_means))
}

make_feature_windows<-function (encode_data, deciles = 10, gene = F) 
{
  library(seqinr)
  library(openxlsx)
  library(data.table)
  windows <- rbindlist(apply(encode_data, 1, function(x) {
    chr = x["chr"]
    body_starts = round(seq(as.numeric(x["start"]), as.numeric(x["stop"]), 
                            length.out = deciles + 1)[-(deciles + 1)])
    body_stops <- round(seq(as.numeric(x["start"]), as.numeric(x["stop"]), 
                            length.out = deciles + 1)[-1])
    upstream_starts <- seq(as.numeric(x["start"]) - 2000, 
                           as.numeric(x["start"]), length.out = deciles + 1)[-(deciles + 
                                                                                 1)]
    upstream_stops <- seq(as.numeric(x["start"]) - 2000, 
                          as.numeric(x["start"]), length.out = deciles + 1)[-1]
    downstream_starts <- seq(as.numeric(x["stop"]), as.numeric(x["stop"]) + 
                               2000, length.out = deciles + 1)[-(deciles + 1)]
    downstream_stops <- seq(as.numeric(x["stop"]), as.numeric(x["stop"]) + 
                              2000, length.out = deciles + 1)[-1]
    out <- data.table(chr = x["chr"], start = c(upstream_starts, 
                                                body_starts, downstream_starts), stop = c(upstream_stops, 
                                                                                          body_stops, downstream_stops), region = c(rep("upstream", 
                                                                                                                                        length(upstream_starts)), rep("gene body", length(body_starts)), 
                                                                                                                                    rep("downstream", length(downstream_starts))), ID = as.numeric(x["ID"]))
    out$pos <- 1:nrow(out)
    out$length <- out$stop - out$start
    if (gene == T) {
      direction = x["direction"]
      if (direction == "-") {
        out$pos <- rev(out$pos)
        out$region <- rev(out$region)
      }
    }
    return(out)
  }))
  setkey(windows, chr, start, stop)
  return(windows)
}


plotgenebodies_QUAL<-function(variants){
  variants<-singles
  setkey(variants, CHROM, POSITION, POSITION2)
  overlaps <- foverlaps(windows, variants[,c("CHROM","POSITION","POSITION2","avg_qual.a","avg_qual"), with=F], type = "any")
  overlaps<-overlaps[!is.na(POSITION)]
  sum<-overlaps[,.(N=.N, length=mean(length), avg_qual.a=mean(avg_qual.a/avg_qual)), by=.(region, pos)]
  P1<-ggplot(sum, aes(x=pos, y=avg_qual.a))+
    geom_bar(stat="identity")+
    scale_y_continuous(name="Variants/b.p.")+
    theme_classic(base_size = 6)
  
  plot(P1)
}


makevariants<-function(singles){
  variants<-singles
  variants<-variants[,c("CHROM","POSITION","POSITION2","base","base.a","counts")]
  setkey(variants, "CHROM", "POSITION","POSITION2")
  variants$start<-variants$POSITION
  variants<-unique(variants)
  return(variants)
}


QCmuts<-function(variants, info="Mutations", windows, breaks=10, fill=T, file=NULL){
  
  if(fill==T){
    variants$dist<-calcdists(variants)
    variants$gene_3k<-features_overlap_mutation(genes_3k, variants)
    variants$genic<-features_overlap_mutation(genes, variants)
    variants$exonic<-features_overlap_mutation(exons, variants)
    
    if(!is.null(file)){fwrite(variants, file)}
   
    plotgenebodies(variants[],info, windows)
  } else {
    plotgenebodies(variants,info, windows)
  }
}

ALTfreqplot<-function(variants){
  
  ALTQUAL=mean(variants$avg_qual.a)
  PCTHP=prop.table(table(variants[TYPE=="SNP"]$potential_homopolymer))[2]
  IE<-exonic(variants)[1]/exonic(variants)[2]
  p<-ggplot(variants, aes(x=AD_DP))+
    geom_density(fill="azure3", col="black", alpha=0.75)+
    theme_classic(base_size = 6)+
    scale_y_continuous(name="Density")+
    scale_x_continuous(name="Alt read / Alt + Ref", limits=c(0,1.05))+
    ggtitle(paste0("ALT Qual=",round(ALTQUAL, 2), "\nPotential HP=", round(PCTHP*100,2),"%" ,"\nIntron/Exon=", round(IE,2)))
  plot(p)
  
  if(all(variants$AD_DP==1)){
    p<-ggplot(variants, aes(x=AD_DP))+
      geom_histogram(fill="azure3", col="black", alpha=0.75)+
      theme_classic(base_size = 6)+
      scale_y_continuous(name="Density")+
      scale_x_continuous(name="Alt read / Alt + Ref", limits=c(0,1.05))+
      ggtitle(paste0("ALT Qual=",round(ALTQUAL, 2), "\nPotential HP=", round(PCTHP*100,2),"%" ,"\nIntron/Exon=",  round(IE,2)))
    plot(p)
    
  }
}


exonic<-function(variants){
  variants$genic<-features_overlap_mutation(genes, variants)
  variants$exonic<-features_overlap_mutation(exons, variants)
  prop.table(table(variants[genic==T]$exonic))
}

writevariants<-function(variants, file){
  
  variants$dist<-calcdists(variants)
  variants$gene_3k<-features_overlap_mutation(genes_3k, variants)
  variants$genic<-features_overlap_mutation(genes, variants)
  variants<-variants[gene_3k==T]
  write.table(variants[,c("CHROM","POS","REF","ALT","sample")], file)
  
}


strReverse <- function(x){
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
}

addmutstogenes<-function(variants, genes){

  variants$POS2<-variants$POS
  variants$CHROM<-gsub("Chr|chr","", as.character(variants$CHROM))
  genes$CHROM<-gsub("Chr|chr","", as.character(genes$CHROM))
  
  setkey(variants, CHROM, POS, POS2)
  setkey(genes, CHROM, START, STOP)
  
  
  overlap<-rbindlist(lapply(as.character(1:5), function(c){
    foverlaps(variants[CHROM ==c], genes[CHROM==c])[!is.na(direction)]
  }))
  overlap$mut2<-paste(overlap$REF, overlap$ALT, sep=">")
  overlap<-overlap[nchar(mut2)==3]
  
  sum<-data.table(table(mut=overlap$mut2, gene=overlap$gene ))
  sumdcast<-dcast(sum, gene~mut, value.var = "N")
  sumdcast$totalmut<-rowSums(sumdcast[,2:ncol(sumdcast)])
  
  merge<-merge(genes[,"gene", with=F], sumdcast, by="gene", all=T)
  setnafill(merge, cols=2:ncol(merge), fill=0)
  merge<-merge(genes, merge, by="gene", all=T)
  
  return(merge)
}

countstrandbias<-function(genes){
  
  strandbias<-genes[totalmut>0,.(N=.N,  
                AT=sum(`A>T`),
                TA=sum(`T>A`),
                CT=sum(`C>T`),
                GA=sum(`G>A`), 
                GT=sum(`G>T`),
                CA=sum(`C>A`),
                
                AC=sum(`A>C`),
                TG=sum(`T>G`),
                CG=sum(`C>G`),
                GC=sum(`G>C`), 
                TC=sum(`T>C`),
                AG=sum(`A>G`),
                
                Tsum=sum(T_s),  
                Asum=sum(A_s), 
                Csum=sum(C_s),  
                Gsum=sum(G_s), 
                Bsum=sum(C_s, G_s, T_s),
                Vsum=sum(C_s, A_s, G_s),
                Dsum=sum(G_s, A_s, T_s),
                Hsum=sum(C_s, A_s, T_s),
                length=sum(sum)), 
             by=.(direction)]
  return(strandbias)
}

countstrandbias_UT<-function(strandbias){
  levels<-c("AT","CT","GT","GC","AC","TC")

  strandbias<-strandbias[order(direction)]
  x<-levels[1]
  UT<-rbindlist(lapply(levels, function(x){
    ref<-substr(x, 1,1)
    alt<-substr(x, 2,2)
    mut<-paste0(ref, ">",alt)
    rev<-paste0(toupper(comp(ref)),toupper(comp(alt)))
    tmp<-strandbias[,c(x, rev), with=F]
    template=unlist(tmp[1,2]+tmp[2,1])
    nontemplate=unlist(tmp[1,1]+tmp[2,2])
    
    tmp<-strandbias[,c(paste0(ref,"sum"), paste0(toupper(comp(ref)),"sum")), with=F]
    
    template_ref=unlist(tmp[1,2]+tmp[2,1])
    nontemplate_ref=unlist(tmp[1,1]+tmp[2,2])
    
    chi<-chisq.test(c(template, nontemplate), p=prop.table(c(template_ref, nontemplate_ref)))
    return(data.table(mut, template=template, nontemplate=nontemplate, template_ref, nontemplate_ref, p=chi$p.value))
  }))
  colnames(UT)<-c("Mutation","Template","Non-template","Template Ref","Non-template Ref", "P")
  UT$Mutation<-factor(UT$Mutation, levels=UT$Mutation)
  return(UT)
}

UTplot<-function(UT){
  ggplot(UT, aes(x=Mutation, y=log((Template/`Template Ref`)/(`Non-template`/`Non-template Ref`)), fill=P<0.05))+
    geom_bar(stat="identity")+
    theme_classic(base_size = 6)
}


strandbiaschi<-function(strandbias, mut){
  from1<-substr(mut,1,1)
  from2<-toupper(comp(from1))
  to1<-substr(mut,2,2)
  to2<-toupper(comp(to1))
  obs<-c(unlist(strandbias[,paste0(from1, to1), with=F]), unlist(strandbias[,paste0(from2, to2), with=F]))
  exp<-prop.table(c(unlist(strandbias[,paste0(from1,"sum"), with=F]), unlist(strandbias[,paste0(from2, "sum"), with=F])))
  chi<-chisq.test(obs, p=exp)
  obs/chi$expected
}

strandbiaschi2<-function(strandbias, mut){
  from1<-substr(mut,1,1)
  from2<-toupper(comp(from1))
  to1<-substr(mut,2,2)
  to2<-toupper(comp(to1))
  obs<-(strandbias[,c(paste0(from1, to1), paste0(from2, to2)), with=F])
  chi<-chisq.test(obs)
  data.table(obs/chi$expected)
}

strandbiaschi2all<-function(strandbias){
  out<-do.call(cbind, lapply(c("AT","CT","GT","GC","AC","TC"), function(x){
    strandbiaschi2(strandbias, x)
  }))
  out$dir<-c("+","-")
  return(out)
}

plotstrandbiasOE<-function(strandbiaschi2all){
  strandbiasOEmelt<-melt(strandbiasOE, id.vars = "dir")
  p1<-ggplot(strandbiasOEmelt, aes(x=variable, y=value, col=dir, group=dir))+
    geom_point()
  list(p1)
}

plotstrandbias<-function(strandbias){
  
  strandbiaschi(strandbias, "AT")
  p1<-ggplot(strandbias, aes(x=direction, y=((AT/(Asum))/(TA/(Tsum))), fill=direction, group=direction))+
    geom_bar(stat="identity", position="dodge", col="black")+
    theme_classic()+
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    scale_y_continuous(name="log2(A>T per A /T>A per T)")
 
   strandbiaschi(strandbias, "GT")
  p2<-ggplot(strandbias, aes(x=direction, y=((GT/(Gsum))/(CA/(Csum))), fill=direction, group=direction))+
    geom_bar(stat="identity", position="dodge", col="black")+
    theme_classic()+
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    scale_y_continuous(name="log2(G>T per G / C>A per C)")

  strandbiaschi(strandbias, "CT")
  p3<-ggplot(strandbias, aes(x=direction, y=(CT/(Csum)/(GA/(Gsum))), fill=direction, group=direction))+
    geom_bar(stat="identity", position="dodge", col="black")+
    theme_classic()+
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    scale_y_continuous(name="log2(C>T per C / G>A per G)")
  
  strandbiaschi(strandbias, "AC")
  p4<-ggplot(strandbias, aes(x=direction, y=((AC/(Asum))/(TG/(Tsum))), fill=direction, group=direction))+
    geom_bar(stat="identity", position="dodge", col="black")+
    theme_classic()+
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    scale_y_continuous(name="log2(A>C per A/T>G per T)")
  
  strandbiaschi(strandbias, "GC")
  p5<-ggplot(strandbias, aes(x=direction, y=((GC/(Gsum))/(CG/(Csum))), fill=direction, group=direction))+
    geom_bar(stat="identity", position="dodge", col="black")+
    theme_classic()+
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    scale_y_continuous(name="log2(G>T per G / C>A per C)")
  
  strandbiaschi(strandbias, "TC")
  p6<-ggplot(strandbias, aes(x=direction, y=(TC/(Tsum)/(AG/(Asum))), fill=direction, group=direction))+
    geom_bar(stat="identity", position="dodge", col="black")+
    theme_classic()+
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    scale_y_continuous(name="log2(T>C per T / A>G per A)")
  
  
  list(p1, p2, p3, p4, p5, p6)
}

plotstrandbias_melt<-function(strandbias){
  
  strandbias$direction<-ifelse(strandbias$direction=="-","reverse","forward")
  strandbias$AT_A<-strandbias$AT/strandbias$Asum
  strandbias$TA_T<-strandbias$TA/strandbias$Tsum
  strandbiasdmelt<-melt(strandbias, id.vars="direction", measure.vars = c("AT_A","TA_T"))
  
  p1<-ggplot(strandbiasdmelt, aes(x=variable, y=value, fill=direction, group=direction))+
    geom_bar(stat="identity", position="dodge", col="black")+
    theme_classic(base_size = 6)+
    facet_grid(~direction)+
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    scale_y_continuous(name="Variant/b.p.")+
    scale_x_discrete(name="", labels=c("A>T","T>A"))
  
  strandbias$CT_A<-strandbias$CT/strandbias$Csum
  strandbias$GA_T<-strandbias$GA/strandbias$Gsum
  strandbiasdmelt<-melt(strandbias, id.vars="direction", measure.vars = c("CT_A","GA_T"))
  
  p2<-ggplot(strandbiasdmelt, aes(x=variable, y=value, fill=direction, group=direction))+
    geom_bar(stat="identity", position="dodge", col="black")+
    theme_classic(base_size = 6)+
    facet_grid(~direction)+
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    scale_y_continuous(name="Variant/b.p.")+
    scale_x_discrete(name="", labels=c("C>T","G>A"))
  
  strandbias$GT_A<-strandbias$GT/strandbias$Gsum
  strandbias$CA_T<-strandbias$CA/strandbias$Csum
  strandbiasdmelt<-melt(strandbias, id.vars="direction", measure.vars = c("GT_A","CA_T"))
  
  p3<-ggplot(strandbiasdmelt, aes(x=variable, y=value, fill=direction, group=direction))+
    geom_bar(stat="identity", position="dodge", col="black")+
    theme_classic(base_size = 6)+
    facet_grid(~direction)+
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    scale_y_continuous(name="Variant/b.p.")+
    scale_x_discrete(name="", labels=c("G>T","C>A"))
  list(p2)
}

genecontents<-function(genes){
  
  genes$start<-genes$START
  genes$stop<-genes$STOP
  contents<-rbindlist(apply(genes, 1, function(x) {
    chr<-genome[[x["CHROM"]]]
    counts<-data.table(table(toupper(chr[x["start"]:x["stop"]])))
    counts<-counts[V1 %in% c("A","C","G","T")]
    counts<-dcast(counts, .~V1, value.var="N")
    counts$gene=x["gene"]
    counts$direction=x["direction"]
    counts$CHROM=x["CHROM"]
    return(counts[,2:8])
  }))
  
  contents<-contents[direction %in% c("+","-")]
  contents$sum<-rowSums(contents[,1:4])
  colnames(contents)[1:4]<-paste(colnames(contents)[1:4],'s', sep="_")
  return(contents)
}
