dd=read.delim("../U133A_NFvsWD/u133A_WDrs_vs_NF_Pfilter_ALL_v3.txt")
data(miRNATargets)

data(miRNATargets)
ii=which(miRNATargets$gene.miRNA.map$mature.miRNA %in% miRNATargets$chr12.miRNA)
chr12GenesTargets=unique(miRNATargets$gene.miRNA.map$GENE[ii])

dd$chr12.miRNA.Target=""
dd$chr12.miRNA.Target.Consistent=0
consistentTargeting=dd$FC<0 & dd$SYMBOL %in% chr12GenesTargets
dd[consistentTargeting,]$chr12.miRNA.Target="X"
dd[consistentTargeting,]$chr12.miRNA.Target.Consistent=dd[consistentTargeting,]$FC

ans=data.frame(
    SYMBOL=dd$SYMBOL,
    chr12.miRNA.Target=dd$chr12.miRNA.Target,
    chr12.miRNA.Target.Consistent=dd$chr12.miRNA.Target.Consistent
    )

write.xls(ans,"chr12_miRNA_ConsistentTargets.txt",row.names=F)
