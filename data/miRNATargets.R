chr12.miRNA=scan("db/chr12q_mirna.txt","",quiet=T)
load("db/targetsUnion.Rdata")
miRNATargets=list(chr12.miRNA=chr12.miRNA,gene.miRNA.map=gene.miRNA.map)
rm("chr12.miRNA")
rm("gene.miRNA.map")

