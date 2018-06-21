library(apcluster)

# load gavin pathogenic variants
v <- read.table(gzfile("/Users/joeri/github/g2files/gavin1_r0.5/calibrationvariants_r0.5.withcaddannot.tsv.gz"), sep="\t", header=T)
vp <- subset(v, group == "PATHOGENIC")

# select relevant columns for clustering based on CADD annotations
# for details see: http://cadd.gs.washington.edu/static/ReleaseNotes_CADD_v1.3.pdf
# note that the index is offset by +9 relative to this document (ie. Length at 7 now 16, ConsScore at 12 now 21 etc.)
vpcol <- vp[,c(1,16,21,23:80,82:84,87:90,100,101,106:111,113,119,121,123:125)]

# calculate aggregate values for gene and set gene names as the row names
vpsagr <- aggregate(vpcol[,2:length(vpcol)], by=list(gene=vpcol$gene), FUN=mean, na.rm=TRUE)
rownames(vpsagr) <- vpsagr[,1]
vpsagr[,1] <- NULL

# perform affinity propagation clustering
genes_apclus <- apcluster(negDistMat(r=2), vpsagr, q=0.1, seed=1)
genes_apclus

# write results
sink("~/gavin1_r0.5_patho_genes_apclus_q0.1.csv")
genes_apclus
sink()
