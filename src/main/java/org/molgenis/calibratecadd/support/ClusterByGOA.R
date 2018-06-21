library(apcluster)

# load GO Annotation data
URL <- "http://geneontology.org/gene-associations/goa_human.gaf.gz"
download.file(URL, destfile = "goa_human.gaf.gz", method="curl")
goaFull <- read.table(gzfile("goa_human.gaf.gz"), sep="\t", header=F, comment.char = "!",quote="")

# load GAVIN variants and get unique genes
variants <- read.table(gzfile("/Users/joeri/github/gavin/data/other/calibrationvariants_r0.5.tsv.gz"), sep="\t", header=T)
genes <- unique(variants$gene)

# select only gene and GO columns from data, condense by unique, and subselect GAVIN genes
goa <- goaFull[,c(3,5)]
colnames(goa) <- c("gene", "goterm")
goa <- unique(goa)
goaGav <- goa[goa$gene %in% genes,]

# construct binary presence/absence matrix of GO terms per gene
pamatrix <- dcast(goaGav,gene~goterm,fun.aggregate = function(x){as.integer(length(x) > 0)})
rownames(pamatrix) <- pamatrix[,1]
pamatrix[,1] <- NULL

# cluster the genes based on GO terms and write results to file
genes_apclus <- apcluster(negDistMat(r=2), pamatrix, q=0.1, seed=1)
sink("/Users/joeri/github/g2files/gavin1_r0.5/genes_apclus_GOterms_q0.1.csv")
genes_apclus
sink()

# quit, good practice when running as screened process
q("no")

#run: nohup R CMD BATCH ClusterByGOA.R > nohup.out &
