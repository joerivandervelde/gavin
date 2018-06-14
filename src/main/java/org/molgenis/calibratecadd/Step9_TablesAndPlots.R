######################################################
# Performance benchmark result processing and plotting
######################################################

version <- "r0.5"
pathToGavinGitRepo <- "/Users/joeri/github/gavin"

library(ggplot2)
lightgray <- "#cccccc"; gray <- "#999999"; orange <- "#E69F00"; skyblue <- "#56B4E9"; blueishgreen <- "#009E73"
yellow <- "#F0E442"; blue <-"#0072B2"; vermillion <- "#D55E00"; reddishpurple <- "#CC79A7"

# Full performance results, Supplementary Table 1 in paper
source(paste(pathToGavinGitRepo,"/data/other/step9_panels_out_",version,".R",sep=""))
df$TotalToolClsfNonVOUS <- 0
df$TotalExpertClsfNonVOUS <- 0
df$Coverage <- 0
df$MCCadj <- 0
df$Sensitivity <- 0
df$Specificity <- 0
df$PPV <- 0
df$NPV <- 0
df$ACC <- 0
for(i in 1:nrow(df)) {
  df[i,]$TotalToolClsfNonVOUS <- df[i,]$TP+df[i,]$TN+df[i,]$FP+df[i,]$FN
  df[i,]$TotalExpertClsfNonVOUS <- df[i,]$TotalToolClsfNonVOUS + df[i,]$ExpBenignAsVOUS + df[i,]$ExpPathoAsVOUS
  df[i,]$Coverage <- df[i,]$TotalToolClsfNonVOUS / df[i,]$TotalExpertClsfNonVOUS
  df[i,]$MCCadj <- df[i,]$MCC * df[i,]$Coverage
  df[i,]$Sensitivity <- df[i,]$TP/(df[i,]$TP+df[i,]$FN+df[i,]$ExpPathoAsVOUS)
  df[i,]$Specificity <- df[i,]$TN/(df[i,]$TN+df[i,]$FP+df[i,]$ExpBenignAsVOUS)
  df[i,]$PPV <- df[i,]$TP/(df[i,]$TP+df[i,]$FP)
  df[i,]$NPV <- df[i,]$TN/(df[i,]$TN+df[i,]$FN)
  df[i,]$ACC <- (df[i,]$TP+df[i,]$TN)/(df[i,]$TP+df[i,]$TN+df[i,]$FP+df[i,]$FN+df[i,]$ExpPathoAsVOUS+df[i,]$ExpBenignAsVOUS)
}
df.selectedcolumns <- df[,c(1,2,5,4,6,7,9,8,16,17)]
#save: write.table(df.selectedcolumns, file = "~/gavinbenchmark.tsv",row.names=FALSE, sep="\t")

# Aggregated performance stats per tool, Table 3 in paper
stats <- data.frame()
for(i in unique(df$Tool)) {
  df.sub <- subset(df, Tool == i)
  medianSens <- median(df.sub$Sensitivity)
  medianSpec <- median(df.sub$Specificity)
  medianPPV <- median(df.sub$PPV)
  medianNPV <- median(df.sub$NPV)
  row <- data.frame(Tool = i, medianSens = medianSens, medianSpec = medianSpec, medianPPV = medianPPV, medianNPV = medianNPV)
  stats <- rbind(stats, row)
}
stats


# Figure 1 in paper
df.notincgd <- subset(df, Data == "NotInCGD")
ggplot() +
  theme_bw() + theme(panel.grid.major = element_line(colour = "black"), axis.text=element_text(size=12),  axis.title=element_text(size=14,face="bold")) +
  # cleaner look: theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=12),  axis.title=element_text(size=14,face="plain")) +
  geom_point(data = df, aes(x = Sensitivity, y = Specificity, shape = Tool, colour = Tool), size=3, stroke = 3, alpha=0.75) +
  geom_text(data = df, aes(x = Sensitivity, y = Specificity, label = Data), hjust = 0, nudge_x = 0.01, size = 3, check_overlap = TRUE) +
  geom_text(data = df.notincgd , aes(x = Sensitivity, y = Specificity, label = Data), hjust = 0, nudge_x = 0.01, size = 3, colour="red",check_overlap = TRUE) +
  scale_colour_manual(values=rainbow(16)) +
  scale_shape_manual(values = c(1, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)) +
  scale_x_continuous(lim=c(0.377,1), breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme(legend.key = element_blank())
ggsave("Figure1.pdf", width = 25, height = 32, units = "cm") #thesis: 25 x 32


###########################################################
# Miscellaneous numbers and plots, used in paper or website
###########################################################

# Percentage variants not classified across total benchmark variant count
df.ponp2 <- subset(df, Tool == "PONP2")
(sum(df.ponp2$TotalExpertClsfNonVOUS)-sum(df.ponp2$TotalToolClsfNonVOUS)) / sum(df.ponp2$TotalExpertClsfNonVOUS)
# GAVIN per panel
df.gavin <- subset(df, Tool == "GAVIN")
df.gavin.sub <- df.gavin[,c(1,2,5,4,6,7,9,8,16,17)]
#save: write.table(df.gavin.sub, file = "~/gavinsub.tsv",row.names=FALSE, sep="\t")

# Basic calculation of pathogenic MAF threshold, CADD means, etc.
# Uses calibration results and not the benchmark output
calibrations <- read.table(paste(pathToGavinGitRepo,"/data/predictions/GAVIN_calibrations_",version,".tsv",sep=""),header=TRUE,sep='\t',quote="",comment.char="",as.is=TRUE)
mean(calibrations$PathoMAFThreshold, na.rm = T)
mean(calibrations$MeanPopulationCADDScore, na.rm = T)
mean(calibrations$MeanPathogenicCADDScore, na.rm = T)
mean(calibrations$MeanDifference, na.rm = T)
sd(calibrations$MeanDifference, na.rm = T)
table(calibrations$Category)

# Some plots/stats on the variants used in calibration
variants <- read.table(gzfile(paste(pathToGavinGitRepo,"/data/other/calibrationvariants_",version,".tsv.gz",sep="")), sep="\t", header=T)
ggplot() +
  theme_bw() + theme(panel.grid.major = element_line(colour = "black"), axis.text=element_text(size=12),  axis.title=element_text(size=14,face="bold")) +
  geom_jitter(data = variants, aes(x = cadd, y = effect, colour = group, alpha=group), stroke = 2, size=2, position = position_jitter(width = .5, height=.5)) +
  scale_colour_manual(values=c(vermillion, skyblue)) +
  scale_alpha_discrete(range = c(.75, .25))

# load rule guide
rules <- read.table(paste(pathToGavinGitRepo,"/data/predictions/GAVIN_ruleguide_",version,".tsv",sep=""),header=TRUE,sep='\t',quote="",comment.char="",as.is=TRUE, skip=33)
rules$PathogenicIfCADDScoreGreaterThan <- as.numeric(rules$PathogenicIfCADDScoreGreaterThan)
rules$BenignIfCADDScoreLessThan <- as.numeric(rules$BenignIfCADDScoreLessThan)


## Gene plots, for single selected genes or all genes in a for-loop

# refseq exons used in new style plot
refeqExons <- read.table(gzfile(paste(pathToGavinGitRepo,"/data/other/refseq-exons-list.tsv.gz", sep="")), sep="\t", header=T)

for (selectGene in unique(rules$Gene)) {
  #selectGene <- "TTN"
  rules.selectedGene <- subset(rules, Gene == selectGene)
  if(startsWith(rules.selectedGene$CalibrationCategory, "C")) {
    variants.selectedGene <- subset(variants, gene == selectGene)
    variants.selectedGene.path <- subset(variants.selectedGene, group == "PATHOGENIC")
    variants.selectedGene.popul <- subset(variants.selectedGene, group == "POPULATION")
    refeqExonsGene <- subset(refeqExons, gene == selectGene)
    if(dim(refeqExonsGene)[[1]]==0){ refeqExonsGene <- data.frame(gene=selectGene,start=-1,end=-1) }
    xmin <- min(variants.selectedGene.popul$pos, variants.selectedGene.path$pos)
    xmax <- max(variants.selectedGene.popul$pos, variants.selectedGene.path$pos)
    ymax <- max(variants.selectedGene.popul$cadd, variants.selectedGene.path$cadd)
    pth <- rules.selectedGene$PathogenicIfCADDScoreGreaterThan
    bth <-rules.selectedGene$BenignIfCADDScoreLessThan
    p <- ggplot() +
      xlim(xmin,xmax) +
      geom_rect(aes(xmin = -Inf, xmax = +Inf, ymin = 15, ymax = ymax), fill = "#FFE0E1", alpha = 1) + # lighter red
      geom_rect(aes(xmin = -Inf, xmax = +Inf, ymin = 0, ymax = 15), fill = "#DFE2FD", alpha = 1) + # lighter blue
      geom_rect(aes(xmin = -Inf, xmax = +Inf, ymin = 0, ymax =  bth), fill = "#C9CEFC", alpha = 1) + # blue
      geom_rect(aes(xmin = -Inf, xmax = +Inf, ymin = pth, ymax = ymax), fill = "#FFCACC", alpha = 1) + # red
      geom_abline(intercept = 15, slope = 0, colour = "darkgrey", size = 2, alpha = 1) +
      geom_abline(intercept = bth, slope = 0, colour = "#6681F7", size = 2, alpha = 1) + # dark blue
      geom_abline(intercept = pth, slope = 0, colour = "#FF6E77", size = 2, alpha = 1) + # dark red
      geom_rect(data = refeqExonsGene, aes(xmin = start, xmax = end, ymin = 0, ymax = ymax), linetype = 0, fill="green3", alpha = 1) +
      geom_point(data = variants.selectedGene.popul, aes(x=pos, y=cadd), colour="blue", pch=19, alpha = 1) +
      geom_point(data = variants.selectedGene.path, aes(x=pos, y=cadd), colour="red", pch=19, alpha = 1) +
      theme_bw() +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major.y = element_line(colour = "black"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()
      ) +
      labs(title=paste("GAVIN ", version, ". Thresholds for ", selectGene, ": pathogenic > ", rules.selectedGene$PathogenicIfCADDScoreGreaterThan, ", benign < ", rules.selectedGene$BenignIfCADDScoreLessThan, "\nMAF benign > ",rules.selectedGene$BenignIfMAFGreaterThan,", cat: ",rules.selectedGene$CalibrationCategory,"\nred: ClinVar pathogenic variants, blue: rarity & impact matched gnomAD\nvariants, green: RefSeq exons, grey: genome-wide fallback threshold", sep="")) +
      ylab("CADD scaled C-score") +
      xlab(paste("Genomic position [",selectGene,", chr. ",sort(unique(variants.selectedGene$chr)),"]", sep="")) +
      theme(legend.position = "none")
  } else {
    p <- ggplot() +
      geom_text(aes(label = paste("There is no variant CADD calibration for ", selectGene, ".\nHowever, estimated impact or allele frequency may\nbe predictive. Please consult the GAVIN rule guide.", sep=""), x = .5, y = .75)) +
      theme_bw() +
      theme(axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()
      ) + 
      labs(title=paste("GAVIN ", version, ". Thresholds for ", selectGene, ": pathogenic > ", rules.selectedGene$PathogenicIfCADDScoreGreaterThan, ", benign < ", rules.selectedGene$BenignIfCADDScoreLessThan, "\nMAF benign > ",rules.selectedGene$BenignIfMAFGreaterThan,", cat: ",rules.selectedGene$CalibrationCategory,"\nred: ClinVar pathogenic variants, blue: rarity & impact matched gnomAD\nvariants, green: RefSeq exons, grey: genome-wide fallback threshold", sep="")) +
      ylab("CADD scaled C-score") +
      xlab(paste("Genomic position [",selectGene,", chr. -]", sep="")) +
      theme(legend.position = "none") +
      xlim(0,1) + ylim(0,1)
  }
  #p
  ggsave(paste("/Users/joeri/Desktop/gavin-r0.5/plots/",selectGene,".png", sep=""), width=8, height=4.5)
}

###################################################
# Bootstrap analysis result processing and plotting
###################################################

bootStrapResults <- paste(pathToGavinGitRepo,"/data/other/performancebootstrap_output_usedinpaper_",version,".R",sep="")
df <- read.table(bootStrapResults,header=TRUE)
df$Acc <- as.double(as.character(df$Acc))

# Figure 2 in paper
ggplot() + 
  geom_boxplot(data = df, aes(x = Label, y = Acc, fill = Calib, colour=Tool), size=1.1) +
  theme_bw() + theme(panel.grid.major = element_line(colour = "black"), axis.text=element_text(size=12),  axis.title=element_text(size=14,face="bold")) +
  ylab("Accuracy") + xlab("GAVIN classification") +
  scale_x_discrete(limits=c("C3_GAVINnocal","C3_GAVIN","C4_GAVINnocal","C4_GAVIN", "C1_C2_GAVINnocal", "C1_C2_GAVIN"),
                   labels = c("C1_C2_GAVIN" = "",
                              "C1_C2_GAVINnocal" = "",
                              "C4_GAVIN" = "", 
                              "C4_GAVINnocal" = "", 
                              "C3_GAVIN" = "",
                              "C3_GAVINnocal"="" )) +
  scale_fill_manual(values=c(blueishgreen, yellow, skyblue), 
                    name="Selected gene group",
                    breaks=c("C1_C2", "C4", "C3"),
                    labels=c("CADD predictive genes (681)", "CADD less predictive genes (732)", "Scarce training data genes (774)")) +
  scale_colour_manual(values=c("black", vermillion), name="GAVIN classification", breaks=c("GAVIN", "GAVINnocal"), labels=c("Gene-specific", "Genome-wide")) +
  coord_flip() +
  guides(colour = guide_legend(override.aes = list(size = .75)), fill = guide_legend(override.aes = list(size = .5)))+
  theme(legend.key = element_blank()) +
  annotate("text", x = 1.2:6.2, y = .75, label = rep(c("Genome-wide", "Gene-specific"), 3))
ggsave("Figure2.pdf", width = 24.5, height = 12.25, units = "cm")

# mann-whitney-wilcoxon test and median accuracy values used in paper
calib <- subset(df, Tool == "GAVIN")
uncalib <- subset(df, Tool == "GAVINnocal")

C1_calib <- subset(calib, Calib == "C1_C2")
C1_uncalib <- subset(uncalib, Calib == "C1_C2")
median(C1_calib$Acc)
median(C1_uncalib$Acc)
wilcox.test(C1_calib$Acc, C1_uncalib$Acc)

C4_calib <- subset(calib, Calib == "C4")
C4_uncalib <- subset(uncalib, Calib == "C4")
median(C4_calib$Acc)
median(C4_uncalib$Acc)
wilcox.test(C4_calib$Acc, C4_uncalib$Acc)

C3_calib <- subset(calib, Calib == "C3")
C3_uncalib <- subset(uncalib, Calib == "C3")
median(C3_calib$Acc)
median(C3_uncalib$Acc)
wilcox.test(C3_calib$Acc, C3_uncalib$Acc)

### r0.3 vs r0.5 ###
variantsv3 <- read.table(gzfile(paste(pathToGavinGitRepo,"/data/other/calibrationvariants_r0.3.tsv.gz",sep="")), sep="\t", header=T)
variantsv5 <- read.table(gzfile(paste(pathToGavinGitRepo,"/data/other/calibrationvariants_r0.5.tsv.gz",sep="")), sep="\t", header=T)
d3 <- as.data.frame(table(variantsv3$gene))
d4 <- as.data.frame(table(variantsv5$gene))
m <- merge(x = d3, y = d4, by.x = "Var1", by.y = "Var1")
plot(m$Freq.x,m$Freq.y, log="yx", ylab="GAVIN r0.5", xlab="GAVIN r0.3", title("Variants used per calibrated gene"))
abline(a = 0, b = 1, col = "red")
m$abs <- m$Freq.y-m$Freq.x
m$rel <- m$Freq.y/m$Freq.x
plot(sort(m$abs))
top50 <- head(m[order(m$abs, decreasing = TRUE),], 50)
#write.table(top50rel, file="~/Desktop/top50rel.tsv", sep="\t")
#genes with variants unique to d3
lost <- d3[which(!(d3$Var1 %in% d4$Var)),]$Var1
length(lost)
#genes with variants  unique to d4
gained <- d4[which(!(d4$Var1 %in% d3$Var)),]$Var1
length(gained)
write.table(lost, file="~/lost.txt", row.names = F, quote=F)
write.table(gained, file="~/gained.txt", row.names = F, quote=F)

### r0.3 vs r0.5 based on rule guide ### 
rules3 <- read.table(paste(pathToGavinGitRepo,"/data/predictions/GAVIN_ruleguide_r0.3.tsv",sep=""),header=TRUE,sep='\t',quote="",comment.char="",as.is=TRUE, skip=33)
rules4 <- read.table(paste(pathToGavinGitRepo,"/data/predictions/GAVIN_ruleguide_r0.5.tsv",sep=""),header=TRUE,sep='\t',quote="",comment.char="",as.is=TRUE, skip=33)
dim(rules3)
dim(rules4)
#genes unique to release 3 (incl. no variants)
lost <- rules3[which(!(rules3$Gene %in% rules4$Gene)),]$Gene
length(lost)
#genes unique to release 4 (incl. no variants)
gained <- rules4[which(!(rules4$Gene %in% rules3$Gene)),]$Gene
length(gained)
gainedInterestingGenes <- subset(rules4[which(rules4$Gene %in% gained),], CalibrationCategory == "C1")
write.table(gainedInterestingGenes[,1], file="~/gainedInterestingGenes.txt", row.names = F, quote=F)

## variant compare ##
ttn3 <- variantsv3[which(variantsv3$gene=="TTN"),]
ttn5 <- variantsv5[which(variantsv5$gene=="TTN"),]
ttnM <- merge(x = ttn3, y = ttn5, by.x = c("chr", "pos", "ref", "alt"), by.y = c("chr", "pos", "ref", "alt"))
plot(ttnM$cadd.x, ttnM$cadd.y)

## which in rules but not in refseq?
notinrefseq <- rules4[which(!(rules4$Gene %in% refeqExons$gene)),]$Gene
## which in variants but not in refseq?
notinrefseq <- d3[which(!(d3$Var1 %in% refeqExons$gene)),]$Var1

#do C-cat genes match genes with variants ?
Cgenes <- rules4[which(startsWith(rules4$CalibrationCategory, "C")),]
dim(Cgenes)[[1]]==dim(d4)[[1]]
# artifacts?
Agenes <- rules4[which(rules4$CalibrationCategory=="C5"),]


## trackview?
source("https://bioconductor.org/biocLite.R")
biocLite("trackViewer")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
biocLite("org.Hs.eg.db")
library(trackViewer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

gr <- GRanges("chr11", IRanges(122929275, 122930122), strand="-")
trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                         org.Hs.eg.db,
                         gr=gr)
viewerStyle <- trackViewerStyle()
setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .05, .02, .02))
vp <- viewTracks(trackList(trs), 
                 gr=gr, viewerStyle=viewerStyle, 
                 autoOptimizeStyle=TRUE)

dtTrack <- DataTrack(start=seq(1,1000, len=100), width=10, data=dat,
                     chromosome=1, genome="mm9", name="random data")
