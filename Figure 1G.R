library(GSVA)
library(reshape2)
library(ggplot2)
library(readr)

#Table of cell lines indicating which are sensitive and resistant, and whether they are HGSC
sampleinfo <- read_tsv("sampleinfo.txt") 

#Read in expression data from cell lines
expression <- read_csv("NoAKTSGKPfizer Expression Normalised Batch Adjusted Intensity Data.csv")

# keep just the cell lines from the OVCA spreadsheet, as discussed with Elaine 
expression <- expression[, c("Probe.Set.ID","Gene.Symbol","Gene.name", sampleinfo$Cell.line)]

expression <- as.data.frame(expression)
expression <- expression[expression$Gene.Symbol != "---",]
expression <- expression[, -c(1,3)]
genes <- expression[,1]
expression <- log2(expression[,2:ncol(expression)])
expression <- cbind(genes, expression)
expression_ag <- aggregate(. ~ genes, expression, mean)
rownames(expression_ag) <- expression_ag$genes
expression_ag <- expression_ag[,-1]

sampleinfo <- as.data.frame(sampleinfo)

###HRD signature from Peng et al. 2014
HRD_signature <- read.delim("HRD signature Peng et al. 2014.txt", header=FALSE)
HRD_signature <- as.character(HRD_signature[,1])
HRD_signature <- list(HRD = HRD_signature)

ssGSEA_results <- gsva(as.matrix(expression_ag), HRD_signature, method="ssgsea", verbose=FALSE, kcdf="Gaussian", ssgsea.norm=TRUE)

ssGSEA_results_melt <- melt(ssGSEA_results)
colnames(ssGSEA_results_melt) <- c("Pathway", "Cell_line", "ssGSEA")
ssGSEA_results_melt$Sample_type <- sampleinfo[match(ssGSEA_results_melt$Cell_line, sampleinfo$Cell.line), "Phenotype"]


ssGSEA_results_melt$Sample_type <- factor(ssGSEA_results_melt$Sample_type, levels=c("normal", "sensitive", "resistant"))

ssGSEA_results_melt <- ssGSEA_results_melt[ssGSEA_results_melt$Sample_type != "normal",]

pdf("HRD signature.pdf", width=2.7, height=4)
g <- ggplot(ssGSEA_results_melt, aes(x=Sample_type, y = ssGSEA))+
  geom_boxplot(aes(fill=Sample_type))+
  scale_fill_manual(values=c(sensitive="steelblue", resistant="brown2"))+
  facet_grid(.~Pathway)+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(g)
dev.off()

wilcox.test(ssGSEA_results_melt$ssGSEA[ssGSEA_results_melt$Sample_type == "sensitive"],
            ssGSEA_results_melt$ssGSEA[ssGSEA_results_melt$Sample_type == "resistant"])

