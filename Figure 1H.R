library(GSVA)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(pheatmap)

##Data from GDS
##Cell line details from GDS
cell_line_details <- read.delim("Cell_Lines_Details.txt")

##Drug sensitivity from GDS
sensitivity <- read.delim("v17.3_fitted_dose_response.txt")

order_paper <- c("OVCAR3","FUOV1","COLO704","TOV112D","CH1",
                 "A2780","OVCAR5","2008","RMG1","ES2","OVCA432",
                 "OAW42","SKOV3","TOV21G","MCAS","Caov3","Colo720E",
                 "EFO27","IGROV1","OVCAR4","JHOC5","JHOC9","OVCAR8",
                 "JHOM1","KURAMOCHI","59M","EFO21","JHOS3","JHOC7",
                 "RMGII","OAW28","OV90")

order_paper <- toupper(order_paper)

cell_line_details$Sample.Name <- toupper(cell_line_details$Sample.Name)

cell_line_details$Sample.Name <- gsub("-", "", cell_line_details$Sample.Name)

cell_line_details_ov <- cell_line_details[cell_line_details$Sample.Name %in% order_paper,]

##Load in MYC and BRCAm signatures
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
MYC_signature <- trim(toupper(as.character(read.delim("MYC_signature.txt", header=FALSE)$V1)))
BRCAm_signature <- trim(toupper(as.character(read.delim("BRCA1_mutated_signature.txt", header=FALSE)$V1)))

##ssGSEA to estimate activity
load("CCLE_exp_agg_ovarian.Rdata")  ##Columns are the cell lines. Rows are the genes
load("human_coding_genes_GRCh37.Rdata")  ##Data frame of human coding genes

CCLE_exp_agg <- CCLE_exp_agg[rownames(CCLE_exp_agg) %in% coding_genes$Gene_name,]

MYC_signature[which(!(MYC_signature %in% rownames(as.matrix(CCLE_exp_agg))))]
BRCAm_signature[which(!(BRCAm_signature %in% rownames(as.matrix(CCLE_exp_agg))))]

##Matching gene IDs
MYC_signature[MYC_signature == "C1ORF51"] <- "CIART"
MYC_signature[MYC_signature == "FAM100A"] <- "UBALD1"
MYC_signature[MYC_signature == "C11ORF48"] <- "LBHD1"
MYC_signature[MYC_signature == "C1ORF107"] <- "UTP25"
MYC_signature[MYC_signature == "C10ORF2"] <- "TWNK"
MYC_signature[MYC_signature == "DKFZP686O24166"] <- "NCR3LG1"
BRCAm_signature[BRCAm_signature == "DDX39"] <- "DDX39A"

##Based on experimental results
sensitive <- c("OVCAR3", "FUOV1", "COLO704", "TOV112D", "CH1", "A2780", "OVCAR5",
               "2008", "RMG1", "ES2", "OVCA432")

resistant <- c("JHOC9", "JHOM1", "KURAMOCHI", "59M", "EFO21", "JHOS3",
               "JHOC7", "RMGII", "OAW28", "OV90")

CCLE_exp_agg <- CCLE_exp_agg[, colnames(CCLE_exp_agg) %in% c(sensitive, resistant)]

sig_ssGSEA <- gsva(as.matrix(CCLE_exp_agg), list(MYC_sig = MYC_signature, BRCAm_sig=BRCAm_signature), method='ssgsea', kcdf='Gaussian', ssgsea.norm=TRUE)

annot_cat <- data.frame(Cell_lines= c(sensitive, resistant),
                        CX5461=NA)
annot_cat$CX5461 <- ifelse(annot_cat$Cell_lines %in% sensitive, "S", "R")
rownames(annot_cat) <- annot_cat$Cell_lines
annot_cat$Cell_lines <- NULL

#save(sig_ssGSEA, file="cell_line_ssGSEA.Rdata")
#save(annot_cat, file="cell_line_annotation.Rdata")

##Convert to 0-1 range
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
sig_ssGSEA <- range01(sig_ssGSEA)


##Keep only cell lines with annotation

pdf("Heatmap_of_MYC_BRCAm_signatures_CCLE.pdf", 
    width=4, height=2.5)
pheatmap(sig_ssGSEA, annotation_col=annot_cat, border_color = NA)
dev.off()

#Comparing levels of the MYC and BRCAm signature
sig_ssGSEA_df <- data.frame(Cell_line =colnames(sig_ssGSEA),
                            MYC_sig = sig_ssGSEA[1,],
                            BRCAm_sig = sig_ssGSEA[2,],
                            CX5461 = annot_cat[match(colnames(sig_ssGSEA), rownames(annot_cat)), "CX5461"])

sig_ssGSEA_df_S <- sig_ssGSEA_df[sig_ssGSEA_df$CX5461 == "S",]
sig_ssGSEA_df_R <- sig_ssGSEA_df[sig_ssGSEA_df$CX5461 == "R",]

wilcox.test(sig_ssGSEA_df_S$MYC_sig, sig_ssGSEA_df_R$MYC_sig, alternative="greater")
wilcox.test(sig_ssGSEA_df_S$BRCAm_sig, sig_ssGSEA_df_R$BRCAm_sig, alternative="greater")


g_brcam <- glm(formula = CX5461 ~ BRCAm_sig, 
               family = binomial(link="logit"), data = sig_ssGSEA_df)

g_brcam_myc <- glm(formula = CX5461 ~ BRCAm_sig+MYC_sig, 
                   family=binomial(link="logit"), data = sig_ssGSEA_df)


anova(g_brcam_myc, g_brcam, test="Chisq")


