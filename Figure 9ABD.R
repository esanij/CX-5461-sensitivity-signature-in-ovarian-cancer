library(GSVA)
library(pheatmap)
library(reshape2)

load("human_coding_genes_GRCh37.Rdata")  ##Data frame with human coding genes

##Metadata from ICGC
donor_AU <- read.delim("donor.tsv")   
specimen_AU <- read.delim("specimen.tsv")
sample_AU <- read.delim("sample.tsv")
specimen_AU$donor_tumour_state_at_diagnosis <- donor_AU[match(specimen_AU$submitted_donor_id, donor_AU$submitted_donor_id), "donor_tumour_stage_at_diagnosis"]

##Expression data from ICGC
exp_AU <- read.delim("exp_seq.relevant_columns.tsv", sep=" ")

exp_AU$submitted_specimen_id <- specimen_AU[match(exp_AU$icgc_specimen_id, specimen_AU$icgc_specimen_id), "submitted_specimen_id"]
exp_AU$submitted_donor_id <- specimen_AU[match(exp_AU$icgc_donor_id, specimen_AU$icgc_donor_id), "submitted_donor_id"]
exp_AU$specimen_type <- specimen_AU[match(exp_AU$submitted_specimen_id, specimen_AU$submitted_specimen_id), "specimen_type"]

exp_AU <- exp_AU[,c("submitted_donor_id", "submitted_specimen_id", "specimen_type",
                    "gene_id", "normalized_read_count")]

coding_genes$Gene_name <- toupper(coding_genes$Gene_name)
exp_AU <- exp_AU[exp_AU$gene_id %in% coding_genes$Gene_ensembl, ]
exp_AU$gene_id <- coding_genes$Gene_name[match(exp_AU$gene_id, coding_genes$Gene_ensembl)]

#Select only primary tumours from solid tissue
exp_AU_primary <- exp_AU[exp_AU$specimen_type == "Primary tumour - solid tissue",]
exp_AU_primary <- dcast(exp_AU_primary, gene_id ~ submitted_specimen_id, value.var = "normalized_read_count", fun.aggregate = median)
rownames(exp_AU_primary) <- exp_AU_primary$gene_id
exp_AU_primary$gene_id <- NULL


##Select recurrent tumours. Prefix -14 means a second specimen from the patient (first one was -13)
exp_AU_recurrent <- exp_AU[exp_AU$specimen_type == "Recurrent tumour - other",]
exp_AU_recurrent <- dcast(exp_AU_recurrent, gene_id ~ submitted_specimen_id, value.var = "normalized_read_count", fun.aggregate = median)
rownames(exp_AU_recurrent) <- exp_AU_recurrent$gene_id
exp_AU_recurrent$gene_id <- NULL

##Load in MYC and BRCAm signatures
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
MYC_signature <- trim(toupper(as.character(read.delim("MYC_signature.txt", header=FALSE)$V1)))

BRCAm_signature <- trim(toupper(as.character(read.delim("BRCA1_mutated_signature.txt", header=FALSE)$V1)))

MYC_signature[which(!(MYC_signature %in% rownames(exp_AU_primary)))]
BRCAm_signature[which(!(BRCAm_signature %in% rownames(exp_AU_primary)))]

##Matching gene IDs
MYC_signature[MYC_signature == "FAM100A"] <- "UBALD1"
MYC_signature[MYC_signature == "C1ORF107"] <- "DIEXF"
MYC_signature[MYC_signature == "FAM86C"] <- "FAM86C1"
MYC_signature[MYC_signature == "DKFZP686O24166"] <- "NCR3LG1"

BRCAm_signature[BRCAm_signature == "SFRS11"] <- "SRSF11"
BRCAm_signature[BRCAm_signature == "SFRS13A"] <- "SRSF10"
BRCAm_signature[BRCAm_signature == "DDX39"] <- "DDX39A"

primary_ssGSEA <- gsva(as.matrix(exp_AU_primary), 
                       list(MYC_sig = MYC_signature, BRCAm_sig=BRCAm_signature), 
                       method='ssgsea', kcdf='Gaussian', ssgsea.norm=TRUE)

recurrent_ssGSEA <- gsva(as.matrix(exp_AU_recurrent), 
                         list(MYC_sig = MYC_signature, BRCAm_sig=BRCAm_signature), 
                         method='ssgsea', kcdf='Gaussian', ssgsea.norm=TRUE)

save(primary_ssGSEA, file="primary_ssGSEA.Rdata")
save(recurrent_ssGSEA, file="recurrent_ssGSEA.Rdata")

#Convert ssGSEA to 0-1 range
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

primary_ssGSEA <- range01(primary_ssGSEA)
recurrent_ssGSEA <- range01(recurrent_ssGSEA)

###Adding clinical data to the patient heatmap
clinical_data <- donor_AU[,c("submitted_donor_id", "disease_status_last_followup", "donor_tumour_stage_at_diagnosis")]
rownames(clinical_data) <- paste(clinical_data$submitted_donor_id, "-1", sep="")
clinical_data$submitted_donor_id <- NULL

clinical_data$disease_status_last_followup <- as.character(clinical_data$disease_status_last_followup)
clinical_data$disease_status_last_followup[clinical_data$disease_status_last_followup == "complete remission"] <- "complete_remission"
colnames(clinical_data) <- c("Disease_status", "Tumour_stage")

anno_colors <- list(Tumour_stage = c(III="goldenrod1", IV="royalblue3"),
                    Disease_status = c(progression="darkgreen", relapse="firebrick1", complete_remission="darkviolet"))

pdf("Heatmap_MYC_BRCAm_primary.pdf", width=25, height=2.5)
pheatmap(primary_ssGSEA, show_colnames=FALSE, annotation_col=clinical_data, 
         annotation_colors = anno_colors)
dev.off()


##Save clinical data
clinical_data_recurrent <- clinical_data
rownames(clinical_data_recurrent) <- paste(rownames(clinical_data_recurrent), 3, sep="")
clinical_data_recurrent2 <- clinical_data
rownames(clinical_data_recurrent2) <- paste(rownames(clinical_data_recurrent2), 4, sep="")
clinical_data_recurrent <- rbind(clinical_data_recurrent, clinical_data_recurrent2)
patients <- rownames(clinical_data_recurrent)
clinical_data_recurrent <- data.frame(Disease_status = "relapse",
                                      Tumour_stage=clinical_data_recurrent$Tumour_stage)
rownames(clinical_data_recurrent) <- patients

save(clinical_data, file="clinical_data.Rdata")
save(clinical_data_recurrent, file="clinical_data_recurrent.Rdata")

pdf("Heatmap_MYC_BRCAm_relapse.pdf", width=6.7, height=1.9)
pheatmap(recurrent_ssGSEA, show_colnames=FALSE, annotation_col=clinical_data_recurrent, 
         annotation_colors = anno_colors, border_color = NA)
dev.off()


##Collating data of patients with primary and relapse data
both_ssGSEA <- cbind(primary_ssGSEA, recurrent_ssGSEA)
relapse_patients <- colnames(both_ssGSEA)
relapse_patients <- relapse_patients[(grep(paste(c("-13$","-14$"),collapse="|"), relapse_patients))]
relapse_patients <- substr(relapse_patients, 1, 8)

primary_patients <- colnames(both_ssGSEA)
primary_patients <- primary_patients[(grep("-1$", primary_patients))]
primary_patients <- substr(primary_patients, 1, 8)

relapse_patients <- relapse_patients[relapse_patients %in% primary_patients]

both_ssGSEA_patients_relapsed <- both_ssGSEA[,grep(paste(relapse_patients,collapse="|"), colnames(both_ssGSEA))]

#Based on previous heatmaps
patients_relapse_sensitive <- c("AOCS-093-13","AOCS-064-13","AOCS-137-13","AOCS-088-13","AOCS-150-13","AOCS-150-14",
                                "AOCS-135-13","AOCS-141-13","AOCS-086-13","AOCS-134-13","AOCS-117-13","AOCS-138-13",
                                "AOCS-167-13","AOCS-135-13","AOCS-095-13","AOCS-091-13")

patients_relapse_sensitive <- c("AOCS-088", "AOCS-064", "AOCS-093")


both_ssGSEA_patients_sensitive <- both_ssGSEA[,grep(paste(patients_relapse_sensitive,collapse="|"), colnames(both_ssGSEA))]

temp <- clinical_data_recurrent
clinical_data_recurrent$Disease_status <- "relapse"
clinical_data_recurrent$Sample_type <- "Relapse"
clinical_data$Sample_type <- "Primary"

clinical_data_both <- rbind(clinical_data, clinical_data_recurrent)
colnames(clinical_data_both) <- c("disease_status_last_followup",
                                  "donor_tumour_stage_at_diagnosis",
                                  "Sample_type")


pdf("Relapse_patients_sensitive_heatmap.pdf", height=4, width=8)
pheatmap(both_ssGSEA_patients_sensitive, annotation_col=clinical_data_both, border_color = NA)
dev.off()

