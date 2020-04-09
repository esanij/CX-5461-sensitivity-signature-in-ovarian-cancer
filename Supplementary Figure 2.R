###CX5461 and PARPi sensitivity as single agents

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

unique(cell_line_details$GDSC.Tissue.descriptor.1)
unique(cell_line_details$GDSC.Tissue.descriptor.2)

cell_line_details_ov <- cell_line_details[cell_line_details$GDSC.Tissue.descriptor.2 == "ovary",] #43

sensitivity_ov <- sensitivity[sensitivity$CELL_LINE_NAME %in% cell_line_details_ov$Sample.Name,]
unique(sensitivity_ov$CELL_LINE_NAME)  #43

sensitivity_CX54641 <- sensitivity_ov[which(sensitivity_ov$DRUG_NAME == "CX-5461"),]
sensitivity_olaparib <- sensitivity_ov[which(sensitivity_ov$DRUG_NAME == "Olaparib"),]
sensitivity_rucaparib <- sensitivity_ov[which(sensitivity_ov$DRUG_NAME == "Rucaparib"),]
sensitivity_veliparib <- sensitivity_ov[which(sensitivity_ov$DRUG_NAME == "Veliparib"),]
sensitivity_talazoparib <- sensitivity_ov[which(sensitivity_ov$DRUG_NAME == "Talazoparib"),]

ov_cell_lines_GDS <- as.character(unique(sensitivity_ov$CELL_LINE_NAME))

drugs_CX_parpi <- data.frame(Cell_lines = ov_cell_lines_GDS,
                             CX5461 = sensitivity_CX54641[match(ov_cell_lines_GDS, sensitivity_CX54641$CELL_LINE_NAME), "LN_IC50"],
                             Olaparib = sensitivity_olaparib[match(ov_cell_lines_GDS, sensitivity_olaparib$CELL_LINE_NAME), "LN_IC50"],
                             Rucaparib = sensitivity_rucaparib[match(ov_cell_lines_GDS, sensitivity_rucaparib$CELL_LINE_NAME), "LN_IC50"],
                             Veliparib = sensitivity_veliparib[match(ov_cell_lines_GDS, sensitivity_veliparib$CELL_LINE_NAME), "LN_IC50"],
                             Talazoparib = sensitivity_talazoparib[match(ov_cell_lines_GDS, sensitivity_talazoparib$CELL_LINE_NAME), "LN_IC50"])


HGSOV <- c("59M","CAOV3", "CAOV4","COV318",
           "COV362","FUOV1","HEYA8","JHOS2",
           "JHOS4","KURAMOCHI","OVCAR3","OAW28",
           "OVCAR4","OVCAR8","OVKATE",
           "OVSAHO","SNU119","TYKNU")



drugs_CX_parpi$Cell_lines <- gsub("-", "", drugs_CX_parpi$Cell_lines)
drugs_CX_parpi$Cell_lines <- toupper(drugs_CX_parpi$Cell_lines)

drugs_CX_parpi$HGSOV <- ifelse(drugs_CX_parpi$Cell_lines %in% HGSOV, "HGSOV", "NOT")

drugs_CX_parpi <- drugs_CX_parpi[drugs_CX_parpi$HGSOV == "HGSOV",]

cor(drugs_CX_parpi$CX5461, drugs_CX_parpi$Olaparib, use="pairwise.complete.obs")
cor(drugs_CX_parpi$CX5461, drugs_CX_parpi$Rucaparib, use="pairwise.complete.obs")
cor(drugs_CX_parpi$CX5461, drugs_CX_parpi$Talazoparib, use="pairwise.complete.obs")

cor.test(drugs_CX_parpi$CX5461, drugs_CX_parpi$Olaparib, use="pairwise.complete.obs")$p.value
cor.test(drugs_CX_parpi$CX5461, drugs_CX_parpi$Rucaparib, use="pairwise.complete.obs")$p.value
cor.test(drugs_CX_parpi$CX5461, drugs_CX_parpi$Talazoparib, use="pairwise.complete.obs")$p.value

#drugs_CX_parpi <- drugs_CX_parpi[drugs_CX_parpi$Cell_lines != "OVCAR8",]

ggplot(drugs_CX_parpi, aes(x=CX5461, y = Olaparib))+
  geom_point()+
  geom_smooth(method='lm',formula=y~x)+
  theme_bw()

ggplot(drugs_CX_parpi, aes(x=CX5461, y = Rucaparib))+
  geom_point()+
  geom_smooth(method='lm',formula=y~x)+
  theme_bw()

ggplot(drugs_CX_parpi, aes(x=CX5461, y = Veliparib))+
  geom_point()+
  geom_smooth(method='lm',formula=y~x)+
  theme_bw()

ggplot(drugs_CX_parpi, aes(x=CX5461, y = Talazoparib))+
  geom_point()+
  geom_smooth(method='lm',formula=y~x)+
  theme_bw()

