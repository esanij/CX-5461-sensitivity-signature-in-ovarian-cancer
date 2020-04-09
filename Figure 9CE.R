library(pheatmap)
library(ggplot2)
library(ggalt)
library(ggrepel)

#Generated in "Figure 9ABD.R" script
load("primary_ssGSEA.Rdata")
load("recurrent_ssGSEA.Rdata")
load("clinical_data.Rdata")
load("clinical_data_recurrent.Rdata")

#Data generated in "Figure 1H.R"
load("cell_line_ssGSEA.Rdata")
cell_line_ssGSEA <- sig_ssGSEA
load("cell_line_annotation.Rdata")
cell_line_annot <- annot_cat

##Converting to 0-1 range
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

primary_ssGSEA_01 <- range01(primary_ssGSEA)
recurrent_ssGSEA_01 <- range01(recurrent_ssGSEA)
cell_line_ssGSEA_01 <- range01(cell_line_ssGSEA)

ssGSEA_df <- t(cbind(primary_ssGSEA_01, recurrent_ssGSEA_01, cell_line_ssGSEA_01))
ssGSEA_df <- as.data.frame(ssGSEA_df)
ssGSEA_df$Sample_type <- NA

ssGSEA_df$Sample_type[grep("-1$", rownames(ssGSEA_df))] <- "Primary"
ssGSEA_df$Sample_type[grep("-13$", rownames(ssGSEA_df))] <- "Relapse"
ssGSEA_df$Sample_type[grep("-14$", rownames(ssGSEA_df))] <- "Relapse"

sensitive_cl <- rownames(cell_line_annot)[cell_line_annot$CX5461 == "S"]
resistant_cl <- rownames(cell_line_annot)[cell_line_annot$CX5461 == "R"]

ssGSEA_df$Sample_type[grep(paste(sensitive_cl, collapse="|"), rownames(ssGSEA_df))] <- "Sensitive_cell_line"
ssGSEA_df$Sample_type[grep(paste(resistant_cl, collapse="|"), rownames(ssGSEA_df))] <- "Resistant_cell_line"

##K-means clustering
ssGSEA_df_kmeans <- subset(ssGSEA_df, Sample_type != "Primary")
ssGSEA_df_kmeans$Sample_type <- NULL

#Determine the number of clusters
wss <- (nrow(ssGSEA_df_kmeans)-1)*sum(apply(ssGSEA_df_kmeans,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(ssGSEA_df_kmeans, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
set.seed(123)
fit <- kmeans(ssGSEA_df_kmeans, 4, nstart=50) # 5 cluster solution
# get cluster means 
aggregate(ssGSEA_df_kmeans,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(ssGSEA_df_kmeans, Cluster=fit$cluster, Sample_name = rownames(ssGSEA_df_kmeans))
mydata$Cluster <- as.factor(mydata$Cluster)

mydata$Sample_type <- ssGSEA_df$Sample_type[match(rownames(mydata), 
                                                  rownames(ssGSEA_df))]

pdf("Relapse_patients_scatterplot.pdf", width=9, height=7)
g <- ggplot(mydata, aes(x=MYC_sig, y=BRCAm_sig))+
  geom_point(aes(color=Sample_type), size=4)+
  geom_encircle(data = subset(mydata, Cluster==1), aes(x=MYC_sig, y=BRCAm_sig), alpha=0.4) +
  geom_encircle(data = subset(mydata, Cluster==2), aes(x=MYC_sig, y=BRCAm_sig), alpha=0.4) +
  geom_encircle(data = subset(mydata, Cluster==3), aes(x=MYC_sig, y=BRCAm_sig), alpha=0.4) +
  geom_encircle(data = subset(mydata, Cluster==4), aes(x=MYC_sig, y=BRCAm_sig), alpha=0.4) +
  theme_bw()
print(g)
dev.off()


ssGSEA_df_relapse_sensitive <- ssGSEA_df

patients_relapse_sensitive <- c("AOCS-088", "AOCS-064", "AOCS-093", "AOCS-137")
ssGSEA_df_relapse_sensitive <- ssGSEA_df_relapse_sensitive[grep(paste(patients_relapse_sensitive, collapse="|"),
                                                                rownames(ssGSEA_df_relapse_sensitive)),]
ssGSEA_df_relapse_sensitive <- rbind(ssGSEA_df_relapse_sensitive, subset(ssGSEA_df, Sample_type %in% c("Resistant_cell_line", "Sensitive_cell_line")))
ssGSEA_df_relapse_sensitive$Sample_name <- rownames(ssGSEA_df_relapse_sensitive)

ssGSEA_df_relapse_sensitive$Sample_name[!(substr(ssGSEA_df_relapse_sensitive$Sample_name,1,8) %in% patients_relapse_sensitive)] <- NA

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

gg_color_hue(3)

pdf("Sensitive_primary_and_relapse_patients_scatterplot.pdf", width=7, height=5)
g <- ggplot(ssGSEA_df_relapse_sensitive, aes(x=MYC_sig, y=BRCAm_sig))+
  geom_point(aes(colour=Sample_type), size=3)+
  geom_text_repel(aes(label=Sample_name))+
  scale_color_manual(values=c(Relapse=gg_color_hue(3)[1],
                              Resistant_cell_line=gg_color_hue(3)[2],
                              Sensitive_cell_line=gg_color_hue(3)[3],
                              Primary=gg_color_hue(4)[4]))+
  theme_bw()
print(g)
dev.off()

