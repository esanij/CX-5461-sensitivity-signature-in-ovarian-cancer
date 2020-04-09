library(affy)
library(limma)
targets <- readTargets("./HuGene-1_0-st-v1/all_chips.txt")
data.raw <- ReadAffy(filenames=targets$FileName, celfile.path="./HuGene-1_0-st-v1", cdfname="hugene10stv1.r3.xdacdf")
data.norm <- rma(data.raw)

cont <- F
skip <- 0
while(!cont) {
  txt <- scan("HuGene-1_0-st-v1.CURRENT.transcript.csv", nlines=1, skip=skip, comment="", what="character")
  if(regexpr("##|#%", txt[1]) == -1) {
    cont <- T
  } else skip <- skip +1
}
full.annot <- read.csv("./Annotation/HuGene-1_0-st-v1/HuGene-1_0-st-v1.CURRENT.transcript.csv", as.is=T, comment="", skip=skip, head=T, quote="\"")
design <- model.matrix(~ -1+factor(targets$CellLine)+factor(targets$Batch))
colnames(design) <- sub("factor\\(targets\\$CellLine\\)", "", colnames(design))
colnames(design) <- make.names(sub("factor\\(targets\\$Batch\\)", "Batch", colnames(design)))
fit <- lmFit(data.norm, design)