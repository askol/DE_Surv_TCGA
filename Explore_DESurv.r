library(dplyr)
library(reshape2)
library(annotables)
library(RColorBrewer)
library(scales) ## for muted
library(cowplot)
library(ggnewscale)

source("/gpfs/data/stranger-lab/askol/TCGA2/DE_Surv/Code/Explore_DESurv_funcs.r")
PlotDir <- "/gpfs/data/stranger-lab/askol/TCGA2/DE_Surv/Results/Plots/"

## READ IN RESULTS FROM DE_SURV_PLOTS.R ##

geneInfo <- grch38
d <- as.data.frame(table(geneInfo$symbol))
names(d)[2] <- "GeneFreq"
geneInfo <- left_join(geneInfo, d, by=c("symbol" = "Var1"))
geneInfo <- geneInfo[-which(duplicated(geneInfo$symbol)),]

ResultDir <- "/gpfs/data/stranger-lab/askol/TCGA2/DE_Surv/Results/"

data <- collect.data(ResultDir)
data <- mutate(data, DEPop = ifelse(grepl("TCGA",tissue), "TCGA" , "GTEx"))
data <- data %>% arrange(emppmn, Psurv)

## GENES WITH SIGNIFICANT SURVIVAL (<.01) AND QUANT = 0.01 AND NO X AND (BOTH OR NOT)##
plot.sharing.all(data, PlotDir)


## CREATE A HEATMAP THAT SHOWS WHICH TISSUE/PROJECT FOUND SEX DE GENES THAT SHOW
## ENRICHMENT FOR SMALL SURVIVAL P-VALUE IN CANCER ##
data <- collect.data.enrichment(ResultDir)
plot.enrich.all(data, PlotDir)



