##  Create QQ plots, Volcano Plots, HeatMap of sign of logFC,
##  tSNPE or PCA Cluster based on sign of log FC
##  SHARED SIGNIFICANT GSs
library(qvalue)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tibble)

source("/gpfs/data/stranger-lab/askol/TCGA2/DE_Surv/Code/DE_Surv_Plots_funcs.r")
    
DE_Surv_Plots <- function(project) {           


   
    ## GET GTEX DE RESULTS ##
    tmp <- get.GTEx.results()
    
    logFC.gtex <- tmp$logFC
    logPs.gtex <- tmp$logPs
    Qs.gtex <- tmp$Qs
    geneInfo.gtex <- tmp$GeneInfo
    tissues <- tmp$tissues
    

    ##
    ## GET TCGA DE RESULTS 
    ##
    ## GET TCGA PROJECTS ##   
    
    tmp <- get.TCGA.results(project)
    
    logFC.tcga <- tmp$logFC
    logPs.tcga <- tmp$logPs
    Qs.tcga <- tmp$Qs
    geneInfo.tcga <- tmp$GeneInfo
    projects <- tmp$projects

    ## #########
    source("/gpfs/data/stranger-lab/askol/TCGA2/DE_Surv/Code/DE_Surv_Plots_funcs.r")
    
    RsltDir <- "/gpfs/data/stranger-lab/askol/TCGA2/Survival/Results/"

    ## GET SURVIVAL DATA ##
    rslts <- compile.surv.results(project, RsltDir)

    ps  <- rslts[[1]]
    qs <- rslts[[2]]
    coef <- rslts[[3]]
    
    ## remove genes with no symbols
    for (i in 1:length(ps)){

        ps[[i]] <- filter(ps[[i]], !is.na(SYMBOL))
        qs[[i]] <- filter(qs[[i]], !is.na(SYMBOL))
        coef[[i]] <- filter(coef[[i]], !is.na(SYMBOL))
    }

    ## PLOT SUMMARIZIES OF SURVIVAL P-VALUES FOR DE GENES (GTEX AND TCGA)
    PlotDir <- paste0("/gpfs/data/stranger-lab/askol/TCGA2/DE_Surv/Results/Plots/")
    
    out <- summary.de.surv.plot(tissues, project, 
                                ps, qs,
                                logPs.gtex, logFC.gtex,
                                logPs.tcga, logFC.tcga,
                                geneInfo.gtex, geneInfo.tcga, 
                                PlotDir,
                                niter = 10000)
    
    
    out <- summary.de.surv.plot(projects, project,
                                ps, qs,
                                logPs.tcga, logFC.tcga,
                                logPs.tcga, logFC.tcga,
                                geneInfo.tcga, geneInfo.tcga,
                                PlotDir,
                                niter=10000)
      

}
## ---- MAIN -----##

args <- commandArgs(TRUE)
project <- args[1]
DE_Surv_Plots(project)


    
