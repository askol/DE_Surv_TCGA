
check.sample.size <- function(tissue){

    OrigDataDir <- paste0("/gpfs/data/gtex-group/sex_biased_regulation_v8/",
                          "sexDE_v8_final/meri/data/")

    count_file <- paste0(OrigDataDir, "Phenotypes/",tissue,
                         ".v8.gene_counts.txt.gz")
    covs_file <- paste0(OrigDataDir, "/Covariates/covs.basalcovs.",tissue,
        ".txt")
    sex.info <- read.table(file = covs_file, as.is=T, header=F, quote="",
                           comment.char = "", sep="\t")
    colnames(sex.info) = c('SUBJID', 'SEX','SMTSISCH','SMRIN','AGE')
    samps <- read.table(file = count_file, nrow=1, header=F, as.is=T)
    samps <- gsub("\\.","-",samps)
    
    males <- samps[samps %in% sex.info$SUBJID[sex.info$SEX=="M"]]
    females <- samps[samps %in% sex.info$SUBJID[sex.info$SEX=="F"]]

    wcs <- list(males = length(males), females=length(females))

    return(wcs)
}

get.GTEx.results <- function(){
    
    source("/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Code/Summary_funcs.r")
    
    ## GET GTEX TISSUES ##
    tissues <- read.table(file = paste0("/gpfs/data/gtex-group/sex_biased_regulation_v8/",
                              "data/support_files/all_v8_tissues_both_sexes.txt"), as.is=T)
    tissues <- unlist(tissues)
    tissues.keep <- c()
    Ns <- c()
    for (tissue in tissues){
        samp.size <- check.sample.size(tissue)
        
        if (min(unlist(samp.size)) < 40){next }
        tissues.keep <- c(tissues.keep, tissue)
        Ns <- rbind(Ns, c(tissue, samp.size, round(mean(unlist(samp.size))),
                          round(sum(sqrt(unlist(samp.size))))))
    }
    Ns <- as.data.frame(Ns, stringAsFactors=FALSE)
    names(Ns) <- c("tissue", "males", "females", "mean", "mean.sqrt")
    tissues <- tissues.keep

    ## REMOVE LATER !!!! ##
    ##                   ##
    ## tissues <- tissues[1:4]
    
    ##
    ## BRING IN GTEX DE RESULTS
    ##

    ResultDir <- "/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Results/"
    
    tmp <- collect.results(tissues, ResultDir)
    tmp[["tissues"]] <- tissues
    
    return(tmp)
}
    
get.TCGA.results <- function(project){

    projects = read.table(file = "/gpfs/data/stranger-lab/askol/TCGA/TCGA_Target_projects.txt",
        header=TRUE, as.is=TRUE, sep="\t")
    projects = projects[grep("TCGA",projects[,1]),1]
    
    skip.cancers <- c("TCGA-BRCA","TCGA-OV", "TCGA-CESC", "TCGA-PRAD", "TCGA-TGCT",
                      "TCGA-UCS", "TCGA-UCEC")
    projects <- projects[!projects %in% skip.cancers]

    ## !! REMOVE LATER !! ##
    ## projects <- projects[1:4]
    ## projects <- unique(c(projects, project))
    
    ResultDir <- "/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Results/"
    SummaryDir <- "/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Summary/"
    
    setwd(SummaryDir)
    source("/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Code/Summary_funcs.r")
    
    tmp <- collect.results(projects, ResultDir)

    tmp[["projects"]] <- projects
    return(tmp)
}


collect.results <- function(projects, ResultDir){
    
    logFC <- c()
    logPs <- c()
    logQs <- c()
    geneInfo <- c()
    for (project in projects){
        
        print(paste0("Working on project : ",project))
        rslt <- load.results(project, ResultDir)
        dupe.ind <- which(duplicated(rslt$SYMBOL))
        if (length(dupe.ind) > 0){
            rslt <- rslt [-dupe.ind,]
        }
        print(paste0("Removing ",length(dupe.ind), " duplicate gene"))
        logFCtmp <- rslt[,c("SYMBOL", "logFC")]
        logPstmp <- rslt[,c("SYMBOL", "P.Value")]
        names(logFCtmp)[2] = names(logPstmp)[2] = project        

        gi <-  rslt[, c("ENSEMBL","SYMBOL", "chr", "start","end","biotype")]
        if (project == projects[1]){
            logFC <- logFCtmp
            logPs <- logPstmp
            geneInfo <- gi
        }else{            
            logFC <- merge(logFC, logFCtmp, by="SYMBOL", all=T)
            logPs <- merge(logPs, logPstmp, by="SYMBOL", all=T)
            geneInfo <- rbind(geneInfo, gi[!gi$SYMBOL %in% geneInfo$SYMBOL,])
        }
    }

    ## CREATE FDR ##
    Qs <- make.Qs(logPs)

    ## take log10 of Ps
    logPs[,-1] = -log10(logPs[,-1])

    geneInfo <- geneInfo[match(logFC$SYMBOL, geneInfo$SYMBOL),]
    return(list(logFC = logFC, logPs = logPs, Qs = Qs, GeneInfo = geneInfo))
}


compile.surv.results <- function(project, RsltDir){

    sexes <- c("all","male","female")
    ps <- coef <- list()
  
    print(paste("Working on project", project))
    tmp <- get.surv.rslts(project, RsltDir)
        
    for (sex in sexes){
        
        coef[[sex]] <- tmp[[1]][[sex]][,1:6]
        names(coef[[sex]])[ncol(coef[[sex]])] <- project
        ps[[sex]] <- tmp[[2]][[sex]][,1:6]
        names(ps[[sex]])[ncol(ps[[sex]])] <- project
    }
    
    qs <- ps
    for (i in 1:length(qs)){
        
        qs[[i]][,project] <- make.q(qs[[i]][,project])
    }
    
    return(list(ps, qs, coef))
}

get.surv.rslts <- function(project, RsltDir){
    ps <- list()
    coef <- list()
    for (sex in c("all","male","female")){
        file <- paste0(RsltDir, project,"_lcpm.invnorm.covs_",sex,"_COEF.txt")
        if (!file.exists(file)){
            print(paste("No results for ", project, " : ",sex))
            ps[[sex]] <- NA
            coef[[sex]] <- NA            
        }else{

            coef[[sex]] <- read.table(file, as.is=T, header=T)

            file <- paste0(RsltDir, project,"_lcpm.invnorm.covs_",sex,"_PS.txt")
            ps[[sex]] <- read.table(file, as.is=T, header=T)
        }
    }

    return(list(coef, ps))
}


summary.de.surv.plot <- function(tissues, project, 
                                 ps, qs,
                                 logPs.gtex, logFC.gtex,
                                 logPs.tcga, logFC.tcga,
                                 geneInfo.gtex, geneInfo.tcga,
                                 PlotDir, niter){

    no.x.ind.p <- which(ps[[1]]$chr %in% c(1:22))
    no.x.ind.gtex <- which(geneInfo.gtex$chr %in% c(1:22))
    no.x.ind.tcga <- which(geneInfo.tcga$chr %in% c(1:22))

    de.study <- "_DEgtex_"
    if (sum(tissues %in% projects)){
        de.study = "_DEtcga_"
    }

    plots.per.page <- 3
    tissues.per.plot <- 9
    no.plots <- length(tissues)%/%tissues.per.plot
    no.plots <- no.plots + 1*(length(tissues)%%tissues.per.plot > 0)
    no.blanks <- plots.per.page - no.plots %% plots.per.page
    if (no.blanks == 3){ no.blanks = 0}
    blank.plots <- ggplot(df=data.frame()) + geom_point() + xlim(0, 10) + ylim(0, 100)
    blank.plots <- rep(list(blank.plots), no.blanks)
    file.pre <- paste0(PlotDir, project,de.study,"Surv")

    ResultDir <- gsub("Plots/","", PlotDir)
    ## DE studies are either TCGA or GTEx. Include that in name ##
    DEstudy <- "_GTEx_"
    if (grepl("TCGA",tissues[1])){ DEstudy <- "_TCGA_"}
    out.file <- paste0(ResultDir,project,DEstudy,"DESurv_Best.txt")
    if (file.exists(out.file)){
        cmd <- paste0("mv ",out.file, " ", out.file,".",gsub(" ","",date()))
        system(cmd)
    }
    system(paste0("touch ",out.file))
    file.init = TRUE
    
    for (quant in c(0.001, 0.01)){

        file.pmed <- paste0(file.pre,"_quant",quant,"_PMED.pdf")
        file.pmin <- paste0(file.pre,"_quant",quant,"_PMIN.pdf")        
              
        ps.med <- list()
        ps.min <- list()
        
        plot.title <- paste0(project,", DE_quant=",quant,", All Chrs")
        tmp <- analyze.ranked.survival(tissues, project,
                                       ps, qs,
                                       logPs.gtex, logFC.gtex,  
                                       logPs.tcga, logFC.tcga,
                                       quant=quant,
                                       niter=niter, both=FALSE,
                                       plot.title, no.x=FALSE)      
        ps.med <- c(ps.med, tmp[[1]], blank.plots)
        ps.min <- c(ps.min, tmp[[2]], blank.plots)
        write.table(file = out.file, tmp[[3]][["DEGenes"]],
                    quote=F, row.names=F, col.names=file.init, append=(file.init!=TRUE))
        file.init=FALSE
        
        plot.title <- paste0(project,", DE_quant=",quant, ", No X")
        tmp <- analyze.ranked.survival(tissues, project,
                                       ps, qs,
                                    logPs.gtex[no.x.ind.gtex,], logFC.gtex[no.x.ind.gtex,],  
                                    logPs.tcga[no.x.ind.tcga,], logFC.tcga[no.x.ind.tcga,],
                                    quant=quant,
                                    niter=niter, both=FALSE,
                                    plot.title, no.x=TRUE)
        
        ps.med <- c(ps.med, tmp[[1]], blank.plots)
        ps.min <- c(ps.min, tmp[[2]], blank.plots)
        write.table(file = out.file, tmp[[3]][["DEGenes"]],
                    quote=F, row.names=F, col.names=file.init, append=T)
        
        if (quant == .001 ){ next }

        plot.title <- paste0(project," DE_quant=",quant,
                             ", All Chrs", ", DE in both")
        
        tmp <- analyze.ranked.survival(tissues, project,
                                       ps, qs,
                                       logPs.gtex, logFC.gtex,
                                       logPs.tcga, logFC.tcga,
                                       quant=quant,
                                       niter=niter, both=TRUE,
                                       plot.title, no.x=FALSE)
        
        ps.med <- c(ps.med, tmp[[1]], blank.plots)
        ps.min <- c(ps.min, tmp[[2]], blank.plots)
        write.table(file = out.file, tmp[[3]][["DEGenes"]],
                    quote=F, row.names=F, col.names=file.init, append=T)
        
        ## NO X ##
        plot.title <- paste0(project,", DE_quant=",quant,
                             ", No X", ", DE in both")
        tmp <- analyze.ranked.survival(tissues, project,
                                       ps, qs,
                                       logPs.gtex[no.x.ind.gtex,], logFC.gtex[no.x.ind.gtex,],
                                       logPs.tcga[no.x.ind.tcga,], logFC.tcga[no.x.ind.tcga,],
                                       quant=quant,
                                       niter=niter, both=TRUE,
                                       plot.title, no.x=TRUE)
        
        ps.med <- c(ps.med, tmp[[1]], blank.plots)
        ps.min <- c(ps.min, tmp[[2]], blank.plots)
        write.table(file = out.file, tmp[[3]][["DEGenes"]],
                    quote=F, row.names=F, col.names=file.init, append=T)
        
        ggexport(plotlist = ps.med, filename = file.pmed,
                 nrow = plots.per.page, ncol = 1, width=20, height=8, res=300)
        ggexport(plotlist = ps.min, filename = file.pmin,
                 nrow = plots.per.page, ncol = 1, width=20, height=8, res=300)
        
    }

    print(paste0("Wrote results in file ",out.file))

}
        
  

analyze.ranked.survival <- function(tissues, project, ps, qs,
                                    logPs.gtex, logFC.gtex,
                                    logPs.tcga, logFC.tcga,
                                    quant, niter, both=FALSE,
                                    plot.title, no.x)
    {

        print(paste("Working on analysis using quant =",quant," both= ",both,
                    " Exclude X = ",no.x, 
                    "and ",niter, "iterations to deterine 2.5% survival quantiles"))
        if (both == TRUE){
            print("Requiring gene to be significantly DE in both GTEx and TCGA")
        }

        surv.summary <- summarize.surv.sign(tissues, project,
                                            ps, qs,
                                            logPs.gtex, logFC.gtex,
                                            logPs.tcga, logFC.tcga,
                                            quant=quant,
                                            niter, both=both, no.x)
        mid <- surv.summary[[1]]
        lb <- surv.summary[[2]]
        ub <- surv.summary[[3]]      

        ## CREATE PLOT OF MEDIAN P VALUES AND EXPECTED MEDIAN P VALUES WITH +/- ##
        print("plotting for p.med. . .")
        
        pmed <- plot.tcga.gtex.surv(mid, lb, ub, val="p.med", ylab="Median P-value", plot.title)

        print("plotting for p.min. . .")
        ## REPEAT AGAIN FOR MIN P VALUE ##
        pmin <- plot.tcga.gtex.surv(mid, lb, ub, val="p.min",
                            ylab = "-log10(Minimum P-value)", plot.title)
        

        return(list(pmed, pmin, surv.summary))
    }
        
summarize.surv.sign <- function(tissues, project,
                                ps, qs,
                                logPs.gtex, logFC.gtex,
                                logPs.tcga, logFC.tcga,
                                quant,
                                niter, both, no.x){

    ## TISSUES are the studies from which the genes with DE in the quant
    ## quantile will be identified

    ## PROJECTS are the studies from which the survival statistics will be extracted
    ## from for the DE identified genes

    ## PS are the survival p-value summaries for the PROJECTS studies
    ## QS are the survival q-value summaries for the PROJECTS studies

    ## LOGPS.GTEX are the -log10(p-values) of DE from the TISSUE studies
    ## LOGFS.GTEX are the log(FC) of male versus female expression
    
    ## QUANT is the DE p-value quantile from logPs.gtex from which genes will be
    ## selected

    ## NITER is the number of iterations used to determine the 2.5% tail of the
    ## survival p-value distribution

    ## BOTH is true if the genes have to be DE in both gtex and tcga
    ##

    sexes <- c("all", "male", "female")
    out <- c()
    lb <- lb.alt <- c()
    ub <- ub.alt <- c()

    survPs <- c() ## contain the median, min p-values and their significance
    DEGenes <- c() ## contain the de genes and survival p-values
    
    for (tissue in tissues){               
        
        P.t <- dplyr::select(logPs.tcga, SYMBOL, eval(project)) %>% filter(!is.na(eval(project)))
        FC.t <- dplyr::select(logFC.tcga, SYMBOL, eval(project)) %>% filter(!is.na(eval(project)))
        
        P.g <- dplyr::select(logPs.gtex, SYMBOL, eval(tissue)) %>% filter(!is.na(eval(tissue)))
        FC.g <- dplyr::select(logFC.gtex, SYMBOL, eval(tissue)) %>% filter(!is.na(eval(tissue)))       
        if (sum(tissue == project)){
            logPs <- P.t
            logFC <- FC.t
        }else{
            logPs <- inner_join(P.t, P.g, by="SYMBOL")
            logFC <- inner_join(FC.t, FC.g, by="SYMBOL")
        }
        
        ## NUMBER OF GENES TO EXTRACT ##    
        thresh.val <- ceiling(nrow(logFC) * quant)

        ## INDICES OF GENES SORTED BY DE SIGNIFICANCE
        ind <- order(logPs[,tissue], decreasing=T)
        
        male.ind <- ind[which(ind %in% which(logFC[,tissue] > 0))][1:thresh.val]
        female.ind <- ind[which(ind %in% which(logFC[,tissue] < 0))][1:thresh.val]
        
        
        genes.m <- logPs$SYMBOL[male.ind]
        genes.f <- logPs$SYMBOL[female.ind]
        genes <- unique(c(genes.m, genes.f))

        if (both == TRUE){
                
            ## in tcga and gtex 
                
            ind <- order(logPs[,project], decreasing=T)
            male.ind <- ind[which(ind %in% which(logFC[,project] > 0))][1:thresh.val]
            female.ind <- ind[which(ind %in% which(logFC[,project] < 0))][1:thresh.val]
            
            genes.m <- intersect(logPs$SYMBOL[male.ind], genes.m)
            genes.f <- intersect(logPs$SYMBOL[female.ind], genes.f)
            genes <- unique(c(genes.m, genes.f))            

        }         
            
        DE.male <- data.frame(logPs[male.ind, c("SYMBOL",tissue)], logFC[male.ind,tissue])
        names(DE.male)[2:3] <- c("logP", "logFC")
        DE.female <- data.frame(logPs[female.ind, c("SYMBOL",tissue)], logFC[female.ind,tissue])
        names(DE.female)[2:3] <- c("logP", "logFC")
        names(DE.male)[2:3] <- c("logP", "logFC")

        DEGenes.tmp <- c()
        for (sex in sexes){  ## TCGA survival subset (all, male, female)      
            
            print(paste(project, tissue, sex))
            ## gtex genes in quantile all, male and female upregulated seperately ##
            gene.ind <- which(ps[[sex]]$SYMBOL %in% genes)
            gene.ind.m <- which(ps[[sex]]$SYMBOL %in% genes.m)
            gene.ind.f <- which(ps[[sex]]$SYMBOL %in% genes.f)
            
            P <- filter(ps[[sex]], SYMBOL %in% logPs$SYMBOL) 
            Q <- filter(qs[[sex]], SYMBOL %in% logFC$SYMBOL)
            
            ## CALCULATE NULL VALUES (RANDOMLY CHOOSE THE SAME NUMBER OF GENES ##         
            surv.eith <- surv.summary(P, Q, project, gene.ind)
            surv.m <-  surv.summary(P, Q, project, gene.ind.m)
            surv.f <-  surv.summary(P, Q, project, gene.ind.f)
            
            sim.either <- sim.surv.summary(P, Q, project, size = length(gene.ind),
                                        obs = surv.eith, niter=niter)
            sim.m <-  sim.surv.summary(P, Q, project, size = length(gene.ind.m),
                                       obs = surv.m, niter=niter)
            sim.f <- sim.surv.summary(P, Q, project, size = length(gene.ind.f),
                                      obs = surv.f, niter=niter)         

            out <- rbind(out,
                         c(project, tissue, sex, "either", surv.eith),
                         c(project, tissue, sex, "null.either", sim.either[[1]]),
                         c(project, tissue, sex, "m", surv.m),
                         c(project, tissue, sex, "null.m", sim.m[[1]]),
                         c(project, tissue, sex, "f", surv.f),
                         c(project, tissue, sex, "null.f", sim.f[[1]])
                         )
            
            lb <- rbind(lb, c(project, tissue, sex, "null.either", sim.either[[2]]),
                        c(project, tissue, sex, "null.m", sim.m[[2]]),
                        c(project, tissue, sex, "null.f", sim.f[[2]]))
            
            ub <- rbind(ub, c(project, tissue, sex, "null.either", sim.either[[3]]),
                        c(project, tissue, sex, "null.m", sim.m[[3]]),
                        c(project, tissue, sex, "null.f", sim.f[[3]]))
            
            lb.alt <- rbind(lb.alt, c(project, tissue, sex, "null.either", sim.either[[4]]),
                            c(project, tissue, sex, "null.m", sim.m[[4]]),
                            c(project, tissue, sex, "null.f", sim.f[[4]]))
            
            ub.alt <- rbind(ub.alt, c(project, tissue, sex, "null.either", sim.either[[5]]),
                            c(project, tissue, sex, "null.m", sim.m[[5]]),
                            c(project, tissue, sex, "null.f", sim.f[[5]]))

            ## UPDATE ENRICHMENT SUMMARIES AND P-VALUES ##
            survPs <- rbind(survPs, c(project, tissue, sex, "either",  surv.eith, sim.either[[6]]),
                            c(project, tissue, sex, "male",  surv.m, sim.m[[6]]),
                            c(project, tissue, sex, "female",  surv.f, sim.f[[6]]))

            ## UPDATE SURVIVAL P-VALUES AND ENRICHMENT P-VALUES OF GENES
            ## DE IN EITHER ##
            if (length(gene.ind) > 0){
                DEGenes.tmp <- rbind(DEGenes.tmp,
                                     cbind(project, tissue, sex, desex="either", quant=quant,
                                           excl.x = no.x, both=both, 
                                           ps[[sex]][gene.ind,],
                                           pmn=surv.eith[1], pmd=surv.eith[2], pmin=surv.eith[3],
                                           emppmn=sim.either[[6]][1], emppmd=sim.either[[6]][2],
                                           emppmin=sim.either[[6]][3]))
            }
            ## DE IN MALES ##
            if (length(gene.ind.m) > 0){
                DEGenes.tmp <- rbind(DEGenes.tmp,
                                     cbind(project, tissue, sex, desex="male", quant=quant,
                                           excl.x = no.x, both=both,
                                           ps[[sex]][gene.ind.m,],
                                           pmn=surv.m[1], pmd=surv.m[2], pmin=surv.m[3],
                                           emppmn=sim.m[[6]][1], emppmd=sim.m[[6]][2],
                                           emppmin=sim.m[[6]][3]))
            }

            ## DE IN FEMALES
            if (length(gene.ind.f) > 0){
                DEGenes.tmp <- rbind(DEGenes.tmp,
                                     cbind(project, tissue, sex, desex="female", quant=quant,
                                           excl.x = no.x, both=both,
                                           ps[[sex]][gene.ind.f,],
                                           pmn=surv.f[1], pmd=surv.f[2], pmin=surv.f[3],
                                           emppmn=sim.f[[6]][1], emppmd=sim.f[[6]][2],
                                           emppmin=sim.f[[6]][3]))
            }

           
        }
        if (length(DEGenes.tmp) > 0){
            ## MERGE IN DE P-VALUES AND FC
            lp <- logPs[,c("SYMBOL", tissue)]
            names(lp) <- gsub(tissue,"logP", names(lp))
            DEGenes.tmp <- merge(DEGenes.tmp, lp, by = "SYMBOL", all.x=T, sort=F)
            
            fc <- logFC[,c("SYMBOL", tissue)]
            names(fc) <- gsub(tissue, "logFC", names(fc))
            DEGenes.tmp <- merge(DEGenes.tmp, fc, by = "SYMBOL", all.x=T, sort=F)
            
        }
        
        if (length(DEGenes.tmp) > 0){
            names(DEGenes.tmp)[grep(project, names(DEGenes.tmp))] <- "Psurv"
            DEGenes <- rbind(DEGenes, DEGenes.tmp)
        }
    }
        
    out <- as.data.frame(out)
    lb <- as.data.frame(lb)
    ub <- as.data.frame(ub)
    lb.alt <- as.data.frame(lb.alt)
    ub.alt <- as.data.frame(lb.alt)
    names(out) <- names(lb) <- names(ub) <- names(lb.alt) <- names(ub.alt) <-
        c("project", "gtex.tis", "tcga.set",
          "gtex.upreg.in","p.mean",
          "p.med","p.min", "q.mean", "q.med", "q.min")


    return(list(out, lb, ub, lb.alt, ub.alt, DEGenes=DEGenes))
}

surv.summary <- function(P, Q, project, gene.ind){
    
    if (length(gene.ind) == 0){
        out <- rep(NA, 6)
    }else{
        
        P <- P[gene.ind , project]
        if (sum(!is.na(P)) > 0){
            Q <- Q[gene.ind , project]
            mean.p <- mean(P, na.rm=T)
            med.p <- median(P, na.rm=T)
            min.p <- min(P, na.rm=T)
            
            mean.q <- mean(Q, na.rm=T)
            med.q <- median(Q, na.rm=T)
            min.q <- min(Q, na.rm=T)
            
            out <- c(mean.p, med.p, min.p, mean.q, med.q, min.q)      
        }else{

            out <- rep(NA, 6)
        }
    }            
    return(out)
}

sim.surv.summary <- function(P, Q, project, size, obs, niter, ci=0.025){

    inds <- which(!is.na(P[,project]))
    mid <- lb <- ub <- c()
    if (size == 0){
        mid <- rep(NA, 6)
        lb <- lb.alt <- rep(NA, 6)
        ub <- ub.alt <- rep(NA, 6)
        ps.emp <- rep(NA,6)

    }else{

        store <- c()
        f <- function(){ ind <- sample(inds, size=size, replace=F);
                         surv.summary(P, Q, project, ind)}
        store <- replicate(niter, f())           

        lb <- apply(store, 1, quantile, probs = ci) 
        ub <- apply(store, 1, quantile, probs = 1 - ci)
        mid <- rowMeans(store)

        sds <- sqrt(apply(store, 1, var))
        adj <- qnorm(1-ci)
        lb.alt <- mid - sds*adj
        ub.alt <- mid + sds*adj

        ps.emp <- (rowSums(apply(store, 2, function(x) x < obs)) + 1) / (rowSums(!is.na(store)) + 1)
        
    }
    return(list(mid, lb, ub, lb.alt, ub.alt, ps.emp))
}
                                           


plot.tcga.gtex.surv <- function(mid, lb, ub, val="p.med",
                                ylab = "Median P-value", plot.title){

    tcga.labs <- c(all = "Both", male = "M", female= "F")

    gtex.labs <- c(Adipose_Subcutaneous = " AdSub",
                   Adipose_Visceral_Omentum = "AdVsc",
                   Adrenal_Gland = "AGld",
                   Artery_Aorta = "Aort",
                   Artery_Coronary = "ArtCor",
                   Artery_Tibial = "ArtTib",
                   Brain_Amygdala = "BrAmy",
                   Brain_Anterior_cingulate_cortex_BA24 = "BrAnt",
                   Brain_Caudate_basal_ganglia = "BrCaud",
                   Brain_Cerebellar_Hemisphere = "BrCerH",
                   Brain_Cerebellum = "BrCer",
                   Brain_Cortex = "BrCtx",
                   Brain_Frontal_Cortex_BA9 = "BrFCtx",
                   Brain_Hippocampu = "BrHip",
                   Brain_Hypothalamus = "BrHyp",
                   Brain_Nucleus_accumbens_basal_ganglia = "BrNuc",
                   Brain_Putamen_basal_ganglia = "BrPut",
                   Brain_Spinal_cord_cervical_c_1 = "BrSp",
                   Breast_Mammary_Tissue = "Brst",
                   Cells_Cultured_fibroblasts = "Fblst",
                   Cells_EBV_transformed_lymphocytes = "Lcl",
                   Colon_Sigmoid = "ColSg",
                   Colon_Transverse = "ColTv",
                   Esophagus_Gastroesophageal_Junction = "EsGJt",
                   Esophagus_Mucosa = "EsMuc",
                   Esophagus_Muscularis = "EsMsc",
                   Heart_Atrial_Appendage = "HtAt",
                   Heart_Left_Ventricle = "HtLV",
                   Liver = "Liv",
                   Lung = "Lng",
                   Minor_Salivary_Gland = "SalG",
                   Muscle_Skeletal = "MSk",
                   Nerve_Tibial = "NvTb",
                   Pancreas = "Panc",
                   Pituitary = "Pit",                   
                   Skin_Not_Sun_Exposed_Suprapubic = "SkNE",
                   Skin_Sun_Exposed_Lower_leg = "SkE",
                   Stomach = "Stm",
                   Thyroid = "Thr",
                   Whole_Blood = "Bld",                  
                   TCGA_GBM = "GBM", TCGA_SKCM = "SKCM",  TCGA_THCA = "THCA",
                   TCGA_LIHC = "LIHC", TCGA_LUAD = "LUAD", TCGA_KIRC="KIRC",
                   TCGA_LAML = "LAML", TCGA_HNSC = "HNSC", TCGA_LGG = "LGG",
                   TCGA_LUSC = "LUSC", TCGA_STAD="STAD", TCGA_COAD="COAD",
                   TCGA_BLCA="BLCA", TCGA_KIRP="KIRP", TCGA_SARC="SARC",
                   TCGA_PAAD="PAAD", TCGA_ESCA="ESCA", TCGA_PCPG="PCPG",
                   TCGA_READ="READ", TCGA_THYM="THYM", TCGA_KIRH="KIRH",
                   TCGA_ACC="ACC", TCGA_MESO="MESO", TCGA_UVM="UVM",
                   TCGA_DLBC="DLBC", TCGA_CHOL="CHOL")
    
    d <- as.tibble(mid)
    d <- dplyr::select(d, project, gtex.tis, tcga.set, gtex.upreg.in, eval(val))

    ## d <- dplyr::rename(d, "VAL"=eval(val))
    d$VAL <- unlist(d[, which(names(d) == val)])
    d <- mutate(d, VAL = as.numeric(as.character(VAL)))
    d$null <- grepl("null", d$gtex.upreg.in)
    d$gtex.upreg.in <- gsub("null\\.","",d$gtex.upreg.in)

    lb <- as.tibble(lb)
    lb <- dplyr::select(lb, project, gtex.tis, tcga.set, gtex.upreg.in, eval(val))
    ## lb <- dplyr::rename(lb, "lb" = val)
    lb$lb <- unlist(lb[,val])
    lb <- mutate(lb, lb = as.numeric(as.character(lb)))
    lb$gtex.upreg.in <- gsub("null\\.","",lb$gtex.upreg.in)
    
    ub <- as.tibble(ub)
    ub <- dplyr::select(ub, project, gtex.tis, tcga.set, gtex.upreg.in, eval(val))
    ## ub <- dplyr::rename(ub, "ub" = eval(val))
    ub$ub <- unlist(ub[,val])
    ub <- mutate(ub, ub = as.numeric(as.character(ub)))
    ub$gtex.upreg.in <- gsub("null\\.","",ub$gtex.upreg.in)

    if (length(grep("min", val))==1){
        d <- mutate(d, VAL = -log10(VAL))
        lb <- mutate(lb, lb = -log10(lb))
        ub <- mutate(ub, ub = -log10(ub))
    }
    
    d <- full_join(d, lb, by=c("project","gtex.tis","tcga.set","gtex.upreg.in"))
    d <- full_join(d, ub, by=c("project","gtex.tis","tcga.set","gtex.upreg.in"))
    d$gtex.tis <- gsub("\\-", "_", d$gtex.tis)
    ind <- which(d$null==FALSE)
    d$lb[ind] <- d$ub[ind] <- d$VAL[ind]

    tissues <- unique(d$gtex.tis)
    no.tis <- length(unique(d$gtex.tis))
    tis.p.plt = 9
    nplots<- no.tis %/% tis.p.plt
  
    if (tis.p.plt*nplots < no.tis){
        nplots <- nplots+1
    }

    p <- list()
    for (i in 1:nplots){

        ind <- intersect(c(1:9) + (i-1)*tis.p.plt, 1:no.tis)
        tss <- gsub("-","_",tissues[ind])
        tmp <- d %>% filter(project == project & gtex.tis %in% tss)
    
        p[[i]] <- ggplot(tmp, aes(x=gtex.upreg.in, y=VAL,color=null)) + 
            geom_pointrange(aes(ymin=lb, ymax=ub), size = .3, fatten=.5) +
                facet_grid(.~gtex.tis + tcga.set, labeller = 
                           labeller(tcga.set=tcga.labs, gtex.tis = gtex.labs)) +
                               theme(strip.text.x=element_text(angle=90, size=4)) +
                               ggtitle(plot.title) +
                                   ylab(ylab) + xlab("Sex of upregulated GTEx genes") +
                                       theme(axis.text.x=element_text(angle = 90, 
                                                 hjust = 0, size=6), legend.position="none") 
                                ## scale_color_discrete(name="", labels=c("Obsv", "Null"))
    }
    
    return(p)    

}
        

make.Qs <- function(logPs){
    Qs <- logPs
    for (i in 2:ncol(logPs)){
        p <- logPs[,i]
        miss.ind <- which(is.na(p))
        if (length(miss.ind) > 0){
            p <- p[-miss.ind]
        }
        small.ind<- which(p < 10^-30)
        p[small.ind] <- 10^-30
        q <- qvalue(p)$qvalues
        if (length(miss.ind) > 0){
            Qs[-miss.ind,i] <- q
        }else{
            Qs[,i] <- q
        }
    }
    return(Qs)
}
