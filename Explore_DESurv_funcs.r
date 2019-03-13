collect.data <- function(ResultDir){

    patrn <- "TCGA_DESurv_Best\\.txt$|GTEx_DESurv_Best\\.txt$"
    files <- dir(ResultDir, pattern=patrn)
    data <- c()
    
    for (file in files){

        if (file.exists(paste0(ResultDir,file)) == FALSE){ next }
        print(paste0("Working on file: ",file))
        data.tmp <- read.table(file = paste0(ResultDir,file), header=T, as.is=T)        
        data.int <- filter(data.tmp, emppmn<.05 | emppmd<.05 | emppmin<.05) %>%
            filter(Psurv < 0.05)

        data <- rbind(data, data.int)

    }

    return(data)
}

find.shared <- function(data, geneInfo){
    
    tcga <- unique(c(data$tissue, data$project))
    tcga <- tcga[grep("TCGA", tcga)]
    gtex <- unique(data$tissue)
    gtex <- gtex[grepl("TCGA", gtex) == FALSE]
    
    ## DE DEFINED BY TISSUES
    ##   detiss <- tcga
    ##   if (DEtissue == "GTEx"){ detiss <- gtex }
    
    ## DE DEFINED BY TCGA CANCERS
    ##   data <- filter(data, tissue %in% detiss)
    
    ## FIND GENES THAT OCCUR MORE THAN ONCE ##
    genesDupe <- unique(data$SYMBOL[duplicated(data$SYMBOL)])
    genesDupe <- filter(data, SYMBOL %in% genesDupe)
    genesDupe <- dplyr::select(genesDupe, SYMBOL, project, tissue, Psurv) %>% 
        dcast(SYMBOL + project ~ tissue, mean,  value.var = "Psurv")
    Psurv <- rowMeans(genesDupe[,-c(1:2)], na.rm=T)
    NoDETiss <- rowSums(!is.na(genesDupe[,-c(1:2)]))
    genesDupe$Psurv <- Psurv; genesDupe$NoDETiss <- NoDETiss
    genesDupe <- dcast(genesDupe, SYMBOL ~ project, mean, value.var = "Psurv")
    
    genesDupe$NoStd <- rowSums(!is.na(genesDupe[,which(names(genesDupe) %in% tcga)]), na.rm=T)
    ind <- which(names(genesDupe) %in% tcga)   
    genesDupe$MedP <- apply(genesDupe[,ind], 1, median, na.rm=T)
    genesDupe$MeanP <- apply(genesDupe[,ind], 1, mean, na.rm=T)
    genesDupe$minP <- apply(genesDupe[,ind], 1, min, na.rm=T)
    genesDupe$maxP <- apply(genesDupe[,ind], 1, max, na.rm=T)
    
    ## GET NUMBER OF DE TISSUES FOR EACH GENE
    DE <-  unique(data$SYMBOL[duplicated(data$SYMBOL)])
    DE <- filter(data, SYMBOL %in% DE)
    DE  <- dplyr::select(DE, SYMBOL, tissue, Psurv) %>% 
        dcast(SYMBOL ~ tissue, length,  value.var = "Psurv")
    DE[,-1] <- 1*(DE[,-1] >0)
    names(DE)[-1] <- paste0("DE_",names(DE)[-1])
    DE$NoDE <- rowSums(DE[,-1])
    
    genesDupe <- left_join(genesDupe, DE, by="SYMBOL")
    genesDupe <- arrange(genesDupe, desc(NoStd), desc(NoDE))
    genesDupe <- genesDupe %>% filter(NoStd > 1)
    
    ## add annotation
    genesDupe <- right_join(geneInfo[,c("symbol","chr","start","biotype", "GeneFreq")],
                            genesDupe, by=c("symbol" = "SYMBOL"))
    
    return(genesDupe)
}


join.shared <- function(data, qnt, psurv, empmd, exX, DETiss){
    
    ## DETERMINE WHICH GENES ARE PROGNOSTIC ACROSS CANCERS ##
    ## all/either male/male male/female female/male female/female
    
    j <- c()
    
    for (srvsex in c("all","male","female")){
        for (desex in c("either","male","female")){
            
            out <- filter.data(data, sx=srvsex, desx=desex, qnt=qnt,
                               exX=exX, psurv=psurv, empmd=empmd,
                               DETiss=DETiss)
            out <- find.shared(out, geneInfo)
            out <-  melt(out[,c("symbol",names(out)[grep("^TCGA", names(out))])],
                         value.name="Psurv", variable.name="project", id.vars="symbol")
            out <- filter(out, !is.na(Psurv))
            out$smp.de <- paste0(substring(srvsex,1,1), substring(desex,1,1))
            j <- rbind(j, out)
        }
    }
    
    return(j)
}

filter.data <- function(data, sx, desx, qnt, psurv, empmd, exX, DETiss) {
    return(
        data %>% filter(sex==sx, desex==desx, quant == qnt,
                        excl.x==exX, Psurv <= psurv,
                        emppmd <= empmd, DEPop == DETiss)
        )
}

plot.sharing.all <- function(data, PlotDir){    
    
    ## DE BY TCGA : NO X CHROMOSOME GENES
    print("Working on DE by TCGA, No X")
    shared.join <- join.shared(data, qnt=0.01, psurv=0.01, empmd=0.01,
                               exX=TRUE, DETiss = "TCGA")
    PlotFile <- paste0(PlotDir, "Shared_DeTCGA_NoX_Psurv.pdf")
    plot.sharing(shared.join, PlotFile)
    
    ## DE BY TCGA : X CHROMOSOME GENES INCL
    print("Working on DE by TCGA, X genes included")
    shared.join <- join.shared(data, qnt=0.01, psurv=0.01, empmd=0.01,
                               exX=FALSE, DETiss = "TCGA")
    PlotFile <- paste0(PlotDir, "Shared_DeTCGA_InclX_Psurv.pdf")
    plot.sharing(shared.join, PlotFile)
    
    ## DE BY GTEX : NO X
    print("Working on DE by GTEx, no X genes")
    shared.join <- join.shared(data, qnt=0.01, psurv=0.01, empmd=0.01,
                               exX=TRUE, DETiss= "GTEx")
    PlotFile <- paste0(PlotDir, "Shared_DeGTEx_NoX_Psurv.pdf")
    plot.sharing(shared.join, PlotFile)
    
    ## DE BY GTEX : X INCL
    print("Working on DE by GTEx, X genes included")
    shared.join <- join.shared(data, qnt=0.01, psurv=0.01, empmd=0.01,
                               exX=FALSE, DETiss = "GTEx")
    PlotFile <- paste0(PlotDir, "Shared_DeGTEx_InclX_Psurv.pdf")
    plot.sharing(shared.join, PlotFile)
    
}    



plot.sharing <- function(shared.join, PlotFile){
    
    file <- PlotFile
    
    shared.join <- mutate(shared.join, SPop = substr(smp.de,1,1))
    shared.join <- mutate(shared.join, DEPop = substr(smp.de,2,2))
    POP <- list(a = "All", m = "Male", f= "Female")

    pdf(file = file, width = 18, height=8)
    
    for (sp in c("a","m","f")){

        
       
        j <- filter(shared.join, grepl("DE", project)==FALSE)
        j <- filter(j, !is.na(Psurv))
        
        genes <- filter(j, SPop == sp) %>% select(symbol)
        genes <- unlist(unique(genes))
        
        j <- filter(j, symbol %in% genes)
                
        j <- mutate(j, proj = paste0(project, ".", SPop)) %>% mutate(proj = gsub("TCGA-","",proj))
        ## CREATE ID TO IDENTIFY GENE-PROJ-SPOP GROUPS ##
        j <- mutate(j, id = paste(symbol, project, SPop, sep="."))
        dupes.ind <- j$id %in% j$id[duplicated(j$id)]
        
        ## REMOVE EITHERS IF GENE IS ALSO IN MALE OR FEMALE ##
        rm.ind <- which(dupes.ind & j$DEPop == "e")
        j <- j[-rm.ind,]       
        
        j<- mutate(j, logP = -log10(Psurv))
        t <- rowSums(table(j$symbol, j$proj))
        gene.ord <- names(sort(t, decreasing=F))
        j$symbol <- factor(j$symbol, levels=gene.ord)
        lvl <- c(sapply(gsub("TCGA-", "", names(sort(table(j$project), decreasing=T))),
                        paste, c("a","m","f"), sep="."))
        j$proj <- factor(j$proj, levels = lvl)
        p <- ggplot(j, aes(x = symbol, y = proj)) +  xlab("Gene") +
                theme(
                    ##remove plot background
                    plot.background=element_blank(),
                    ##remove plot border
                    panel.border=element_blank(),
                    axis.text.y = element_text(angle = 0, hjust = 1, size=8),
                    axis.text.x = element_text(angle = 90, hjust = 1, size=6))
        
        p <-  p + geom_tile(aes(fill=logP), data=subset(j, DEPop=="e")) + 
            scale_fill_gradient("DE Either Only", high = "green", low="white", limits=c(0,10))
        
        p <- p + new_scale("fill")
        p <-  p + geom_tile(aes(fill=logP), data=subset(j, DEPop=="m")) +
            scale_fill_gradient2("UpReg in Males", high = "blue", low="white", limits=c(0,10))
        
        p <- p + new_scale("fill")
        p <- p +geom_tile(aes(fill=logP), data=subset(j, DEPop=="f")) +
            scale_fill_gradient2("UpReg in Females", high = "red", low="white", limits=c(0,10))

        ggtitle(paste0("Genes with Significant Survival P-Values in ", POP[[sp]]))
        print(p)
        
    }
    dev.off()
    print(paste0("Wrote to plots to ",file))
}
