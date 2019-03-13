ScriptDir <- "/gpfs/data/stranger-lab/askol/TCGA2/DE_Surv/Scripts/"
CodeDir <- "/gpfs/data/stranger-lab/askol/TCGA2/DE_Surv/Code/"

projects = read.table(file = "/gpfs/data/stranger-lab/askol/TCGA/TCGA_Target_projects.txt",
    header=TRUE, as.is=TRUE, sep="\t")
projects = projects[grep("TCGA",projects[,1]),1]

skip.cancers <- c("TCGA-BRCA","TCGA-OV", "TCGA-CESC", "TCGA-PRAD", "TCGA-TGCT",
                  "TCGA-UCS", "TCGA-UCEC")
projects <- projects[!projects %in% skip.cancers] 

create.pbs <- function(ScriptDir, project){
    
    file <- paste0(ScriptDir, project, "_summary.pbs")
    sh.txt <- rbind(
        "#!/bin/bash",
        "#PBS -l  nodes=1:ppn=1,mem=8gb",
        "#PBS -l walltime=96:00:00",
        paste0("#PBS -o ", ScriptDir),
        "#PBS -j oe",
        paste0("#PBS -N ",gsub("TCGA-","",project)),
        "module load R")
    
    write.table(file = file, sh.txt, quote=F, row.names=F, col.names=F)

    return(file=file)
}

job.files <- c()
for (project in projects){


    cmd <-  paste0("R CMD BATCH  '--args ",project,
                   "' ", CodeDir, 
                   "DE_Surv_Plots.r ",
                   ScriptDir,
                   project,".out")      
        
    file <- create.pbs(ScriptDir, project)
    write.table(file = file, cmd, quote=F, row.names=F, col.names=F, append=T)

    job.files <- c(job.files, file)
    
}

for (job in job.files){
    system(paste("qsub",job))
}
