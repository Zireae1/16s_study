rm(list=ls(all=TRUE)); gc()
################################################
# 16S rRNA metagenomes (Base Rep)
################################################

#--- choose your flavor: ---
current_directory <- "~/metagenome/AVVA"
#current_directory <- ""
setwd(current_directory)

source("lib/functions.R")
source("lib/visualization.R")

#loading libraries
library(stringr)
library(ecodist)
library(stringr)
library(ecodist)
library(data.table)
library(ape)
library(ggplot2)
library(scales)
library(MASS)
library(stringr)
library(matrixStats)
library(gridExtra)
library(grid)
library(phyloseq) #source("https://bioconductor.org/biocLite.R"), biocLite("phyloseq")
library(reshape)
library(gplots)
require(MASS)
library(futile.logger)


flog.info("start")

path.to.config <- "/home/anna/metagenome/AVVA/src/pathway.ini"
pathway <- ReadIni(path.to.config) 
totalTable <- LoadSixProject(pathway)

family.avva <- totalTable$family[which(rownames(totalTable$family) 
                                       %in% totalTable$meta[which(totalTable$meta[,"Type.1"] 
                                                                  %in% "EKOK_1"),"samples_name"]),]
genus.avva <- totalTable$genus[which(rownames(totalTable$family) 
                                      %in% totalTable$meta[which(totalTable$meta[,"Type.1"] 
                                                                 %in% "EKOK_1"),"samples_name"]),]
species.avva <- totalTable$species[which(rownames(totalTable$family) 
                                     %in% totalTable$meta[which(totalTable$meta[,"Type.1"] 
                                                                %in% "EKOK_1"),"samples_name"]),]

message("start working with alpha diversity")
alphaDivCase <- (totalTable$AlphaDiv[which(rownames(totalTable$AlphaDiv) 
                                           %in% totalTable$Meta[which(totalTable$Meta[,"Type.1"] 
                                                                      %in% "case"),"samples_name"]),])]

WriteTable(alphaDivCase, pathway$Case$OutdirCase, "AlphaDivCaseTbl"))
AlphaDivMean <- lapply(AlphaDivCaseL[match('AlphaDiversity', names(AlphaDivCaseL))], mean)
AlphaDivSd <- lapply(AlphaDivCaseL[match('AlphaDiversity', names(AlphaDivCaseL))], sd)
AlphaDivMeanSd <- c(AlphaDivMean, AlphaDivSd)
names(AlphaDivMeanSd) <- c("mean", "sd")

WriteTable (AlphaDivMeanSd, OutdirCase, "AlphaDivMeanAndSd") 

cairo_pdf((paste(OutdirCase, '/Graphs/alpha_boxplot.pdf', sep = "/")), width = 10,  height = 10)
boxplot(alphaDivCase$AlphaDiversity)
dev.off()

message("end start working with alpha diversity")

message("start sequencing statistics")
#ToDo: statistics for project samples
seqStatTable <- read.table(pathway$statistics$SeqStatTblFile, comment.char = "", quote ="", as.is = T, header = T)

seqStatList <- as.list(seqStatTable)
seqStatMean <- lapply(seqStatList[-match('samp_name', names(seqStatList))], mean)
seqStatSd <- lapply(seqStatList[-match('samp_name', names(seqStatList))], sd)
names(seqStatMean) <- c("mean_init", "mean_AF", "mean_AR", "mean_Mapp")
names(seqStatSd) <- c("sd_init", "sd_AF", "sd_AR", "sd_Mapp")
seqStatMeanSd <- c(seqStatMean, seqStatSd)
seqStatMeanSd.df <- t(as.data.frame(seqStatMeanSd))

WriteTable(seqStatMeanSd.df, pathway$Case$OutdirCase, "SeqStatMeanAndSd")
message("end sequencing statistics")
###########################################################################3 test code
######################################################################################

#ToDO: запись таблиц топовых значений для неограниченного числа групп по family, genus, species, otu
flog.info("writing")

fam.gen.spe.list <- list(family.avva, genus.avva, species.avva)
TopFeatures <- function(FeatList)
{
  TopFeatures <- lapply(FeatList, chooseTOPfeature, perc = 85)
  
}

fam.gen.spe.listTOP <- TopFeatures(fam.gen.spe.list)
#### Print Mean and Stdev for group of samples 
PrintMeanAndStd(fam.gen.spe.listTOP, min.trsh=1)

#--- draw MDS for different taxonomy levels---
for (i in names(totalTable)[1:3]){
  dist<-bcdist(totalTable[[i]][which(totalTable$meta[,"Type.1"] %in% "EKOK_1"),])
  MDS(dist, pathway$group.one$GraphsDir, totalTable$meta[which(totalTable$meta[,"Type.1"] %in% "EKOK_1"),],
      "Type.1", i)
}

# выдает имена образцов из таблицы метаданных которым соответствует EKOK_1 
#totalTable$meta[which(totalTable$meta[, "Type.1"] %in% "EKOK_1"), "samples_name"]
