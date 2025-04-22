
###Transcriptional noise quantification using BASiCS###


library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(BASiCS)


##Load and prepare data##


#Load Seurat or SingleCellExperiment object containing pre-processed scRNA-seq data.
seu<-readRDS("Roux_SKNSH.rds")

#If Seurat loaded, conver to SingleCellExperiment.
sce<-as.SingleCellExperiment(seu)

#Separate different experiments into individual objects (e.g. Untreated, Cisplatin course).
sce_UT <- subset (sce, ,sce$Condition == c("Untreated"))
sce_Cis <- subset (sce, ,sce$Condition == c("Untreated", "Cisplatin(1)_ON", "Cisplatin(2)_ON", "Cisplatin_1weeksOFF", "Cisplatin_4weeksOFF"))

#From here showing example pipeline for sce_Cis for simplicity.

#Add replicate information under "BatchInfo" 
sce_Cis$BatchInfo<- recode(as.character(sce_Cis$description),
                           "UT A" = "1",
                           "UT B" = "2", 
                           "UT C" = "3",
                           "Cisplatin_1weeksOFF A" = "1",
                           "Cisplatin_1weeksOFF B" = "2", 
                           "Cisplatin_1weeksOFF C" = "3",
                           "Cisplatin_4weeksOFF A" = "1",
                           "Cisplatin_4weeksOFF B" = "2", 
                           "Cisplatin_4weeksOFF C" = "3",
                           "Cisplatin(1)_ON A" = "1",
                           "Cisplatin(1)_ON B" = "2", 
                           "Cisplatin(1)_ON C" = "3",
                           "Cisplatin(2)_ON A" = "1",
                           "Cisplatin(2)_ON B" = "2", 
                           "Cisplatin(2)_ON C" = "3",
)

#Eliminate genes that have 0 expression (not expressed in any cells). 
sce_Cis <- sce_Cis[rowSums(logcounts(sce_Cis)) != 0, ]

#Separate the groups to be compared in individual objects. For intrinsic plasticity analysis this would be ADRN, BRI, MES cells in the untreated condition. For cisplatin course, this would be the MES cells at the different timepoints.
sce.ut<- subset (sce_Cis, ,sce_Cis$Condition == "Untreated" & sce_Cis$AMT.state == "MES")
sce.cis<- subset (sce_Cis, ,sce_Cis$Condition == "Cisplatin(1)_ON" & sce_Cis$AMT.state == "MES")
sce.r1<- subset (sce_Cis, ,sce_Cis$Condition == "Cisplatin_1weeksOFF" & sce_Cis$AMT.state == "MES")
sce.r4<- subset (sce_Cis, ,sce_Cis$Condition == "Cisplatin_4weeksOFF" & sce_Cis$AMT.state == "MES")

#Randomly subset the same number of cells from each object (use the maximum number of cells in the smallest object)
sce.ut <- sce.ut[ ,sample(ncol(sce.ut), 130)] 
sce.cis <- sce.cis[ ,sample(ncol(sce.cis), 130)] 
sce.r1 <- sce.r1[ ,sample(ncol(sce.r1), 130)] 
sce.r4 <- sce.r4[ ,sample(ncol(sce.r4), 130)] 


##Parameter estimation##

#Create the BASiCS Markov Chains
set.seed(100)
chain.ut= BASiCS_MCMC(sce.UT, N = 4000, Thin = 10, Burn = 2000, WithSpikes = FALSE, Regression = F, PrintProgress = TRUE, StoreChains = TRUE, StoreDir = tempdir(), RunName = "UT")
chain.cis = BASiCS_MCMC(sce.cis, N = 4000, Thin = 10, Burn = 2000, WithSpikes = FALSE, Regression = F, PrintProgress = TRUE, StoreChains = TRUE, StoreDir = tempdir(), RunName = "Cisplatin")
chain.r1 = BASiCS_MCMC(sce.r1, N = 4000, Thin = 10, Burn = 2000, WithSpikes = FALSE, Regression = F, PrintProgress = TRUE, StoreChains = TRUE, StoreDir = tempdir(), RunName = "1weekOFF")
chain.r4 = BASiCS_MCMC(sce.r4, N = 4000, Thin = 10, Burn = 2000, WithSpikes = FALSE, Regression = F, PrintProgress = TRUE, StoreChains = TRUE, StoreDir = tempdir(), RunName = "4weekOFF")

#Saved chains can be reloaded if needed
chain.ut<- BASiCS_LoadChain("UT", StoreDir = "")
chain.cis<- BASiCS_LoadChain("Cisplatin", StoreDir = "")
chain.r1<- BASiCS_LoadChain("1weekOFF", StoreDir = "")
chain.r4<- BASiCS_LoadChain("4weekOFF", StoreDir = "")

#Visually assess chain convergence for the parameters of interest: mean (mu) and overdispersion (delta)
plot(chain.ut, Param = "mu", Gene = 1, log = "y")
plot(chain.ut, Param = "delta", Gene = 1, log = "y")
plot(chain.cis, Param = "mu", Gene = 1, log = "y")
plot(chain.cis, Param = "delta", Gene = 1, log = "y")
plot(chain.r1, Param = "mu", Gene = 1, log = "y")
plot(chain.r1, Param = "delta", Gene = 1, log = "y")
plot(chain.r4, Param = "mu", Gene = 1, log = "y")
plot(chain.r4, Param = "delta", Gene = 1, log = "y")


##Compare transcriptional noise levels between groups##

#For differential expression and varialbility analysis use BASiCS_TestDE
de_uvc <- BASiCS_TestDE(Chain1 = chain.ut, Chain2 = chain.cis, GroupLabel1 = "Untreted", GroupLabel2 = "Cisplatin", EpsilonM = log2(2), EpsilonD = log2(1.5), EpsilonR = 0.41, ProbThresholdM = 0.85, Plot = TRUE)
de_uvr1 <- BASiCS_TestDE(Chain1 = chain.ut, Chain2 = chain.r1, GroupLabel1 = "Untreted", GroupLabel2 = "1-week recovery", EpsilonM = log2(2), EpsilonD = log2(1.5), EpsilonR = 0.41, ProbThresholdM = 0.85, Plot = TRUE)
de_uvr4 <- BASiCS_TestDE(Chain1 = chain.ut, Chain2 = chain.r4, GroupLabel1 = "Untreated", GroupLabel2 = "4-week recovery", EpsilonM = log2(2), EpsilonD = log2(1.5), EpsilonR = 0.41, ProbThresholdM = 0.85, Plot = TRUE)
de_cvr1 <- BASiCS_TestDE(Chain1 = chain.cis, Chain2 = chain.r1, GroupLabel1 = "Cisplatin", GroupLabel2 = "1-week recovery", EpsilonM = log2(2), EpsilonD = log2(1.5), EpsilonR = 0.41, ProbThresholdM = 0.85, Plot = TRUE)
de_cvr4 <- BASiCS_TestDE(Chain1 = chain.cis, Chain2 = chain.r4, GroupLabel1 = "Cisplatin", GroupLabel2 = "4-week recovery", EpsilonM = log2(2), EpsilonD = log2(1.5), EpsilonR = 0.41, ProbThresholdM = 0.85, Plot = TRUE)

de_uvc
de_uvr1
de_uvr4
de_cvr1
de_cvr4

#Differentially noisy genes for each comparison can be saved in a dataframe
summary_uvc <- as.data.frame(de_uvc, Parameter = "Disp")
summary_uvr1 <- as.data.frame(de_uvr1, Parameter = "Disp")
summary_uvr4 <- as.data.frame(de_uvr4, Parameter = "Disp")
summary_cvr1<- as.data.frame(de_cvr1, Parameter = "Disp")
summary_cvr4<- as.data.frame(de_cvr4, Parameter = "Disp")

#Volcano plots of overdispersion (noise)
uvc <- BASiCS_PlotDE(de_uvc, Plots="Volcano", Parameters="Disp") #Remove Parameters if you also want to plot changes in mean
uvr1 <- BASiCS_PlotDE(de_uvr1, Plots="Volcano", Parameters="Disp")
uvr4 <- BASiCS_PlotDE(de_uvr4, Plots="Volcano", Parameters="Disp")
cvr1 <- BASiCS_PlotDE(de_cvr1, Plots="Volcano", Parameters="Disp")
cvr4 <- BASiCS_PlotDE(de_cvr4, Plots="Volcano", Parameters="Disp")


##Identify highly variable genes within groups##

#To detect highly variable genes (HVG) in individual conditions use BASiCS_DetectHVG
HVGut <- BASiCS_DetectHVG(chain.ut, VarThreshold = 0.6)
HVGcis <- BASiCS_DetectHVG(chain.cis, VarThreshold = 0.6)
HVGr1 <- BASiCS_DetectHVG(chain.r1, VarThreshold = 0.6)
HVGr4 <- BASiCS_DetectHVG(chain.r4, VarThreshold = 0.6)

BASiCS_PlotVG(HVGut, "VG")
BASiCS_PlotVG(HVGcis, "VG")
BASiCS_PlotVG(HVGr1, "VG")
BASiCS_PlotVG(HVGr4, "VG")

#Store HVG as a dataframe
all.ut <- as.data.frame(HVGut@Table) #All genes with HVG=TRUE or HVG=FALSE
high.ut <- as.data.frame(HVGut) #Only HVG=TRUE genes 
all.cis <- as.data.frame(HVGcis@Table) 
high.cis <- as.data.frame(HVGcis)
all.r1 <- as.data.frame(HVGr1@Table) 
high.r1 <- as.data.frame(HVGr1)
all.r4 <- as.data.frame(HVGr4@Table) 
high.r4 <- as.data.frame(HVGr4)


