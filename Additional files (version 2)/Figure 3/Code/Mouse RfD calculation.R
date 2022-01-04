library(dplyr)
library("mouse4302.db") #Mouse Genome 430 2.0 Array annotation data
library(qdapTools)

# HED and RfD calculation of mice
#AUC results are located in the folder of Table 1/Results
AUC <- readRDS(file = "AUC simulation.RDS")

AUC.idv  = as.data.frame(AUC) # make a copy
## Estimated the average serum concentration (ASC)/average liver concentration (ALC)
AUC.idv <- AUC.idv %>% transmute(Avg.CA.mouse   = AUC.CA.mouse/(7*24), # mouse ASC: AUC normalized by exposure duration (7 days or 7*24 hours) 
                                 Avg.CA.human.m = AUC.CA.human.m/(365*25*24), # Human ASC: AUC normalized by the duration reaching the steady state (Based on the simulation results, it reached steady state at ~25 years.)   
                                 Avg.CL.mouse   = AUC.CL.mouse/(7*24),
                                 Avg.CL.human.m = AUC.CL.human.m/(365*25*24))

AUC.ratio.CA <- (AUC.idv$Avg.CA.mouse)/(AUC.idv$Avg.CA.human.m)
AUC.ratio.CL <- (AUC.idv$Avg.CL.mouse)/(AUC.idv$Avg.CL.human.m)

# Calculation by KEGG pathways

BBMD <- read.csv(file ="BBMD result mice.csv") # "BBMD result mice.csv" is in the folder of Figure S1/Results
BBMD.sig <- BBMD %>% filter(BMDU.BMDL < 40 & bmd < 10)

KEGG.m <- readRDS(file = "KEGG mice sig.RDS") # "KEGG mice sig.RDS" is in the folder of Figure S2/Results

OUT <- select(mouse4302.db,BBMD$ï..Probe.ID,c("SYMBOL","ENTREZID"))
OUT.uni <- OUT[!duplicated(OUT),]

HED.path <- list()

for(i in 1:7){

  genes <- KEGG.m$geneID[i]
  gene.ID <- unlist(strsplit(genes,"/"))
  
  gene.symbol <- OUT.uni %>% filter(ENTREZID %in% gene.ID) %>% dplyr::select(SYMBOL,ENTREZID)
  gene.symbol <- apply(gene.symbol,2,as.character)
  BBMD.path <- BBMD.sig %>% filter(Gene.Symbol %in% gene.symbol)
  BBMD.path <- BBMD.path[order(BBMD.path$Gene.Symbol,BBMD.path$Adjusted.P.Value),]
  BBMD.path.uni <- BBMD.path[!duplicated(BBMD.path$Gene.Symbol),]
  
  HED.iter <-list()
  
  for (j in 1:nrow(BBMD.path.uni)){
   
    POD <- BBMD.path.uni$bmd[j]
    HED.iter.CA = POD * AUC.ratio.CA #
    HED.iter.CL = POD * AUC.ratio.CL
    HED.iter[[j]] <- cbind(HED.iter.CA,HED.iter.CL)
    names(HED.iter)[j] <- BBMD.path.uni$Gene.Symbol[j]
    
  }
  
  HED.path[[i]] <- list_df2df(HED.iter)
}

names(HED.path) <-KEGG.m$Description

HED.simulation <-list_df2df(HED.path)
colnames(HED.simulation) <- c("Pathway","Gene","HED.CA","HED.CL")

#Summarize the plasma HED by pathway
HED.summary.CA <- HED.simulation %>% dplyr::group_by(Pathway) %>%
  dplyr::summarise(min = quantile(HED.CA, prob = 0.025),
                   low = quantile(HED.CA, prob = 0.25),
                   mid = quantile(HED.CA, prob = 0.5),
                   top = quantile(HED.CA, prob = 0.75),
                   max = quantile(HED.CA, prob = 0.975),
                   compartment = "CA")
#Summarize the liver HED by pathway
HED.summary.CL <- HED.simulation %>% dplyr::group_by(Pathway) %>%
  dplyr::summarise(min = quantile(HED.CL, prob = 0.025),
                   low = quantile(HED.CL, prob = 0.25),
                   mid = quantile(HED.CL, prob = 0.5),
                   top = quantile(HED.CL, prob = 0.75),
                   max = quantile(HED.CL, prob = 0.975),
                   compartment = "CL")

#Summarize the plasma HED for all mice significant pathways
HED.summary.all.CA <- HED.simulation %>% dplyr::summarise(Pathway = "All",
                                                          min = quantile(HED.CA, prob = 0.025),
                                                          low = quantile(HED.CA, prob = 0.25),
                                                          mid = quantile(HED.CA, prob = 0.5),
                                                          top = quantile(HED.CA, prob = 0.75),
                                                          max = quantile(HED.CA, prob = 0.975),
                                                          compartment = "CA")

#Summarize the liver HED for all mice significant pathways
HED.summary.all.CL <- HED.simulation %>% dplyr::summarise(Pathway = "All",
                                                          min = quantile(HED.CL, prob = 0.025),
                                                          low = quantile(HED.CL, prob = 0.25),
                                                          mid = quantile(HED.CL, prob = 0.5),
                                                          top = quantile(HED.CL, prob = 0.75),
                                                          max = quantile(HED.CL, prob = 0.975),
                                                          compartment = "CL")


#Combine the HEDs to a table
HED.summary.mice <- rbind(HED.summary.CA,HED.summary.CL,HED.summary.all.CA,HED.summary.all.CL)

HED.pathway.mice <- HED.summary.mice[1:14,]
HED.duration.mice <- HED.summary.mice[15:16,]

saveRDS(HED.pathway.mice, file = "Mouse HED for pahtways.RDS")
saveRDS(HED.duration.mice, file = "Mouse duration HED. RDS")

TD.CA <- HED.simulation %>% dplyr::group_by(Pathway) %>%
                         dplyr::summarise(TD = quantile(HED.CA, prob = 0.05),
                                          compartment = "CA")

TD.CL <- HED.simulation %>% dplyr::group_by(Pathway) %>%
                            dplyr::summarise(TD = quantile(HED.CL, prob = 0.05),
                                             compartment = "CL")

TD <-rbind(TD.CA, TD.CL)
RfD <- TD %>% mutate(RfD = TD/30) # RfD = TD/uncertainty factor of 30; unit: mg/kg/day

TD.duration.CA <- HED.simulation %>% dplyr::summarise(TD = quantile(HED.CA, prob = 0.05),
                                                      compartment = "CA")
TD.duration.CL <- HED.simulation %>% dplyr::summarise(TD = quantile(HED.CL, prob = 0.05),
                                                      compartment = "CL")

TD.duration <- rbind(TD.duration.CA, TD.duration.CL)
RfD.duration <- TD.duration %>% mutate (RfD = TD/30)

saveRDS(RfD, file = "Mouse RfD for significant pathways.RDS")
saveRDS(RfD.duration, file = "Mouse duration RfD.RDS")

