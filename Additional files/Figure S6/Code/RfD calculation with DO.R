library(clusterProfiler)
library("org.Hs.eg.db")
library(qdapTools)
library(dplyr)
library(plyr)
library(DOSE)

#Loading the BBMD results in the folder of Figure 2/Results
BBMD.list <- list(BBMD.D1 = read.csv("BBMD result D1.csv", header = TRUE),   #Import the data of BBMD for Day 1
                  BBMD.D4 = read.csv("BBMD result D4.csv", header = TRUE),   #Import the data of BBMD for Day 4
                  BBMD.D10 = read.csv("BBMD result D10.csv", header = TRUE), #Import the data of BBMD for Day 10
                  BBMD.D14 = read.csv("BBMD result D14.csv", header = TRUE)) #Import the data of BBMD for Day 14

#Loading the Css simulation results
Css <- readRDS(file = "Css simulation.RDS")

#Filter the significant BBMD results
BBMD.list.sig <-list()

for(i in 1:4) {

  BBMD <- BBMD.list[[i]]
  BBMD.new <- BBMD %>% mutate(BMDU.BMDL = bmdu/bmdl) #Compute the ratio of BMDU and BMDL
  BBMD.list.sig[[i]] <- BBMD.new %>% filter(BMDU.BMDL < 40 & bmd < 20) #Remove the BMDU/BMDL ratio larger than 40 and the BMD larger than the highest dose (20uM)
  
}
#Loading DO analysis results in the folder of Figure S1-Figure S4
DO.D1 <- readRDS(file = "DO human D1.RDS")
DO.D4 <- readRDS(file = "DO human D4.RDS")
DO.D10 <- readRDS(file = "DO human D10.RDS")
DO.D14 <- readRDS(file = "DO human D14.RDS")

#Get the significant DO results
DO.sig.D1 <- DO.D1 %>% filter(p.adjust < 0.05 & qvalue < 0.05)
DO.sig.D4 <- DO.D4 %>% filter(p.adjust < 0.05 & qvalue < 0.05)
DO.sig.D10 <- DO.D10 %>% filter(p.adjust < 0.05 & qvalue < 0.05)
DO.sig.D14 <- DO.D14 %>% filter(p.adjust < 0.05 & qvalue < 0.05)

DO.results.sig <- list(DO.sig.D1, DO.sig.D4, DO.sig.D10, DO.sig.D14)
names(DO.results.sig) <- c("Human: Day 1", "Human: Day 4", "Human: Day 10", "Human: Day 14")

HED.DO.all <- list()
TD.DO.all  <- list()

for(i in 1:4){

  DO.result <- DO.results.sig[[i]]
  BBMD <- BBMD.list.sig [[i]] %>% mutate (probe = ldply(strsplit(ï..Probe.ID, "_"))[,1])
  
  HED.DO <-list()
  TD.sig <- list()
  
  for (j in 1: nrow(DO.result)){
 
  genes <- DO.result$geneID[j]
  gene.ID <- unlist(strsplit(genes, "/"))
  gene.symbol <- bitr(gene.ID, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")[,2]
  BBMD.DO <- BBMD %>% filter(probe %in% gene.symbol) #Match the geneID in the DO results to the BBMD results
  
  BBMD.DO <-BBMD.DO[order(BBMD.DO$probe,BBMD.DO$Adjusted.P.Value),]
  BBMD.DO2 <-BBMD.DO[!duplicated(BBMD.DO$probe),]
  
  HED <- list()
 
  for(k in 1:nrow(BBMD.DO2)){

    HED[[k]] <- cbind(HED.CA = BBMD.DO2$bmd[k]/Css$Css.CA,
                      HED.CL = BBMD.DO2$bmd[k]/Css$Css.CL) 
    
  }
  
  HED.simulation <- list_df2df(HED)
  
  # Summarize the plasma HED for significant pathways
  HED.sig.CA <- c(min = quantile(HED.simulation$HED.CA,prob = 0.025),
                  low = quantile(HED.simulation$HED.CA,prob = 0.25),
                  mid = quantile(HED.simulation$HED.CA,prob = 0.5),
                  top = quantile(HED.simulation$HED.CA,prob = 0.75),
                  max = quantile(HED.simulation$HED.CA,prob = 0.975),
                  Compartment = "CA")
  
  # Summarize the liver HED for significant pathways
  HED.sig.CL <- c(min = quantile(HED.simulation$HED.CL,prob = 0.025),
                  low = quantile(HED.simulation$HED.CL,prob = 0.25),
                  mid = quantile(HED.simulation$HED.CL,prob = 0.5),
                  top = quantile(HED.simulation$HED.CL,prob = 0.75),
                  max = quantile(HED.simulation$HED.CL,prob = 0.975),
                  Compartment = "CL")
  
  # Create the data frame  for HED summary for a Disease
  HED.DO[[j]] <-rbind(HED.sig.CA,HED.sig.CL) 
  
  TD.sig.CA <- c(Compartment = "CA", TD = quantile(HED.simulation$HED.CA, prob = 0.05))
  TD.sig.CL <- c(Compartment = "CL", TD = quantile(HED.simulation$HED.CL, prob = 0.05))
  
  TD.sig[[j]] <-rbind(TD.sig.CA, TD.sig.CL)
  }
  names(HED.DO) <- DO.result$Description
  HED.DO.all[[i]] <- list_df2df(HED.DO)
  
  names(TD.sig) <-DO.result$Description
  TD.DO.all[[i]] <- list_df2df(TD.sig)
}
names(HED.DO.all) <- c("Human: Day 1", "Human: Day 4", "Human: Day 10", "Human: Day 14")
names(TD.DO.all) <- c("Day 1", "Day 4", "Day 10", "Day 14")

HED.DO.df <- list_df2df(HED.DO.all)
colnames(HED.DO.df) <- c("Duration","DO","min","low","mid","top","max","Compartment")

TD.DO.df <- list_df2df(TD.DO.all)
colnames(TD.DO.df) <- c("Duration","DO","Compartment","TD")
TD.DO.df$TD <- as.numeric(TD.DO.df$TD)


RfD.DO.df <- TD.DO.df %>% mutate(RfD = TD / 10)  #RfD for each significant DO


HED.DO.duration <-list()
TD.duration <- list()

for(i in 1:4) {

  DO.result <- DO.results.sig[[i]]
  genes <- paste(DO.result$geneID, collapse = "/")
  gene.ID <- unlist(strsplit(genes, "/"))
  gene.ID.uni <- unique(gene.ID)
  gene.symbol <- bitr(gene.ID.uni, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")[,2]
  
  BBMD <- BBMD.list.sig [[i]] %>% mutate (probe = ldply(strsplit(ï..Probe.ID,"_"))[,1])
  BBMD.DO <- BBMD %>% filter(probe %in% gene.symbol)
  BBMD.DO <- BBMD.DO[order(BBMD.DO$probe,BBMD.DO$Adjusted.P.Value),]
  BBMD.DO2 <- BBMD.DO[!duplicated(BBMD.DO$probe),]
  
  HED <- list()
  for(j in 1:nrow(BBMD.DO2)){
    
    HED[[j]] <- cbind(HED.CA = BBMD.DO2$bmd[j]/Css$Css.CA,
                      HED.CL = BBMD.DO2$bmd[j]/Css$Css.CL) 
  }
  
  HED.simulation <-list_df2df(HED)
  
  # Summarize the plasma HED for significant pathways
  HED.duration.CA <- c(min = quantile(HED.simulation$HED.CA,prob = 0.025),
                       low = quantile(HED.simulation$HED.CA,prob = 0.25),
                       mid = quantile(HED.simulation$HED.CA,prob = 0.5),
                       top = quantile(HED.simulation$HED.CA,prob = 0.75),
                       max = quantile(HED.simulation$HED.CA,prob = 0.975),
                       Compartment = "CA")
  
  # Summarize the liver HED for significant pathways
  HED.duration.CL <- c(min = quantile(HED.simulation$HED.CL,prob = 0.025),
                       low = quantile(HED.simulation$HED.CL,prob = 0.25),
                       mid = quantile(HED.simulation$HED.CL,prob = 0.5),
                       top = quantile(HED.simulation$HED.CL,prob = 0.75),
                       max = quantile(HED.simulation$HED.CL,prob = 0.975),
                       Compartment = "CL")
  
  # Create the data frame  for HED summary for a Disease
  HED.DO.duration[[i]] <-rbind(HED.duration.CA, HED.duration.CL) 
  
  TD.duration.CA <- c(Compartment = "CA", TD = quantile(HED.simulation$HED.CA, prob = 0.05))
  TD.duration.CL <- c(Compartment = "CL", TD = quantile(HED.simulation$HED.CL, prob = 0.05))
  TD.duration[[i]] <- rbind(TD.duration.CA,TD.duration.CL)
  
}

names(HED.DO.duration) <- c("Human: Day 1", "Human: Day 4", "Human: Day 10", "Human: Day 14")
names(TD.duration) <- c("Human: Day 1", "Human: Day 4", "Human: Day 10", "Human: Day 14")

HED.duration.df <- list_df2df(HED.DO.duration)
colnames(HED.duration.df) <- c("Duration","min","low","mid","top","max","Compartment")

TD.duration.df <- list_df2df(TD.duration)
colnames(TD.duration.df) <- c("Duration","Compartment","TD")

TD.duration.df$TD <- as.numeric(TD.duration.df$TD)
RfD.duration.df <- TD.duration.df %>% mutate(RfD = TD / 10)  #RfD for each significant DO

saveRDS(HED.DO.df, file = "HED for DO.RDS")
saveRDS(RfD.DO.df, file = "RfD for DO.RDS")
saveRDS(HED.duration.df, file = "HED duration.RDS")
saveRDS(RfD.duration.df, file = "RfD duration.RDS")
