library(clusterProfiler)
library(AnnotationDbi)
library("mouse4302.db") #Affymetrix Mouse Genome 430 2.0 Array annotation data
library(cowplot)
library(dplyr)
library(plyr)
library(qdapTools)

# Import the BBMD results, the location is in the folder of Figure 2/Results
BBMD.mice <- read.csv("BBMD result mice.csv",header = TRUE)
BBMD.D1  <- read.csv("BBMD result D1.csv",header = TRUE)
BBMD.D4  <- read.csv("BBMD result D4.csv", header = TRUE)
BBMD.D10 <- read.csv("BBMD result D10.csv", header = TRUE)
BBMD.D14 <- read.csv("BBMD result D14.csv", header = TRUE)

# Compute the ratio of BMDU and BMDL
BBMD.D1.2  <- BBMD.D1 %>% mutate(BMDU.BMDL = bmdu/bmdl)
BBMD.D4.2  <- BBMD.D4 %>% mutate(BMDU.BMDL = bmdu/bmdl)
BBMD.D10.2 <- BBMD.D10 %>% mutate(BMDU.BMDL = bmdu/bmdl)
BBMD.D14.2 <- BBMD.D14 %>% mutate(BMDU.BMDL = bmdu/bmdl)

# Filter the significant genes with the ratio of BMDU/BMDL <40 and the BMD less than the highest dose
BBMD.mice.sig <- BBMD.mice %>% filter(BMDU.BMDL < 40 & bmd < 10)
BBMD.D1.sig   <- BBMD.D1.2 %>% filter(BMDU.BMDL < 40 & bmd < 20)
BBMD.D4.sig   <- BBMD.D4.2 %>% filter(BMDU.BMDL < 40 & bmd < 20)
BBMD.D10.sig  <- BBMD.D10.2 %>% filter(BMDU.BMDL < 40 & bmd < 20)
BBMD.D14.sig  <- BBMD.D14.2 %>% filter(BMDU.BMDL < 40 & bmd < 20)

# Import the KEGG results, the location is in the folder of Figure 3/Results
pathway.mice <- readRDS(file = "KEGG mice sig.RDS")
pathway.D1 <- readRDS(file = "KEGG human D1.RDS")
pathway.D4 <- readRDS(file = "KEGG human D4.RDS")
pathway.D10 <- readRDS(file = "KEGG human D10.RDS")
pathway.D14 <- readRDS(file = "KEGG human D14.RDS")

#Filter the significant KEGG pathways
path.sig.m <- pathway.mice %>% filter(p.adjust < 0.05 & qvalue <0.05)
path.sig.D1 <- pathway.D1 %>% filter(p.adjust < 0.05 & qvalue <0.05)
path.sig.D4 <- pathway.D4 %>% filter(p.adjust < 0.05 & qvalue <0.05)
path.sig.D10 <- pathway.D10 %>% filter(p.adjust < 0.05 & qvalue <0.05)
path.sig.D14 <- pathway.D14 %>% filter(p.adjust < 0.05 & qvalue <0.05)

######################################### pathway HED & RfD calculation ###############################
# HED and RfD calculation for significant pathways
path.sig.list <-list(path.sig.D1, path.sig.D4, path.sig.D10, path.sig.D14) 
BBMD.list <-list(BBMD.D1, BBMD.D4, BBMD.D10, BBMD.D14)
Css <- readRDS(file = "Css simulation.RDS") # Css results are in the Table 2/Results

TD.pathway.df <- list()

for(i in 1:4){
  
  path.sig <- path.sig.list[[i]]
  BBMD <- BBMD.list[[i]]
  probe <- ldply(strsplit(BBMD$ï..Probe.ID,"_"))[,1]
  BBMD <- cbind(BBMD, probe)
  
  HED.path <- list()
  TD.pathway <- list()
  
  for (j in 1:nrow(path.sig)){
    
    genes <- path.sig[j,]$geneID 
    gene.ID <- unlist(strsplit(genes,"/"))
    gene.symbol <- bitr(gene.ID,fromType = "ENTREZID", toType = "SYMBOL",OrgDb = "org.Hs.eg.db")[,2]
    
    BBMD.cb <- BBMD %>% filter(probe %in% gene.symbol) # Filter the gene in the pathway
    genelist <- BBMD.cb[!duplicated(BBMD.cb$Gene.Symbol),] #Drop the same gene with larger adjusted p-value
    
    HED <-list() 
    
    for(k in 1: nrow(genelist)){
      
      HED[[k]] <- cbind(HED.CA = genelist$bmd[k]/Css$Css.CA,  # The HED calculated with Css in plasma
                       HED.CL = genelist$bmd[k]/Css$Css.CL)  # The HED calculated with Css in the liver
      names(HED)[k] <- genelist$Gene.Symbol[k]
      
    }
    HED.simulation <- list_df2df(HED)
    TD.path.CA <- c(Compartment = "CA", TD = quantile(HED.simulation$HED.CA, prob = 0.05)) # The threshold dose: 5% of the HED simulation
    TD.path.CL <- c(Compartment = "CL", TD = quantile(HED.simulation$HED.CL, prob = 0.05))
    
    TD.pathway[[j]] <- rbind(TD.path.CA,TD.path.CL)
    names(TD.pathway)[j] <- path.sig$Description[j]
  }
  TD.pathway.df[[i]] <- list_df2df(TD.pathway)
}

names(TD.pathway.df) <- c("Day 1", "Day 4", "Day 10", "Day 14")
TD.sigpath.df <-list_df2df(TD.pathway.df)
colnames(TD.sigpath.df) <- c("Duration","Pathway","Compartment","TD")

TD.sigpath.df$TD <- as.numeric(TD.sigpath.df$TD)
RfD.sigpath.df <- TD.sigpath.df %>% mutate(RfD = TD / 10) #Reference dose = TD/UF, a uncertainty factor of 10 is used to account for the intraspecies extrapolation

saveRDS(RfD.sigpath.df,file = "reference dose for significant pathways of Human.RDS")


## Calculate the pathway HED and RfD for the pathways matched to mice significant pathways
pathway.match <- list(path.match.D1 = pathway.D1 %>% filter(Description %in% path.sig.m$Description),
                      path.match.D4 = pathway.D4 %>% filter(Description %in% path.sig.m$Description),
                      path.match.D10 = pathway.D10 %>% filter(Description %in% path.sig.m$Description),
                      path.match.D14 = pathway.D14 %>% filter(Description %in% path.sig.m$Description))

TD.match.df <- list()
HED.matchpath.df <- list()

for(i in 1:4){

  path.match <- pathway.match[[i]]
  BBMD <- BBMD.list[[i]]
  probe <- ldply(strsplit(BBMD$ï..Probe.ID,"_"))[,1]
  BBMD <- cbind(BBMD, probe)
  
  HED.match <- list()
  TD.match  <- list()
  
  for (j in 1:nrow(path.match)){
 
    genes <- path.match[j,]$geneID
    gene.ID <- unlist(strsplit(genes,"/"))
    gene.symbol <- bitr(gene.ID,fromType = "ENTREZID", toType = "SYMBOL",OrgDb = "org.Hs.eg.db")[,2]
    
    BBMD.cb <- BBMD %>% filter(probe %in% gene.symbol) # Filter the gene in the pathway
    genelist <- BBMD.cb[!duplicated(BBMD.cb$probe),] #Drop the same gene with larger adjusted p-value
    
    HED <- list() 
    for(k in 1: nrow(genelist)){
      
      HED[[k]] <- cbind(HED.CA = genelist$bmd[k]/Css$Css.CA,
                       HED.CL = genelist$bmd[k]/Css$Css.CL)
      names(HED)[k] <- genelist$Gene.Symbol[k]
      
    }
    HED.simulation <- list_df2df(HED)
    TD.match.CA <- c(Compartment = "CA", TD = quantile(HED.simulation$HED.CA, prob = 0.05))
    TD.match.CL <- c(Compartment = "CL", TD = quantile(HED.simulation$HED.CL, prob = 0.05))
    
    TD.match[[j]] <- rbind(TD.match.CA,TD.match.CL)
    names(TD.match)[j] <- path.match$Description[j]
    
    HED.CA <- c(Compartment = "CA", HED = median(HED.simulation$HED.CA))
    HED.CL <- c(Compartment = "CL", HED = median(HED.simulation$HED.CL))
    
    HED.match[[j]] <- rbind(HED.CA,HED.CL)
  }
  names(HED.match) <- path.match$Description
  HED.matchpath.df[[i]] <- list_df2df((HED.match))
  
  TD.match.df[[i]] <- list_df2df(TD.match)
}
names(HED.matchpath.df) <- c("Day 1", "Day 4", "Day 10", "Day 14")
HED.matchpath.df <- list_df2df(HED.matchpath.df)
colnames(HED.matchpath.df) <- c("Duration","Pathway","Compartment","HED")

TD.matchpath.df <- list_df2df(TD.match.df)
colnames(TD.matchpath.df) <- c("Duration","Pathway","Compartment","TD")

TD.matchpath.df$TD <- as.numeric(TD.matchpath.df$TD)
RfD.matchpath.df <- TD.matchpath.df %>% mutate(RfD = TD / 10)

saveRDS(RfD.matchpath.df,file = "Human RfD for matched pathways.RDS")
saveRDS(HED.matchpath.df,file = "Human HED for matched pathways.RDS")

## Duration HED & RfD calculation
HED.duration <- list()
TD.duration  <- list()

for(i in 1:4) {

  path.sig <- path.sig.list[[i]]
  genes <- paste(path.sig$geneID,collapse = "/")
  gene.ID <- unlist(strsplit(genes,"/"))
  gene.ID.uni <- unique(gene.ID)
  gene.symbol <- bitr(gene.ID.uni, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")[,2]
  
  BBMD <- BBMD.list[[i]]
  probe <- ldply(strsplit(BBMD$ï..Probe.ID,"_"))[,1]
  BBMD <- cbind(BBMD, probe)
  
  BBMD.duration <- BBMD %>% filter(probe %in% gene.symbol)
  BBMD.duration <- BBMD.duration[order(BBMD.duration$probe,BBMD.duration$Adjusted.P.Value),]
  BBMD.duration2 <- BBMD.duration[!duplicated(BBMD.duration$probe),]
  
  HED <- list()
  for(j in 1:nrow(BBMD.duration2)){
    
    HED[[j]] <- cbind(HED.CA = BBMD.duration2$bmd[j]/Css$Css.CA,
                      HED.CL = BBMD.duration2$bmd[j]/Css$Css.CL) 
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
  HED.duration[[i]] <- rbind(HED.duration.CA,HED.duration.CL) 
  
  TD.duration.CA <- c(Compartment = "CA", TD = quantile(HED.simulation$HED.CA, prob = 0.05))
  TD.duration.CL <- c(Compartment = "CL", TD = quantile(HED.simulation$HED.CL, prob = 0.05))
  TD.duration[[i]] <-rbind(TD.duration.CA,TD.duration.CL)
  
}

names(HED.duration) <- c("Day 1", "Day 4", "Day 10", "Day 14")
names(TD.duration) <- c("Day 1", "Day 4", "Day 10", "Day 14")

HED.duration.df <- list_df2df(HED.duration) # Duration HED
colnames(HED.duration.df) <- c("Duration","min","low","mid","top","max","Compartment") 

TD.duration.df <- list_df2df(TD.duration) #Duration TD
colnames(TD.duration.df) <- c("Duration","Compartment","TD")

TD.duration.df$TD <- as.numeric(TD.duration.df$TD)
RfD.duration.df <- TD.duration.df %>% mutate(RfD = TD / 10)  #RfD for each duration

saveRDS(HED.duration.df, file = "duration HED for human.RDS")
saveRDS(RfD.duration.df, file = "duration RfD for human.RDS")

