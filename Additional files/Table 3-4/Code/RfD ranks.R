library(dplyr)
library(qdapTools)
# The data files are in the folder of Figure 5/Results
RfD.h <- readRDS(file = "Human RfD for significant pathways.RDS")

RfD.D1.h.CA <- RfD.h %>% filter(Duration == "Day 1" & Compartment == "CA")
RfD.D4.h.CA <- RfD.h %>% filter(Duration == "Day 4" & Compartment == "CA")
RfD.D10.h.CA <- RfD.h %>% filter(Duration == "Day 10" & Compartment == "CA")
RfD.D14.h.CA <- RfD.h %>% filter(Duration == "Day 14" & Compartment == "CA")

RfD.D1.h.CL <- RfD.h %>% filter(Duration == "Day 1" & Compartment == "CL")
RfD.D4.h.CL <- RfD.h %>% filter(Duration == "Day 4" & Compartment == "CL")
RfD.D10.h.CL <- RfD.h %>% filter(Duration == "Day 10" & Compartment == "CL")
RfD.D14.h.CL <- RfD.h %>% filter(Duration == "Day 14" & Compartment == "CL")

pathlist.h.D1.CA <- RfD.D1.h.CA[order(RfD.D1.h.CA$RfD),]%>% select(Pathway, RfD)
pathlist.h.D4.CA <- RfD.D4.h.CA[order(RfD.D4.h.CA$RfD),]%>% select(Pathway, RfD)
pathlist.h.D10.CA <- RfD.D10.h.CA[order(RfD.D10.h.CA$RfD),]%>% select(Pathway, RfD)
pathlist.h.D14.CA <- RfD.D14.h.CA[order(RfD.D14.h.CA$RfD),]%>% select(Pathway, RfD)

pathlist.h.D1.CL <- RfD.D1.h.CL[order(RfD.D1.h.CL$RfD),]%>% select(Pathway, RfD)
pathlist.h.D4.CL <- RfD.D4.h.CL[order(RfD.D4.h.CL$RfD),]%>% select(Pathway, RfD)
pathlist.h.D10.CL <- RfD.D10.h.CL[order(RfD.D10.h.CL$RfD),]%>% select(Pathway, RfD)
pathlist.h.D14.CL <- RfD.D14.h.CL[order(RfD.D14.h.CL$RfD),]%>% select(Pathway, RfD)

pathlist.h.D1.CA <- cbind(pathlist.h.D1.CA, rank = seq_along(pathlist.h.D1.CA$Pathway), Duration = "Human: 1 day")
pathlist.h.D4.CA <- cbind(pathlist.h.D4.CA, rank = seq_along(pathlist.h.D4.CA$Pathway), Duration = "Human: 4 day")
pathlist.h.D10.CA <- cbind(pathlist.h.D10.CA, rank = seq_along(pathlist.h.D10.CA$Pathway), Duration = "Human: 10 day")
pathlist.h.D14.CA <- cbind(pathlist.h.D14.CA, rank = seq_along(pathlist.h.D14.CA$Pathway),Duration = "Human: 14 day")

df.CA <- rbind.data.frame(pathlist.h.D1.CA,pathlist.h.D4.CA,pathlist.h.D10.CA,pathlist.h.D14.CA)

path.CA <- unique(df.CA$Pathway)

ranks.CA <- list()

for(i in 1:length(path.CA)){

  pathway = path.CA[i]
  ranks.CA[[i]] <- df.CA %>% filter(Pathway == pathway)
  print(ranks.CA[[i]])
  
}

ranks.CA <-list_df2df(ranks.CA)

ranks.average.CA <- ranks.CA %>% group_by(Pathway) %>% summarise(average.RfD = mean(RfD),
                                                      times = sum(rank < 20))


ranks.average.CA <- ranks.average.CA[order(-ranks.average.CA$times, ranks.average.CA$average.RfD),]
#saveRDS(ranks.average.CA, file = "ranks CA.RDS")
write.csv(ranks.average.CA, "Ranks CA_human.csv",row.names = FALSE)

pathlist.h.D1.CL <- cbind(pathlist.h.D1.CL, rank = seq_along(pathlist.h.D1.CL$Pathway), Duration = "Human: 1 day")
pathlist.h.D4.CL <- cbind(pathlist.h.D4.CL, rank = seq_along(pathlist.h.D4.CL$Pathway), Duration = "Human: 4 day")
pathlist.h.D10.CL <- cbind(pathlist.h.D10.CL, rank = seq_along(pathlist.h.D10.CL$Pathway), Duration = "Human: 10 day")
pathlist.h.D14.CL <- cbind(pathlist.h.D14.CL, rank = seq_along(pathlist.h.D14.CL$Pathway),Duration = "Human: 14 day")

df.CL <- rbind.data.frame(pathlist.h.D1.CL,pathlist.h.D4.CL,pathlist.h.D10.CL,pathlist.h.D14.CL)

path.CL <- unique(df.CL$Pathway)

ranks.CL <- list()

for(i in 1:length(path.CL)){
  
  pathway = path.CL[i]
  ranks.CL[[i]] <- df.CL %>% filter(Pathway == pathway)
  
}

ranks.CL <-list_df2df(ranks.CL)

ranks.average.CL <- ranks.CL %>% group_by(Pathway) %>% summarise(average.RfD = mean(RfD),
                                                                 times = sum(rank < 20))
                                                                 
ranks.average.CL <- ranks.average.CL[order(-ranks.average.CL$times, ranks.average.CL$average.RfD),]
#saveRDS(ranks.average.CL, file = "ranks CL.RDS")
write.csv(ranks.average.CL,"Ranks CL_human.csv",row.names = FALSE)

summary <- RfD.h %>% group_by(Compartment)%>%
                     summarise(lower = quantile(RfD,0.025), 
                               median = median(RfD),
                               upper = quantile(RfD,0.975))

