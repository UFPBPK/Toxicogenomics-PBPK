library(dplyr)
library(qdapTools)

RfD <- readRDS(file = "RfD for DO.RDS") #This data file is in the folder of Figure S10 results
RfD.CA <- RfD %>% filter(Compartment == "CA")

RfD.CA.rank <- RfD.CA %>% group_by(Duration) %>% mutate(ranking = rank(RfD,ties.method = "first")) #Give a ranking number to each disease


DO.sensitive.summary.CA <- RfD.CA.rank %>% group_by(DO) %>%
                                               summarise(averege.RfD = mean(RfD),
                                                         times = sum(ranking <= 20))

DO.sensitive.summary.CA <- DO.sensitive.summary.CA[order(-DO.sensitive.summary.CA$times,DO.sensitive.summary.CA$averege.RfD),]

write.csv(DO.sensitive.summary.CA, "DO rank CA.csv",row.names = FALSE)


RfD.CL <- RfD %>% filter(Compartment ==  "CL")
RfD.CL.rank <- RfD.CL %>% group_by(Duration) %>% mutate(ranking = rank(RfD, ties.method = "first"))

DO.sensitive.summary.CL <- RfD.CL.rank %>% group_by(DO) %>%
                                           summarise(averege.RfD = mean(RfD),
                                           times = sum(ranking <= 20))

DO.sensitive.summary.CL <- DO.sensitive.summary.CL[order(-DO.sensitive.summary.CL$times,DO.sensitive.summary.CL$averege.RfD),]


write.csv(DO.sensitive.summary.CL, "DO rank CL.csv",row.names = FALSE)

