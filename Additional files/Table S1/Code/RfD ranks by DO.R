library(dplyr)
library(qdapTools)

RfD <- readRDS(file = "RfD for DO.RDS") #This data file is in the folder of Figure S7 results
RfD.CA <- RfD %>% filter(Compartment == "CA")

RfD.CA.rank <- RfD.CA %>% group_by(Duration) %>% mutate(ranking = rank(RfD,ties.method = "first")) #Give a ranking number to each disease


DO.sensitive.summary.CA <- RfD.CA.rank %>% group_by(DO) %>%
                                               summarise(averege.RfD = mean(RfD),
                                                         times = sum(ranking <= 20))

DO.sensitive.summary.CA <- DO.sensitive.summary.CA[order(-DO.sensitive.summary.CA$times,DO.sensitive.summary.CA$averege.RfD),]
#DO.CA  <- unique(RfD.CA.rank$DO)

#DO.sensitive.CA <- list()

#for (i in 1:length(DO.CA)){

#  disease <- DO.CA[i]
#  DO.sensitive.CA[[i]] <- RfD.CA.rank %>% filter(DO == disease & ranking <= 20)

#}

#DO.sensitive.CA <- list_df2df(DO.sensitive.CA)

#DO.sensitive.summary.CA <- DO.sensitive.CA %>% group_by(DO) %>%
#                                               summarise(times = sum(ranking <= 20),
#                                                         average.rank = mean(ranking),
#                                                         average.RfD = mean(RfD))
write.csv(DO.sensitive.summary.CA, "DO rank CA.csv",row.names = FALSE)
#saveRDS(DO.sensitive.summary.CA, file = "DO rank CA.RDS")

RfD.CL <- RfD %>% filter(Compartment ==  "CL")
RfD.CL.rank <- RfD.CL %>% group_by(Duration) %>% mutate(ranking = rank(RfD, ties.method = "first"))

DO.sensitive.summary.CL <- RfD.CL.rank %>% group_by(DO) %>%
                                           summarise(averege.RfD = mean(RfD),
                                           times = sum(ranking <= 20))

DO.sensitive.summary.CL <- DO.sensitive.summary.CL[order(-DO.sensitive.summary.CL$times,DO.sensitive.summary.CL$averege.RfD),]

#DO.CL <- unique(RfD.CL.rank$DO)

#DO.sensitive.CL <- list()

#for (i in 1: length(DO.CL)) {

#  disease <- DO.CL[i]
#  DO.sensitive.CL[[i]] <- RfD.CL.rank %>% filter(DO == disease & ranking <= 20)
#}

#DO.sensitive.CL <- list_df2df(DO.sensitive.CL)
#DO.sensitive.summary.CL <- DO.sensitive.CL %>% group_by(DO) %>% 
#                           summarise(times = sum(ranking <=20),
#                                     average.rank = mean(ranking),
#                                     average.RfD = mean(RfD))

#saveRDS(DO.sensitive.summary.CL, file  = "DO rank CL.RDS")
write.csv(DO.sensitive.summary.CL, "DO rank CL.csv",row.names = FALSE)

