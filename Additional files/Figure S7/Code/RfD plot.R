library(ggplot2)
library(dplyr)

#Open the Result files in the folder of Figure S7/Results
df.RfD.duration <- readRDS("RfD duration.RDS")
df.RfD.DO <- readRDS("RfD for DO.RDS")

df.RfD.DO$RfD <- (df.RfD.DO$RfD)*1e6
df.RfD.duration$RfD <- (df.RfD.duration$RfD)*1e6

df.RfD.DO.log <- df.RfD.DO %>% mutate(RfD.log = log(RfD,10))
df.RfD.duration.log <-df.RfD.duration %>% mutate(RfD.log = log(RfD,10))
df.RfD.duration.log[df.RfD.duration.log == "Human: Day 1"] <- "Day 1"
df.RfD.duration.log[df.RfD.duration.log == "Human: Day 4"] <- "Day 4"
df.RfD.duration.log[df.RfD.duration.log == "Human: Day 10"] <- "Day 10"
df.RfD.duration.log[df.RfD.duration.log == "Human: Day 14"] <- "Day 14"

df.RfD.DO.log[df.RfD.DO.log == "CA"] <-"Plasma"
df.RfD.DO.log[df.RfD.DO.log == "CL"] <-"Liver"


RfD.EPA <- 20 #ng/kg/day
RfD.EFSA.2018 <- 1.8 #ng/kg/day
RfD.EFSA.2020 <- 0.63 #ng/kg/day
RfD.EPA.log <-log(RfD.EPA, 10)
RfD.EFSA.2018.log <- log(RfD.EFSA.2018,10)
RfD.EFSA.2020.log <- log(RfD.EFSA.2020,10)

plot.RfD.DO <- ggplot()+
  geom_point(data = df.RfD.DO.log,aes(x = Duration, y= RfD.log, color = Compartment),
             alpha = 0.8,position = position_jitterdodge(jitter.width = 0.1), size = 1.7,show.legend = TRUE)+
  scale_x_discrete(limits = c("Day 1", "Day 4", "Day 10", "Day 14"))+
  scale_fill_manual(values = c("yellow","green"),name = "Duration RfD",labels = c("Plasma","Liver"))+
  geom_point(data = df.RfD.duration.log,aes(x = Duration, y = RfD.log,fill = Compartment), 
             pch = 23, size =3, color = "Black", show.legend = NA,position = position_jitterdodge(jitter.width = 0.1))+
  ylab(expression(paste(log[10],"RfD")))+ #y-axis labeling
  geom_hline(aes(yintercept = RfD.EPA.log,linetype = "EPA RfD"),color = "red")+ #Dashline for EPA RfD
  geom_hline(aes(yintercept = RfD.EFSA.2018.log,linetype = "EFSA TDI, 2018"),color = "blue")+ #Dashline for 2018 EFSA RfD
  geom_hline(aes(yintercept = RfD.EFSA.2020.log,linetype = "EFSA TDI, 2020"),color = "purple")+ #Dashline for 2020 EFSA RfD
  scale_linetype_manual(name = "Limit", values = c(2,2,2),guide = guide_legend(override.aes = list(color = c("blue","purple","red"))))+
  theme_bw()+ #white background
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line())+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "serif",size = 12),
        axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 12),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"))+
  labs(fill = "Duration RfD")+
  guides(color = guide_legend(order = 1),
         fill = guide_legend(order = 2),
         line = guide_legend(order = 3))
  

ggsave("Figure S7.tiff", scale = 1,
       plot = plot.RfD.DO,
       width = 12, height = 6, dpi = 320)
