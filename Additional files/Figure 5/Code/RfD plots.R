library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)

#Loading human RfD (The location is Figure 5/Results)

RfD.sig.h <- readRDS(file = "Human RfD for significant pathways.RDS")
RfD.duration.h <- readRDS(file = "Human duration RfD.RDS")
RfD.matched.h <- readRDS(file = "Human RfD for matched pathways.RDS")

#Loading mouse RfD

RfD.sig.m <- readRDS(file = "Mouse RfD for significant pathways.RDS")
RfD.duration.m <- readRDS(file = "Mouse duration RfD.RDS")

RfD.sig.m <- cbind(Duration = "Mouse: Day 7",RfD.sig.m)
RfD.duration.m <- cbind(Duration = "Mouse: Day 7", RfD.duration.m)

RfD.sig.h[RfD.sig.h == "Day 1"]  <- "Human: Day 1"
RfD.sig.h[RfD.sig.h == "Day 4"]  <- "Human: Day 4"
RfD.sig.h[RfD.sig.h == "Day 10"] <- "Human: Day 10"
RfD.sig.h[RfD.sig.h == "Day 14"] <- "Human: Day 14"

RfD.duration.h[RfD.duration.h == "Day 1"]  <- "Human: Day 1"
RfD.duration.h[RfD.duration.h == "Day 4"]  <- "Human: Day 4"
RfD.duration.h[RfD.duration.h == "Day 10"] <- "Human: Day 10"
RfD.duration.h[RfD.duration.h == "Day 14"] <- "Human: Day 14"

RfD.matched.h[RfD.matched.h == "L1"] <- "Human: Day 1"
RfD.matched.h[RfD.matched.h == "L2"] <- "Human: Day 4"
RfD.matched.h[RfD.matched.h == "L3"] <- "Human: Day 10"
RfD.matched.h[RfD.matched.h == "L4"] <- "Human: Day 14"

RfD.sig.h$RfD <- RfD.sig.h$RfD*1e6 #convert the unit to ng/kg/day
RfD.duration.h$RfD <- RfD.duration.h$RfD*1e6
RfD.matched.h$RfD <- RfD.matched.h$RfD*1e6

RfD.sig.m$RfD <- RfD.sig.m$RfD*1e6
RfD.duration.m$RfD <- RfD.duration.m$RfD*1e6

RfD.EPA <- 20 #ng/kg/day
RfD.EFSA.2018 <- 1.8 #ng/kg/day
RfD.EFSA.2020 <- 0.63 #ng/kg/day
RfD.EPA.log <-log(RfD.EPA, 10)
RfD.EFSA.2018.log <- log(RfD.EFSA.2018,10)
RfD.EFSA.2020.log <- log(RfD.EFSA.2020,10)

#######################################CA plot###################################################
RfD.sig.CA.h <- RfD.sig.h %>% filter(Compartment == "CA")
RfD.duration.CA.h <- RfD.duration.h %>% filter(Compartment == "CA")
RfD.matched.CA.h <- RfD.matched.h %>% filter(Compartment == "CA")
RfD.sig.CA.m <- RfD.sig.m %>% filter(compartment == "CA")
RfD.duration.CA.m <- RfD.duration.m %>% filter(compartment == "CA")

RfD.sig.CA.h.log <- RfD.sig.CA.h %>% mutate(RfD.log = log(RfD,10))
RfD.duration.CA.h.log <- RfD.duration.CA.h %>% mutate(RfD.log = log(RfD,10))
RfD.matched.CA.h.log <- RfD.matched.CA.h %>% mutate(RfD.log = log(RfD,10))
RfD.sig.CA.m.log <- RfD.sig.CA.m %>% mutate(RfD.log = log(RfD,10))
RfD.duration.CA.m.log <- RfD.duration.CA.m %>% mutate(RfD.log = log(RfD,10))

colnames(RfD.sig.CA.m.log)[4] <- "Compartment"
RfD.duration.CA.m.log <- RfD.duration.CA.m.log %>% dplyr::select(Duration,compartment,TD, RfD, RfD.log)
colnames(RfD.duration.CA.m.log) <- colnames(RfD.duration.CA.h.log)


df.sig.CA <- rbind(RfD.sig.CA.h.log,RfD.sig.CA.m.log)
df.duration.CA <- rbind(RfD.duration.CA.h.log,RfD.duration.CA.m.log)
df.matched.CA <- rbind(RfD.matched.CA.h.log, RfD.sig.CA.m.log)


plot.CA <- ggplot()+
  geom_point(data = df.sig.CA, aes(x = Duration, y = RfD.log,fill = "Pathway RfD"), size = 2.5, show.legend = NA, color = "tomato",pch = 16 )+ #Pathway RfD
  scale_fill_manual(values = c("tomato","yellow"),name = "RfD",labels = c("Pathway RfD","Duration RfD"))+ #the color of the pathway RfD and duration RfD
  scale_x_discrete(limits = c("Human: Day 1", "Human: Day 4", "Human: Day 10", "Human: Day 14","Mouse: Day 7"))+ #the order of category of the x-axis
  geom_point(data = df.duration.CA,aes(x = Duration, y = RfD.log, fill= "yellow"), pch = 23, size =3, color = "Black", show.legend = NA)+ #duration RfD
  geom_point(data = df.matched.CA, aes(x = Duration, y = RfD.log, shape = Pathway),size = 1.5, color = "lightseagreen",position = position_dodge(0.3, preserve = "single"))+ #matched RfD
  scale_shape_manual("Matched pathways", values = 0:length(unique(df.matched.CA$Pathway)))+ #the shapes of matched RfD
  ylab(expression(paste(log[10],"RfD")))+ #y-axis labeling
  geom_hline(aes(yintercept = RfD.EPA.log,linetype = "EPA RfD"),color = "red")+ #Dashline for EPA RfD
  geom_hline(aes(yintercept = RfD.EFSA.2018.log,linetype = "EFSA TDI,2018"),color = "blue")+ #Dashline for EFSA RfD
  geom_hline(aes(yintercept = RfD.EFSA.2020.log,linetype = "EFSA TDI,2020"),color = "purple")+ #Dashline for EFSA RfD
  scale_linetype_manual(name = "Limit", values = c(2,2,2),guide = guide_legend(override.aes = list(color = c("blue","purple","red"))))+ #the legend of dash lines
  theme_bw()+ #white background
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(),
        axis.title.x = element_blank())+
  theme(axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 12),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 10, family = "serif"),
        axis.text.x = element_text(size = 10, family = "serif"))

 
#######################################CL plot###################################################
RfD.sig.CL.h <- RfD.sig.h %>% filter(Compartment == "CL")
RfD.duration.CL.h <- RfD.duration.h %>% filter(Compartment == "CL")
RfD.matched.CL.h <- RfD.matched.h %>% filter(Compartment == "CL")
RfD.sig.CL.m <- RfD.sig.m %>% filter(compartment == "CL")
RfD.duration.CL.m <- RfD.duration.m %>% filter(compartment == "CL")

RfD.sig.CL.h.log <- RfD.sig.CL.h %>% mutate(RfD.log = log(RfD,10))
RfD.duration.CL.h.log <- RfD.duration.CL.h %>% mutate(RfD.log = log(RfD,10))
RfD.matched.CL.h.log <- RfD.matched.CL.h %>% mutate(RfD.log = log(RfD,10))
RfD.sig.CL.m.log <- RfD.sig.CL.m %>% mutate(RfD.log = log(RfD,10))
RfD.duration.CL.m.log <- RfD.duration.CL.m %>% mutate(RfD.log = log(RfD,10))

colnames(RfD.sig.CL.m.log)[4] <- "Compartment"
RfD.duration.CL.m.log <- RfD.duration.CL.m.log %>% dplyr::select(Duration,compartment,TD, RfD, RfD.log)
colnames(RfD.duration.CL.m.log) <- colnames(RfD.duration.CL.h.log)


df.sig.CL <- rbind(RfD.sig.CL.h.log,RfD.sig.CL.m.log)
df.duration.CL <- rbind(RfD.duration.CL.h.log,RfD.duration.CL.m.log)
df.matched.CL <- rbind(RfD.matched.CL.h.log, RfD.sig.CL.m.log)

plot.CL <- ggplot()+
  geom_point(data = df.sig.CL, aes(x = Duration, y = RfD.log,fill = "Pathway RfD"), size = 2.5, show.legend = NA, color = "tomato",pch = 16 )+ #Pathway RfD
  scale_fill_manual(values = c("tomato","yellow"),name = "RfD",labels = c("Pathway RfD","Duration RfD"))+ #the color of the pathway RfD and duration RfD
  scale_x_discrete(limits = c("Human: Day 1", "Human: Day 4", "Human: Day 10", "Human: Day 14","Mouse: Day 7"))+ #the order of category of the x-axis
  geom_point(data = df.duration.CL,aes(x = Duration, y = RfD.log, fill= "yellow"), pch = 23, size =3, color = "Black", show.legend = NA)+ #duration RfD
  geom_point(data = df.matched.CL, aes(x = Duration, y = RfD.log, shape = Pathway),size = 1.5, color = "lightseagreen",position = position_dodge(0.3,preserve = 'single'))+ #matched RfD
  scale_shape_manual(values = 0:length(unique(df.matched.CL$Pathway)))+ #the shapes of matched RfD
  ylab(expression(paste(log[10],"RfD")))+ #y-axis labeling
  geom_hline(aes(yintercept = RfD.EPA.log,linetype = "EPA RfD"),color = "red")+ #Dashline for EPA RfD
  geom_hline(aes(yintercept = RfD.EFSA.2018.log,linetype = "EFSA TDI, 2018"),color = "blue")+ #Dashline for EFSA RfD
  geom_hline(aes(yintercept = RfD.EFSA.2020.log,linetype = "EFSA TDI, 2020"),color = "purple")+ #Dashline for EFSA RfD
  scale_linetype_manual(name = "Limit", values = c(2,2,2),guide = guide_legend(override.aes = list(color = c("blue","purple","red"))))+ #the legend of dash lines
  theme_bw()+ #white background
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line())+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "serif",size = 12),
        axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 12),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 10, family = "serif"),
        axis.text.x = element_text(size = 10, family = "serif"))

legend <- get_legend(plot.CA + theme(legend.box.margin = margin(0,0,0,12)))
plot.RfD <- plot_grid(plot.CA+rremove("legend"), plot.CL+rremove("legend"),
                       nrow = 2, ncol = 1, align = c("hv"),
                       labels = c("A.","B."), label_fontfamily = "serif", rel_widths = 1 )
plot.RfD.parallel <- plot_grid(plot.RfD,legend,nrow = 1, ncol =2,rel_widths = c(4,1.5))



ggsave("Figure 5_0715.tiff",scale = 1,
       plot = plot.RfD.parallel,
       width = 20, height = 15, units = "cm",dpi=320)

