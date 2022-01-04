library(dplyr)
library(ggplot2)

#Loading the HED results for human and mouse,the location is the folder of Figure 2/Results

HED.duration.h <- readRDS(file ="Human duration HED.RDS")
HED.matched.h <- readRDS(file = "Human HED for matched pathways.RDS")

HED.duration.m <- readRDS(file = "Mouse duration HED.RDS")
HED.matched.m <- readRDS(file = "Mouse HED for pathways.RDS")

colnames(HED.duration.m) <- c("Duration","min","low","mid","top","max","Compartment")

df.duration.HED <- rbind(HED.duration.h, HED.duration.m)
df.duration.HED[df.duration.HED == "Day 1"] <- "Human: Day 1"
df.duration.HED[df.duration.HED == "Day 4"] <- "Human: Day 4"
df.duration.HED[df.duration.HED == "Day 10"] <- "Human: Day 10"
df.duration.HED[df.duration.HED == "Day 14"] <- "Human: Day 14"
df.duration.HED[df.duration.HED == "All"] <- "Mouse: Day 7"

HED.matched.m <- HED.matched.m %>% dplyr::select(Pathway,compartment,mid)
HED.matched.m <- cbind(Duration = "Mouse: Day 7", HED.matched.m)

colnames(HED.matched.m)[3:4] <- c("Compartment","HED")

df.matched.HED <- rbind(HED.matched.h,HED.matched.m)
df.matched.HED[df.matched.HED == "Day 1"] <- "Human: Day 1"
df.matched.HED[df.matched.HED == "Day 4"] <- "Human: Day 4"
df.matched.HED[df.matched.HED == "Day 10"] <- "Human: Day 10"
df.matched.HED[df.matched.HED == "Day 14"] <- "Human: Day 14"

df.duration.HED[,2:6]<- apply(df.duration.HED[,2:6],2,as.numeric)
df.matched.HED$HED <- as.numeric(df.matched.HED$HED)

df.duration.HED[,2:6] <- apply(df.duration.HED[,2:6],2,function(x) log((x*1e6),10))
df.matched.HED$HED <- log((df.matched.HED$HED*1e6),10)

plot.HED <- ggplot()+
  geom_boxplot(stat = "identity", data  = df.duration.HED,
               aes(x = Duration, ymin = min, lower = low, middle = mid, upper = top, ymax = max,fill = Compartment),
               position = position_dodge(width = 0.8), width = 0.6, alpha = 0.3)+
  ylab(expression(paste(log[10], "[HED (ng/kg/day)]")))+#y-axis labeling
  scale_x_discrete(limits = c("Human: Day 1", "Human: Day 4", "Human: Day 10", "Human: Day 14","Mouse: Day 7"))+ #the order of category of the x-axis
  geom_point(data = df.matched.HED, # points of HEDs for matched pathways
             aes(x = Duration, y = HED, fill = Compartment,shape = Pathway),size =1.5, 
             position = position_jitterdodge(jitter.width = 1e-6))+
  scale_shape_manual(values = 12:18)+
  scale_fill_discrete(name = "Compartment",labels = c("Plasma","Liver"))+
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

ggsave("Figure 2.tiff",scale = 1,
       plot = plot.HED,
       width = 20, height = 10, units = "cm",dpi=320)

