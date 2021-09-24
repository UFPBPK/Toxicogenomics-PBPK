library(dplyr)
library(ggplot2)

#The data files are located in the folder of Figure S6/Results
df.HED <- readRDS(file = "HED for DO.RDS")
df.duration <- readRDS(file  ="HED duration.RDS")

df.HED <-df.HED %>% dplyr::select(Duration, DO, mid, Compartment)
df.HED[df.HED == "CA"] <- "Plasma"
df.HED[df.HED == "CL"] <- "Liver"

df.HED$mid <-as.numeric(df.HED$mid)
df.HED <- df.HED %>% mutate(HED.log = log(mid*1e6,10))

df.duration[,2:6] <- apply(df.duration[,2:6],2,as.numeric)
df.duration <- df.duration %>% mutate(min.log = log(min*1e6,10),
                                      low.log = log(low*1e6,10),
                                      mid.log = log(mid*1e6,10),
                                      top.log = log(top*1e6,10),
                                      max.log = log(max*1e6,10),)
plot.HED.DO <- ggplot()+
  geom_boxplot(stat = "identity", data = df.duration,
               aes(x = Duration, ymin = min.log, lower = low.log, middle = mid.log, upper = top.log, ymax = max.log,fill = Compartment),
               alpha = 0.3,width = 0.6)+
  ylab(expression(paste(log[10],"HED")))+ #y-axis labeling
  geom_point(data = df.HED,aes(x = Duration, y = HED.log,color = Compartment),
             position = position_jitterdodge(jitter.width = 1e-6),size = 1)+
  theme_bw()+ #white background
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line())+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "serif",size = 12,face = "plain"),
        axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 12),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 10, family = "serif"),
        axis.text.x = element_text(size = 10, family = "serif"))+
  scale_fill_discrete(name = "Compartment",labels = c("Plasma","Liver"))+
  scale_color_manual(values = c("black","black"),name = "", labels = c("Disease HEDs",""))+
  guides(fill = guide_legend(order = 1),
         color = guide_legend(order = 2))

ggsave("Figure S6.tiff", scale = 1,
       plot = plot.HED.DO,
       width = 12, height = 8, dpi = 320)
