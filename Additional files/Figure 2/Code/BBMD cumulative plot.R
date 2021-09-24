library(jcolors)
library(ggplot2)
library(cowplot)
library(dplyr)

#The BBMD results are located in the folder of Figure 2/Results
BBMD.D1  <- read.csv("BBMD result D1.csv", header = TRUE)
BBMD.D4  <- read.csv("BBMD result D4.csv", header = TRUE)
BBMD.D10 <- read.csv("BBMD result D10.csv", header = TRUE)
BBMD.D14 <- read.csv("BBMD result D14.csv", header = TRUE)
BBMD.m   <- read.csv("BBMD result mice.csv", header = TRUE)

BBMD.D1.sig <- BBMD.D1 %>% mutate(BMDU.BMDL = bmdu/bmdl) %>% filter(bmd <20 & BMDU.BMDL < 40)
BBMD.D4.sig <- BBMD.D4 %>% mutate(BMDU.BMDL = bmdu/bmdl) %>% filter(bmd <20 & BMDU.BMDL < 40)
BBMD.D10.sig <- BBMD.D10 %>% mutate(BMDU.BMDL = bmdu/bmdl) %>% filter(bmd <20 & BMDU.BMDL < 40)
BBMD.D14.sig <- BBMD.D14 %>% mutate(BMDU.BMDL = bmdu/bmdl) %>% filter(bmd <20 & BMDU.BMDL < 40)
BBMD.m.sig <- BBMD.m  %>% filter(bmd <10 & BMDU.BMDL < 40)

BBMD.D1.2 <- BBMD.D1.sig[order(BBMD.D1.sig$bmd),]
BBMD.D4.2 <- BBMD.D4.sig[order(BBMD.D4.sig$bmd),]
BBMD.D10.2 <- BBMD.D10.sig[order(BBMD.D10.sig$bmd),]
BBMD.D14.2 <- BBMD.D14.sig[order(BBMD.D14.sig$bmd),]
BBMD.m.2 <- BBMD.m.sig[order(BBMD.m.sig$bmd),]

median.D1 <- sapply(seq_along(BBMD.D1.2$bmd),function(x) median(BBMD.D1.2$bmd[1:x]))
median.D4 <- sapply(seq_along(BBMD.D4.2$bmd),function(x) median(BBMD.D4.2$bmd[1:x]))
median.D10 <- sapply(seq_along(BBMD.D10.2$bmd),function(x) median(BBMD.D10.2$bmd[1:x]))
median.D14 <- sapply(seq_along(BBMD.D14.2$bmd),function(x) median(BBMD.D14.2$bmd[1:x]))
median.m  <- sapply(seq_along(BBMD.m.2$bmd),function(x) median(BBMD.m.2$bmd[1:x]))

df.D1 <- cbind.data.frame(index = seq_along(median.D1), Median = median.D1, Duration = "Human: Day 1")
df.D4 <- cbind.data.frame(index = seq_along(median.D4), Median = median.D4, Duration = "Human: Day 4")
df.D10 <- cbind.data.frame(index = seq_along(median.D10), Median = median.D10, Duration = "Human: Day 10")
df.D14 <- cbind.data.frame(index = seq_along(median.D14), Median = median.D14, Duration = "Human: Day 14")
df.m <- cbind.data.frame(index = seq_along(median.m), Median = median.m, Duration = "Mouse: Day 7")

df.median.h <- rbind.data.frame(df.D1, df.D4, df.D10, df.D14)
df.median.h$Duration <- factor(df.median.h$Duration, levels = c("Human: Day 1","Human: Day 4", "Human: Day 10", "Human: Day 14"))

plot.h <- ggplot()+
  geom_point(df.median.h,mapping = aes(x = Median,y = index, color = Duration),size = 1, pch = 21)+
  theme_bw()+
  scale_color_jcolors(palette = "pal2")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line())+
  xlab(expression(paste("BMD Median (",mu,"M)")))+
  ylab("Gene Acumulation (#)")+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "serif",size = 12),
        axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 12),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 10, family = "serif"),
        axis.text.x = element_text(size = 10, family = "serif"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  scale_x_continuous(expand = expansion(mult = c(0,0.1)))+
  theme(panel.spacing=unit(2, "lines"))

plot.m <- ggplot()+
  geom_point(df.m,mapping = aes(x = Median, y = index, fill = Duration),colour = "darkgreen", size = 1, pch = 21)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line())+
  xlab(c("BMD Median (mg/kg)"))+
  ylab("Gene Acumulation (#)")+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "serif",size = 12),
        axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 12),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 10, family = "serif"),
        axis.text.x = element_text(size = 10, family = "serif"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  scale_x_continuous(expand = expansion(mult = c(0,0.1)))+
  scale_fill_discrete(name = "", labels = c("Mouse: Day 7"))+
  theme(panel.spacing=unit(2, "lines"))
  
BBMD.parellel <- plot_grid(plot.h, plot.m, nrow = 2,ncol = 1, align = c("hv"),
                          labels = "AUTO", label_fontfamily = "serif", label_size = 12)

ggsave("Fig.2.tiff",scale = 1,
       plot = BBMD.parellel,
       width = 20, height = 18, units = "cm",dpi=320)
