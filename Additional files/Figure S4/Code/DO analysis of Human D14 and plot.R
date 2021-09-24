library(DOSE)
library(plyr)
library(dplyr)
library(clusterProfiler)
library("org.Hs.eg.db")
library(cowplot)
library(ggplot2)
library(ggpubr)

#BBMD result D14 is in the folder of Figure 2/Results
BBMD.D14 <- read.csv("BBMD result D14.csv", header = TRUE)

BBMD.D14.2 <- BBMD.D14 %>% mutate(BMDU.BMDL = bmdu/bmdl)

BBMD.D14.sig <- BBMD.D14.2 %>% filter(BMDU.BMDL < 40 & bmd < 20)

# DO analysis of Human D14
probe.D14 <-ldply(strsplit(BBMD.D14.sig$ï..Probe.ID,"_"))[,1]
gene.entrezid.D14 <- bitr(probe.D14,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")[,2]

DO.D14 <- enrichDO(gene.entrezid.D14, ont = "DO",pvalueCutoff=0.05, qvalueCutoff = 0.05)
Result.D14<-DO.D14@result
Result.sig.D14 <- Result.D14 %>% filter(p.adjust < 0.05 & qvalue < 0.05)

# Save DO results
saveRDS(Result.D14, file = "DO human D14.RDS")

#Plot of human D14 data
plot.D14 <- dotplot(DO.D14, showCategory = nrow(Result.sig.D14))+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "serif",size = 15),
        axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 15),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 9, family = "serif"),
        axis.text.x = element_text(size = 10,family = "serif"))

plot.D14

ggsave("Figure S4.tiff", scale = 1,
       plot = plot.D14,
       width = 10, height = 30, dpi = 320)
