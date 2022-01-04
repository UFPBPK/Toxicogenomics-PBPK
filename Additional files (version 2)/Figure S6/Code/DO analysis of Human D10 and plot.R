library(DOSE)
library(plyr)
library(dplyr)
library(clusterProfiler)
library("org.Hs.eg.db")
library(cowplot)
library(ggplot2)
library(ggpubr)

#BBMD result D10 is in the folder of Figure S1/Results
BBMD.D10 <- read.csv("BBMD result D10.csv", header = TRUE)

BBMD.D10.2 <- BBMD.D10 %>% mutate(BMDU.BMDL = bmdu/bmdl)

BBMD.D10.sig <- BBMD.D10.2 %>% filter(BMDU.BMDL < 40 & bmd < 20)

# DO analysis of Human D10
probe.D10 <-ldply(strsplit(BBMD.D10.sig$ï..Probe.ID,"_"))[,1]
gene.entrezid.D10 <- bitr(probe.D10,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")[,2]

DO.D10 <- enrichDO(gene.entrezid.D10, ont = "DO",pvalueCutoff=0.05, qvalueCutoff = 0.05)
Result.D10<-DO.D10@result
Result.sig.D10 <- Result.D10 %>% filter(p.adjust <0.05 & qvalue <0.05)

# Save DO results
saveRDS(Result.D10, file = "DO human D10.RDS")

#Plot of human D10 data
plot.D10 <- dotplot(DO.D10, showCategory = nrow(Result.sig.D10))+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "serif",size = 15),
        axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 15),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 9, family = "serif"),
        axis.text.x = element_text(size = 10,family = "serif"))

ggsave("Figure S6.tiff", scale = 1,
       plot = plot.D10,
       width = 10, height = 20, dpi = 320)
