library(DOSE)
library(plyr)
library(dplyr)
library(clusterProfiler)
library("org.Hs.eg.db")
library(cowplot)
library(ggplot2)
library(ggpubr)

#BBMD result D4 is in the folder of Figure 2/Results
BBMD.D4 <- read.csv("BBMD result D4.csv", header = TRUE)

BBMD.D4.2 <- BBMD.D4 %>% mutate(BMDU.BMDL = bmdu/bmdl)

BBMD.D4.sig <- BBMD.D4.2 %>% filter(BMDU.BMDL < 40 & bmd < 20)

# DO analysis of Human D4
probe.D4 <-ldply(strsplit(BBMD.D4.sig$ï..Probe.ID,"_"))[,1]
gene.entrezid.D4 <- bitr(probe.D4,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")[,2]

DO.D4 <- enrichDO(gene.entrezid.D4, ont = "DO",pvalueCutoff=0.05, qvalueCutoff = 0.05)
Result.D4<-DO.D4@result
Result.sig.D4 <- Result.D4 %>% filter(p.adjust <0.05 & qvalue <0.05)


# Save DO results
saveRDS(Result.D4, file = "DO human D4.RDS")

#Plot of human D4 data
plot.D4 <- dotplot(DO.D4, showCategory = nrow(Result.sig.D4))+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "serif",size = 15),
        axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 15),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 10, family = "serif"),
        axis.text.x = element_text(size = 10,family = "serif"))

plot.D4

ggsave("Figure S2.tiff", scale = 1,
       plot = plot.D4,
       width  = 10, height = 20, dpi = 320)

