library(DOSE)
library(plyr)
library(dplyr)
library(clusterProfiler)
library("org.Hs.eg.db")
library(cowplot)
library(ggplot2)
library(ggpubr)

# BBMD result D1 is in the folder of Figure S1/Results
BBMD.D1 <- read.csv("BBMD result D1.csv",header = TRUE)

BBMD.D1.2 <- BBMD.D1 %>% mutate(BMDU.BMDL = bmdu/bmdl)

BBMD.D1.sig <- BBMD.D1.2 %>% filter(BMDU.BMDL < 40 & bmd < 20)

#DO analysis of Human D1
probe.D1 <-ldply(strsplit(BBMD.D1.sig$ï..Probe.ID,"_"))[,1]
gene.entrezid.D1 <- bitr(probe.D1,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")[,2]

DO.D1 <- enrichDO(gene.entrezid.D1,ont = "DO",pvalueCutoff=0.05, qvalueCutoff = 0.05)
Result.D1<-DO.D1@result
Result.sig.D1 <- Result.D1 %>% filter(p.adjust <0.05 & qvalue <0.05)


# Save DO results
saveRDS(Result.D1, file = "DO human D1.RDS")


#Plot of human D1 data

plot.D1 <- dotplot(DO.D1, showCategory = nrow(Result.sig.D1))+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "serif",size = 15),
        axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 15),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 10, family = "serif"),
        axis.text.x = element_text(size = 10,family = "serif"))

plot.D1

ggsave("Figure S4.tiff",scale = 1,
       plot = plot.D1,
       width = 10, height = 8,dpi = 320)

