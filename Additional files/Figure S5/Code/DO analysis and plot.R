library(DOSE)
library(plyr)
library(dplyr)
library(clusterProfiler)
library("org.Hs.eg.db")
library(cowplot)
library(ggplot2)
library(ggpubr)

#BBMD results are located in the folder of Figure 2/Results
BBMD.D1 <- read.csv("BBMD result D1.csv",header = TRUE)
BBMD.D4 <- read.csv("BBMD result D4.csv", header = TRUE)
BBMD.D10 <- read.csv("BBMD result D10.csv", header = TRUE)
BBMD.D14 <- read.csv("BBMD result D14.csv", header = TRUE)

BBMD.D1.2 <- BBMD.D1 %>% mutate(BMDU.BMDL = bmdu/bmdl)
BBMD.D4.2 <- BBMD.D4 %>% mutate(BMDU.BMDL = bmdu/bmdl)
BBMD.D10.2 <- BBMD.D10 %>% mutate(BMDU.BMDL = bmdu/bmdl)
BBMD.D14.2 <- BBMD.D14 %>% mutate(BMDU.BMDL = bmdu/bmdl)

BBMD.D1.sig <- BBMD.D1.2 %>% filter(BMDU.BMDL < 40 & bmd < 20)
BBMD.D4.sig <- BBMD.D4.2 %>% filter(BMDU.BMDL < 40 & bmd < 20)
BBMD.D10.sig <- BBMD.D10.2 %>% filter(BMDU.BMDL < 40 & bmd < 20)
BBMD.D14.sig <- BBMD.D14.2 %>% filter(BMDU.BMDL < 40 & bmd < 20)

#DO analysis of Human D1
probe.D1 <-ldply(strsplit(BBMD.D1.sig$ï..Probe.ID,"_"))[,1]
gene.entrezid.D1 <- bitr(probe.D1,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")[,2]

DO.D1 <- enrichDO(gene.entrezid.D1,ont = "DO",pvalueCutoff=0.05, qvalueCutoff = 0.05)
Result.D1<-DO.D1@result
Result.sig.D1 <- Result.D1 %>% filter(p.adjust <0.05 & qvalue <0.05)

# DO analysis of Human D4
probe.D4 <-ldply(strsplit(BBMD.D4.sig$ï..Probe.ID,"_"))[,1]
gene.entrezid.D4 <- bitr(probe.D4,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")[,2]

DO.D4 <- enrichDO(gene.entrezid.D4, ont = "DO",pvalueCutoff=0.05, qvalueCutoff = 0.05)
Result.D4<-DO.D4@result
Result.sig.D4 <- Result.D4 %>% filter(p.adjust <0.05 & qvalue <0.05)

# DO analysis of Human D10
probe.D10 <-ldply(strsplit(BBMD.D10.sig$ï..Probe.ID,"_"))[,1]
gene.entrezid.D10 <- bitr(probe.D10,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")[,2]

DO.D10 <- enrichDO(gene.entrezid.D10, ont = "DO",pvalueCutoff=0.05, qvalueCutoff = 0.05)
Result.D10<-DO.D10@result
Result.sig.D10 <- Result.D10 %>% filter(p.adjust <0.05 & qvalue <0.05)

# DO analysis of Human D14
probe.D14 <-ldply(strsplit(BBMD.D14.sig$ï..Probe.ID,"_"))[,1]
gene.entrezid.D14 <- bitr(probe.D14,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")[,2]

DO.D14 <- enrichDO(gene.entrezid.D14, ont = "DO",pvalueCutoff=0.05, qvalueCutoff = 0.05)
Result.D14<-DO.D14@result
Result.sig.D14 <- Result.D14 %>% filter(p.adjust < 0.05 & qvalue < 0.05)

# Save DO results
saveRDS(Result.D1, file = "DO human D1.RDS")
saveRDS(Result.D4, file = "DO human D4.RDS")
saveRDS(Result.D10, file = "DO human D10.RDS")
saveRDS(Result.D14, file = "DO human D14.RDS")


topD1.plot <- dotplot(DO.D1,showCategory = 20)+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "serif",size = 15),
        axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 15),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 10, family = "serif"),
        axis.text.x = element_text(size = 10,family = "serif"))

topD4.plot <- dotplot(DO.D4,showCategory = 20)+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "serif",size = 15),
        axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 15),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 10, family = "serif"),
        axis.text.x = element_text(size = 10,family = "serif"))

topD10.plot <- dotplot(DO.D10,showCategory = 20)+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "serif",size = 15),
        axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 15),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 10, family = "serif"),
        axis.text.x = element_text(size = 10,family = "serif"))

topD14.plot <- dotplot(DO.D14,showCategory = 20)+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "serif",size = 15),
        axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 15),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 10, family = "serif"),
        axis.text.x = element_text(size = 10,family = "serif"))

grobs <- ggplotGrob(topD14.plot)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

DO.parellel <- plot_grid(topD1.plot+rremove("legend"),topD4.plot+rremove("legend"),
                         topD10.plot+rremove("legend"),topD14.plot+rremove("legend"),
                         ncol = 2, align = c("hv"),rel_widths = c(3,3),labels = c("A.","B.","C.","D."),
                         label_fontfamily = "serif")
DO.parellel <-plot_grid(DO.parellel,legend,ncol=2,rel_widths = c(1,.1))

DO.parellel

ggsave("Figure S5.tiff", scale = 1,
       plot = DO.parellel,
       width = 15, height = 10, dpi = 320)
