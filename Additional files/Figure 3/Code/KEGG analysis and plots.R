library(AnnotationDbi)
library(clusterProfiler)
library("mouse4302.db") #Affymetrix Mouse Genome 430 2.0 Array annotation data
library("org.Hs.eg.db") #Human annotation 
library(cowplot)
library(dplyr)
library(ggplot2)
library(plyr)
library(ggpubr)

# The BBMD results are located in the folder of Figure 2/Results
BBMD.mice <- read.csv("BBMD result mice.csv",header = TRUE)
BBMD.D1 <- read.csv("BBMD result D1.csv",header = TRUE)
BBMD.D4 <- read.csv("BBMD result D4.csv", header = TRUE)
BBMD.D10 <- read.csv("BBMD result D10.csv", header = TRUE)
BBMD.D14 <- read.csv("BBMD result D14.csv", header = TRUE)

BBMD.D1.2 <- BBMD.D1 %>% mutate(BMDU.BMDL = bmdu/bmdl)
BBMD.D4.2 <- BBMD.D4 %>% mutate(BMDU.BMDL = bmdu/bmdl)
BBMD.D10.2 <- BBMD.D10 %>% mutate(BMDU.BMDL = bmdu/bmdl)
BBMD.D14.2 <- BBMD.D14 %>% mutate(BMDU.BMDL = bmdu/bmdl)

BBMD.mice.sig <- BBMD.mice %>% filter(BMDU.BMDL < 40 & bmd < 10)
BBMD.D1.sig <- BBMD.D1.2 %>% filter(BMDU.BMDL < 40 & bmd < 20)
BBMD.D4.sig <- BBMD.D4.2 %>% filter(BMDU.BMDL < 40 & bmd < 20)
BBMD.D10.sig <- BBMD.D10.2 %>% filter(BMDU.BMDL < 40 & bmd < 20)
BBMD.D14.sig <- BBMD.D14.2 %>% filter(BMDU.BMDL < 40 & bmd < 20)

#KEGG analysis for GSE22871 mice data
OUT <- AnnotationDbi::select(mouse4302.db,keys = BBMD.mice.sig$ï..Probe.ID,columns = c("SYMBOL","ENTREZID"))
OUT <- OUT[!duplicated(OUT$PROBEID),]
kegg.mice <- enrichKEGG(OUT$ENTREZID,organism = "mmu", keyType = "kegg",
                        pvalueCutoff = 0.05, pAdjustMethod = "BH",
                        qvalueCutoff = 0.05, use_internal_data = FALSE)
pathway.mice <- kegg.mice@result
path.sig.mice <- pathway.mice %>% filter(p.adjust < 0.05 & qvalue < 0.05)

#KEGG analysis for GSE144775 human data

#D1
probe.D1 <-ldply(strsplit(BBMD.D1.sig$ï..Probe.ID,"_"))[,1]
gene.entrezid.D1 <- bitr(probe.D1,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")[,2]

kegg.D1<-enrichKEGG(gene.entrezid.D1,organism = "hsa", keyType = "kegg",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH",
                 qvalueCutoff = 0.05, use_internal_data = FALSE)

pathway.D1 <- kegg.D1@result
path.sig.D1 <- pathway.D1 %>% filter(p.adjust < 0.05 & qvalue < 0.05) 

#D4
probe.D4 <-ldply(strsplit(BBMD.D4.sig$ï..Probe.ID,"_"))[,1]
gene.entrezid.D4 <- bitr(probe.D4,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")[,2]

kegg.D4<-enrichKEGG(gene.entrezid.D4,organism = "hsa", keyType = "kegg",
                    pvalueCutoff = 0.05, pAdjustMethod = "BH",
                    qvalueCutoff = 0.05, use_internal_data = FALSE)

pathway.D4 <- kegg.D4@result
path.sig.D4 <- pathway.D4 %>% filter(p.adjust < 0.05 & qvalue < 0.05)

#D10
probe.D10 <-ldply(strsplit(BBMD.D10.sig$ï..Probe.ID,"_"))[,1]
gene.entrezid.D10 <- bitr(probe.D10,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")[,2]

kegg.D10<-enrichKEGG(gene.entrezid.D10,organism = "hsa", keyType = "kegg",
                    pvalueCutoff = 0.05, pAdjustMethod = "BH",
                    qvalueCutoff = 0.05, use_internal_data = FALSE)

pathway.D10 <- kegg.D10@result
path.sig.D10 <- pathway.D10 %>% filter(p.adjust < 0.05 & qvalue < 0.05)

#D14
probe.D14 <-ldply(strsplit(BBMD.D14.sig$ï..Probe.ID,"_"))[,1]
gene.entrezid.D14 <- bitr(probe.D14,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")[,2]

kegg.D14<-enrichKEGG(gene.entrezid.D14,organism = "hsa", keyType = "kegg",
                    pvalueCutoff = 0.05, pAdjustMethod = "BH",
                    qvalueCutoff = 0.05, use_internal_data = FALSE)

pathway.D14 <- kegg.D14@result
path.sig.D14 <- pathway.D14 %>% filter(p.adjust < 0.05 & qvalue < 0.05)

#save KEGG results
saveRDS(path.sig.mice, file = "KEGG sig mice.RDS")
saveRDS(pathway.D1, file = "KEGG human D1.RDS")
saveRDS(pathway.D4, file = "KEGG human D4.RDS")
saveRDS(pathway.D10,file = "KEGG human D10.RDS")
saveRDS(pathway.D14, file = "KEGG human D14.RDS")

#KEGG plots
dot.m <- dotplot(kegg.mice)+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "Serif",size = 15),
        axis.title.y = element_text(lineheight = 1.5,family = "Serif",size = 15),
        legend.title = element_text(lineheight = 1.5,family = "Serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 8, family = "serif"),
        axis.text.x = element_text(size = 8, family = "serif"))

dot.D1 <- dotplot(kegg.D1,showCategory = nrow(path.sig.D1))+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "serif",size = 15),
        axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 15),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 8, family = "serif"),
        axis.text.x = element_text(size = 8, family = "serif"))
dot.D4 <- dotplot(kegg.D4,showCategory = nrow(path.sig.D4))+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "serif",size = 15),
        axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 15),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 8, family = "serif"),
        axis.text.x = element_text(size = 8, family = "serif"))
dot.D10 <- dotplot(kegg.D10,showCategory = nrow(path.sig.D10))+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "serif",size = 15),
        axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 15),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 8, family = "serif"),
        axis.text.x = element_text(size = 8, family = "serif"))
dot.D14 <- dotplot(kegg.D14,showCategory = nrow(path.sig.D14))+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "serif",size = 15),
        axis.title.y = element_text(lineheight = 1.5,family = "serif",size = 15),
        legend.title = element_text(lineheight = 1.5,family = "serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 8, family = "serif"),
        axis.text.x = element_text(size = 8, family = "serif"))

legend <- get_legend(dot.D10 + theme(legend.box.margin = margin(0,0,0,12)))
plot.kegg <- plot_grid(dot.D1+rremove("legend"), dot.D4+rremove("legend"), dot.D10+rremove("legend"), dot.D14+rremove("legend"), dot.m+rremove("legend"),
                       nrow = 3, ncol = 2, align = c("hv"),rel_widths = c(5,5),
                       labels = c("A.","B.","C.","D.","E."), label_fontfamily = "serif" )

plot.legend <- plot_grid(legend,nrow = 1,ncol = 1)
ggsave("Fig.3.tiff",scale = 1,
       plot = plot.kegg,
       width = 30, height = 30, units = "cm",dpi=320)

ggsave("legend.tiff",scale = 1,
       plot = plot.legend,
       width = 10, height = 10, units = "cm",dpi=320)

