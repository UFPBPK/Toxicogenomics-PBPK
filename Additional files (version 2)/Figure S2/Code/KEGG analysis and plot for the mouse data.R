library(AnnotationDbi)
library(clusterProfiler)
library("mouse4302.db") #Affymetrix Mouse Genome 430 2.0 Array annotation data
library(cowplot)
library(dplyr)
library(ggplot2)
library(plyr)
library(ggpubr)

# The BBMD results are located in the folder of Figure S1/Results
BBMD.mice <- read.csv("BBMD result mice.csv",header = TRUE)

BBMD.mice.sig <- BBMD.mice %>% filter(BMDU.BMDL < 40 & bmd < 10)

#KEGG analysis for GSE22871 mice data
OUT <- AnnotationDbi::select(mouse4302.db,keys = BBMD.mice.sig$ï..Probe.ID,columns = c("SYMBOL","ENTREZID"))
OUT <- OUT[!duplicated(OUT$PROBEID),]
kegg.mice <- enrichKEGG(OUT$ENTREZID,organism = "mmu", keyType = "kegg",
                        pvalueCutoff = 0.05, pAdjustMethod = "BH",
                        qvalueCutoff = 0.05, use_internal_data = FALSE)
pathway.mice <- kegg.mice@result
path.sig.mice <- pathway.mice %>% filter(p.adjust < 0.05 & qvalue < 0.05)

#save KEGG results
saveRDS(path.sig.mice, file = "KEGG sig mice.RDS")

#KEGG plots
dot.m <- dotplot(kegg.mice)+
  theme(axis.title.x = element_text(lineheight = 1.5,family = "Serif",size = 15),
        axis.title.y = element_text(lineheight = 1.5,family = "Serif",size = 15),
        legend.title = element_text(lineheight = 1.5,family = "Serif"),
        legend.text  = element_text(lineheight = 1.5,family = "serif"),
        axis.text.y = element_text(size = 8, family = "serif"),
        axis.text.x = element_text(size = 8, family = "serif"))

ggsave("Figure S2.tiff",scale = 1,
       plot = dot.m,
       width = 15, height = 10, units = "cm",dpi=320)