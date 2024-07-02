##Analysis of Urine Volume Data - solely figure generation. Stats performed in QIIME on normalized data.
library(phyloseq)
library(vegan)
library(ggplot2)
library(readxl)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(tidyverse)


##all statistical testing here is extraneous; this code was exclusively
##for figure generation, which match patterns seen on qiime emperor plots.
##

setwd('/Users/zachlewis/Desktop/Projects/HLCYG/UrineVolume/Reanalysis')
otu_table <- read.delim('taxa-filtered-table.txt') %>% ##this is a decontammed, taxa-filtered OTU table.
  column_to_rownames("X")
metadata <- read.delim('urinevolume.metadata.rean.v4.txt') %>% ##metadata is available on BioProject.
  column_to_rownames("X.SampleID") %>%
  rename(Dog = dog)

metadata$urine_volume_ml <- as.factor(metadata$urine_volume_ml)

##import to phyloseq
OTU = otu_table(otu_table, taxa_are_rows=TRUE)
sampleData = sample_data(metadata)
ps <- phyloseq(OTU, sampleData)  

##statistics have already been run. ##this is solely for figure-generation. 
##filter negatives and the ArB sample that didn't sequence. 
urineVolume.filtered <- subset_samples(ps, timepoint!="NA") %>%
  subset_samples(PerSampleFrequency > 1)

##make a boxplot of dispersions. 
brayMatrix <- distance(urineVolume.filtered, method = "bray")
disperMeta <- filter(metadata, sample_type != "Negative control") %>%
  filter(BarcodeSequence != "AGAAAGGGTGTG") ##get rid of ArB_3
jaccmatrix <- distance(urineVolume.filtered, method = "jaccard", binary = T)
##make my color palettes
perDogPalette1 <- brewer.pal(n=9, "Blues")[4]
perDogPalette2 <- brewer.pal(n=9, "Blues")[8]
my_palette <- brewer.pal(n=9, "Blues")[3:8]

##permdisp
#bray_dispersions <- betadisper(brayMatrix, disperMeta$Volume_High_Low, type = c("centroid"))
#anova(bray_dispersions)
#boxplot(bray_dispersions, col = c(perDogPalette2, perDogPalette1), xLab = "Volume")

#jacc_dispersions <- betadisper(jaccmatrix, disperMeta$Volume_High_Low, type = c("centroid"))
#anova(jacc_dispersions)
#boxplot(jacc_dispersions, col = c(perDogPalette2, perDogPalette1), xLab = "Volume")

#adonis2(brayMatrix ~ Dog, data = disperMeta)
#adonis2(brayMatrix ~ urine_volume_ml, data = disperMeta)
#adonis2(jaccmatrix ~ Dog, data = disperMeta)
#adonis2(jaccmatrix ~ urine_volume_ml, data = disperMeta)

##plot jaccard and bray curtis. 
#odinations
urine.bray.ord <- ordinate(urineVolume.filtered, method = "PCoA", distance = "bray")
urine.jacc.ord <- ordinate(urineVolume.filtered, method = "PCoA", distance = "jaccard", binary = T)

#plots
plot_ordination(urineVolume.filtered, urine.bray.ord, type='sample', shape = 'Dog', color = 'urine_volume_ml') + 
  scale_color_manual(values = my_palette, name = "Sample Volume (mL)") + theme_bw() + 
  scale_shape_manual(values = c(15, 16, 17, 18, 8)) + ##pick shapes and color!!!! this is a sexy plot
  geom_point(size=4) + theme(legend.text = element_text(size = 18), 
                             legend.title = element_text(size = 20))

plot_ordination(urineVolume.filtered, urine.jacc.ord, type='sample', shape = 'Dog', color = 'urine_volume_ml') + 
  scale_color_manual(values = my_palette, name = "Sample Volume (mL)") + theme_bw() + 
  scale_shape_manual(values = c(15, 16, 17, 18, 8)) + 
  geom_point(size=4) + theme(legend.text = element_text(size = 18), 
                             legend.title = element_text(size = 20))


##generate per-dog ordinations
urineVolumeKH <- subset_samples(urineVolume.filtered, Dog=='KH')
urineVolumeArB <- subset_samples(urineVolume.filtered, Dog=='ArB')
urineVolumeMS <- subset_samples(urineVolume.filtered, Dog=='MS')
urineVolumeFC <- subset_samples(urineVolume.filtered, Dog=='FC')
urineVolumeHF <- subset_samples(urineVolume.filtered, Dog=='HF')
KH.ord <- ordinate(urineVolumeKH, "PCoA", "bray")
ArB.ord <- ordinate(urineVolumeArB, "PCoA", "bray")
MS.ord <- ordinate(urineVolumeMS, "PCoA", "bray")
FC.ord <- ordinate(urineVolumeFC, "PCoA", "bray")
HF.ord <- ordinate(urineVolumeHF, "PCoA", "bray")
KH.jacc.ord <- ordinate(urineVolumeKH, "PCoA", distance = "jaccard", BINARY = T)
ArB.jacc.ord <- ordinate(urineVolumeArB, "PCoA", distance = "jaccard", BINARY = T)
MS.jacc.ord <- ordinate(urineVolumeMS, "PCoA", distance = "jaccard", BINARY = T)
FC.jacc.ord <- ordinate(urineVolumeFC, "PCoA", distance = "jaccard", BINARY = T)
HF.jacc.ord <- ordinate(urineVolumeHF, "PCoA", distance = "jaccard", BINARY = T)

#plot bray
plot_ordination(urineVolumeKH, KH.ord, type='sample', color = 'Volume_High_Low') + 
  scale_color_manual(values = c(perDogPalette2, perDogPalette1), name = "Sample Volume") + theme_bw() +
  ggtitle("Bray Curtis: KH") + geom_point(size=4) 

plot_ordination(urineVolumeArB, ArB.ord, type='sample', color = 'Volume_High_Low') + 
  scale_color_manual(values = c(perDogPalette2, perDogPalette1), name = "Sample Volume") + 
  theme_bw() + ggtitle("Bray Curtis: ArB") + geom_point(size=4)

plot_ordination(urineVolumeMS, MS.ord, type='sample', color = 'Volume_High_Low') + 
  scale_color_manual(values = c(perDogPalette2, perDogPalette1), name = "Sample Volume") + 
  theme_bw() + ggtitle("Bray Curtis: MS") + geom_point(size=4)


plot_ordination(urineVolumeFC, FC.ord, type='sample', color = 'Volume_High_Low') + 
  scale_color_manual(values = c(perDogPalette2, perDogPalette1), name = "Sample Volume") + 
  theme_bw() + ggtitle("Bray Curtis: FC") + geom_point(size=4)

plot_ordination(urineVolumeHF, HF.ord, type='sample', color = 'Volume_High_Low') + 
  scale_color_manual(values = c(perDogPalette2, perDogPalette1), name = "Sample Volume") + 
  theme_bw() + ggtitle("Bray Curtis: HF") + geom_point(size=4)

##plot jaccard
plot_ordination(urineVolumeKH, KH.jacc.ord, type='sample', color = 'Volume_High_Low') + 
  scale_color_manual(values = c(perDogPalette2, perDogPalette1), name = "Sample Volume") + theme_bw() +
  ggtitle("Jaccard: KH") + geom_point(size=4) 

plot_ordination(urineVolumeArB, ArB.jacc.ord, type='sample', color = 'Volume_High_Low') + 
  scale_color_manual(values = c(perDogPalette2, perDogPalette1), name = "Sample Volume") + 
  theme_bw() + ggtitle("Jaccard: ArB") + geom_point(size=4)

plot_ordination(urineVolumeMS, MS.jacc.ord, type='sample', color = 'Volume_High_Low') + 
  scale_color_manual(values = c(perDogPalette2, perDogPalette1), name = "Sample Volume") + 
  theme_bw() + ggtitle("Jaccard: MS") + geom_point(size=4)


plot_ordination(urineVolumeFC, FC.jacc.ord, type='sample', color = 'Volume_High_Low') + 
  scale_color_manual(values = c(perDogPalette2, perDogPalette1), name = "Sample Volume") + 
  theme_bw() + ggtitle("Jaccard: FC") + geom_point(size=4)

plot_ordination(urineVolumeHF, HF.jacc.ord, type='sample', color = 'Volume_High_Low') + 
  scale_color_manual(values = c(perDogPalette2, perDogPalette1), name = "Sample Volume") + 
  theme_bw() + ggtitle("Jaccard: HF") + geom_point(size=4)
