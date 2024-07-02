###making figures for HLCYG - cleaned for submission. No stats performed here. 
library(phyloseq)
library(vegan)
library(ggplot2)
library(readxl)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(pairwiseAdonis)
library(devtools)
library(pairwiseAdonis)
library(tidyverse)

##this is cleaned code used only for generating figures. 
##all figures match patterns observed in depth-adjusted qiime emperor plots. 

setwd('/Users/zachlewis/Desktop/Projects/HLCYG/HLCYG.Reanalysis')
otu_table <- read.csv('HL.decontammed.table.csv') %>%
  column_to_rownames("ASV")

metadata <- read.delim('HL.Dilution.Metadata.Urine.V3.txt') %>%
  column_to_rownames("SampleID")

OTU = otu_table(otu_table, taxa_are_rows=TRUE) ##import to ps
sampleData = sample_data(metadata) ##import to ps
tree = read_tree('rooted.tree.unzip/data/tree.nwk') ##generated in qiime

#import to phyloseq
ps = phyloseq(OTU, sampleData, tree)

##filtering down to relevant samples: 
HLCYG.urine <- subset_samples(ps, KitID != "N") %>%
  subset_samples(Sample_type_negative == "Urine") %>%
  subset_samples(Description != "SJM")

##creating ordinations
wunifrac.ord <- ordinate(HLCYG.urine, method = "PCoA", distance = "wunifrac")
bray.ord <- ordinate(HLCYG.urine, method = "PCoA", distance = "bray")
jacc.ord <- ordinate(HLCYG.urine, method = "PCoA", distance = "jaccard", BINARY = T)
unifrac.ord <- ordinate(HLCYG.urine, method = "PCoA", distance = "unifrac")

##create plots
plot_ordination(HLCYG.urine, wunifrac.ord, type='sample', shape = 'KitID', color = 'Dog') +
  theme_bw() + 
  geom_point(size=3) + scale_shape_manual(values = c(15, 16, 17, 18, 8), name = "Extraction Method") +
  ggtitle("wunifrac")

plot_ordination(HLCYG.urine, bray.ord, type='sample', shape = 'KitID', color = 'Dog') +
  theme_bw() + 
  geom_point(size=3) + scale_shape_manual(values = c(15, 16, 17, 18, 8), name = "Extraction Method") +
  ggtitle("bray")

plot_ordination(HLCYG.urine, jacc.ord, type='sample', shape = 'KitID', color = 'Dog') +
  theme_bw() + 
  geom_point(size=3) + scale_shape_manual(values = c(15, 16, 17, 18, 8), name = "Extraction Method") +
  ggtitle("jacc")

plot_ordination(HLCYG.urine, unifrac.ord, type='sample', shape = 'KitID', color = 'Dog') +
  theme_bw() + 
  geom_point(size=3) + scale_shape_manual(values = c(15, 16, 17, 18, 8), name = "Extraction Method") +
  ggtitle("unifrac")


