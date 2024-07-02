##phyloseq for beta diversity on urine shotgun data, singleM
library(vegan)
library(tidyverse)
library(corrplot)
library(ggplot2)
library(dplyr)
library(phyloseq)

##import metadata
setwd("/Users/zachlewis/Desktop/Projects/Urine_Shotgun/Stats_Files")
metadata <- read_csv("Urine_Shotgun_Metadata.V8.csv")

##import table and format long to generate cpm
table <- read_csv("urine_shotgun_singlem_decontammed.csv") %>%
  column_to_rownames("...1")

#prepare for physeq
metadata_phy <- metadata %>%
  column_to_rownames("Sample")
  
#generate physeq object
OTU = otu_table(table, taxa_are_rows = TRUE)
sampledata = sample_data(metadata_phy)
singleM = phyloseq(OTU, sampledata)  

singleM_norm <- transform_sample_counts(singleM, function(x) 100*x/sum(x))
  
bray.ord <- ordinate(singleM_norm, method = "PCoA", distance = "bray")

plot_ordination(singleM_norm, bray.ord, type='sample', 
                shape = 'Kit', color = 'Subject')

bray.matrix <- distance(singleM_norm, "bray") %>%
  as.matrix()

singleM_spiked_dogs <- singleM_norm %>%
  subset_samples(SampleType == "Urine") %>%
  subset_samples(Spiked == "Y")


singleM_dogs_ord <- ordinate(singleM_spiked_dogs, 
                     method = "PCoA", distance = "bray")
plot_ordination(singleM_spiked_dogs, singleM_dogs_ord, type='sample', 
                shape = 'Kit', color = 'Subject')

spiked_matrix <- distance(singleM_spiked_dogs, "bray") %>%
  as.matrix
spiked_metadata <- as(sample_data(singleM_spiked_dogs), "data.frame")

adonis2(spiked_matrix ~ Subject, data = spiked_metadata)
adonis2(spiked_matrix ~ Kit, data = spiked_metadata)
##dog: 0.001**
#kit: 0.1

##merging for taxa barplot
nonnorm_spiked_dog <- singleM %>%
  subset_samples(SampleType == "Urine") %>%
  subset_samples(Spiked == "Y")
merged_kit_singleM <- merge_samples(nonnorm_spiked_dog, "Kit")
Top30 <- names(sort(taxa_sums(merged_kit_singleM), TRUE)[1:30])
urineT30 <- prune_taxa(Top30, merged_kit_singleM)

merged_table <- as(otu_table(urineT30), "matrix") %>%
  as.data.frame()
write_csv(merged_table, file = "singleM_T30.csv")



  