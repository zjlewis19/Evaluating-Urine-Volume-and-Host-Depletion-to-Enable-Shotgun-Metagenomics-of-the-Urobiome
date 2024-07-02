###starting exploration of mphlan species count table using phyloseq ##table has been decontammed and 
##had zymo taxa manually removed from urine samples. 
library(phyloseq)
library(vegan)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(pairwiseAdonis)
library(tidyverse)

#imports
setwd('/Users/zachlewis/Desktop/Projects/Urine_Shotgun/urine_shotgun_mphlan')
mphlan <- read_csv("mphlan_species_manual.csv") %>%
  column_to_rownames("clade_name")
metadata <- read_csv("/Users/zachlewis/Desktop/Projects/Urine_Shotgun/Stats_Files/Urine_Shotgun_Metadata.V3.csv") %>%
  filter(Sample != "NA") %>%
  column_to_rownames("SampleID") %>%
  filter(Kit != "NA") %>%
  filter(Subject != "NA")

taxonomy <- read_tsv("mphlan_species_taxonomy.txt") %>%
  column_to_rownames("...1")
taxmat <- as.matrix(taxonomy)
sampleData = sample_data(metadata)
OTU = otu_table(mphlan, taxa_are_rows = TRUE)
tax <- tax_table(taxmat)
mphlanps <- phyloseq(OTU, sampleData, tax)

mphlanps.norm <- transform_sample_counts(mphlanps, function(x) 100*x/sum(x))

#filter out zero abundance samples (appears to be the 5 bacteremia spiked samples)
mphlanps.nozero <- subset_samples(mphlanps.norm, sample_sums(mphlanps) != 0) ##spared 47 samples
mphlanps.zeros <- subset_samples(mphlanps.norm, sample_sums(mphlanps) == 0) ##21 samples sum zero
##these include: the spiked bacteremia samples, some negatives, a number of PMA-extracted samples

bray.ord <- ordinate(mphlanps.nozero, method = "PCoA", distance = "bray")

##okay, first mphlan plot
plot_ordination(mphlanps.nozero, bray.ord, type='sample', shape = 'Kit', color = 'Subject') +
  theme_bw() + 
  geom_point(size=3) + scale_shape_manual(values = c(15, 16, 17, 18, 8), name = "Extraction Method") +
  ggtitle("all samples, bray")

mphlanps.urine <- subset_samples(mphlanps.nozero, SampleType == "Urine")
bray.urine <- ordinate(mphlanps.urine, method = "PCoA", distance = "bray")
plot_ordination(mphlanps.urine, bray.ord, type='sample', shape = 'Kit', color = 'Subject') +
  theme_bw() + 
  geom_point(size=3) + scale_shape_manual(values = c(15, 16, 17, 18, 8), name = "Extraction Method") +
  ggtitle("all urine samples, bray")

##testing beta diversity
urine.braymatrix <- distance(mphlanps.urine, method = "bray")
urineMetadata <- as(sample_data(mphlanps.urine), "data.frame")
adonis2(urine.braymatrix ~ Subject, data = urineMetadata)
adonis2(urine.braymatrix ~ Kit, data = urineMetadata)
##effect of dog p=0.001
##effect of kit p=0.97

##subsetting spiked samples ##23 samples
mphlanps.urine.spiked <- subset_samples(mphlanps.urine, Spiked == "Y") #26 samples
bray.urine.spiked <- ordinate(mphlanps.urine.spiked, method = "PCoA", distance = "bray")
plot_ordination(mphlanps.urine.spiked, bray.urine.spiked, type='sample', shape = 'Kit', color = 'Subject') +
  theme_bw() + 
  geom_point(size=3) + scale_shape_manual(values = c(15, 16, 17, 18, 8), name = "Extraction Method") +
  ggtitle("spiked urine, bray")
##testing beta diversity
urine.spiked.braymatrix <- distance(mphlanps.urine.spiked, method = "bray")
urineMetadata.spiked <- as(sample_data(mphlanps.urine.spiked), "data.frame")
adonis2(urine.spiked.braymatrix ~ Subject, data = urineMetadata.spiked)
adonis2(urine.spiked.braymatrix ~ Kit, data = urineMetadata.spiked)
##effect of dog p=0.001
##effect of kit p=0.96
urine.spiked.jaccmatrix <- distance(mphlanps.urine.spiked, method = "jaccard", binary = TRUE)
adonis2(urine.spiked.jaccmatrix ~ Subject, data = urineMetadata.spiked)
adonis2(urine.spiked.jaccmatrix ~ Kit, data = urineMetadata.spiked)
jacc.urine.ord <- ordinate(mphlanps.urine.spiked, method = "PCoA", distance = "jaccard", binary = TRUE)
plot_ordination(mphlanps.urine.spiked, jacc.urine.ord, type='sample', shape = 'Kit', color = 'Subject') +
  theme_bw() + 
  geom_point(size=3) + scale_shape_manual(values = c(15, 16, 17, 18, 8), name = "Extraction Method") +
  ggtitle("spiked urine, jacc")


#look at some alpha diversity ##this doesn't work now that it's been normalized
alpha_meas = c("Observed", "Shannon")
p <- plot_richness(mphlanps.urine.spiked, "Kit", "Subject", measures = alpha_meas)
p

barplot(sort(taxa_sums(mphlanps.urine.spiked), TRUE)[1:20]/nsamples(mphlanps.urine),las=2)
Top20 <- names(sort(taxa_sums(mphlanps.urine.spiked), TRUE)[1:20])
urineT20 <- prune_taxa(Top20, mphlanps.urine.spiked)

plot_bar(urineT20, x= "Genus", fill = "Genus")

##generating alpha diversity tables of all samples to export to the metadata

mphlan_alpha <- estimate_richness(mphlanps, split = TRUE, measures = alpha_meas) %>%
  rownames_to_column("SampleID") %>%
  rename(Mphlan_Observed = Observed, Mphlan_Shannon = Shannon)

write_csv(mphlan_alpha, file = "mphlan_alpha_estimates.csv")



