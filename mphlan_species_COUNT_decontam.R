##decontam of species level count table from metaphlan run
library(phyloseq)
library(tidyverse)
library(readxl)
library(decontam)

##import the  table
setwd('/Users/zachlewis/Desktop/Projects/Urine_Shotgun/urine_shotgun_mphlan')
mphlan_species_counts <- read_tsv("mphlan_counts_species.txt") %>%
  column_to_rownames("clade_name")
#import metadata 
metadata <- read_csv("/Users/zachlewis/Desktop/Projects/Urine_Shotgun/Stats_Files/Urine_Shotgun_Metadata.V3.csv") %>%
  filter(Sample != "NA")

#import to phyloseq
metadata <- metadata %>%
  column_to_rownames("Sample")
sampleData = sample_data(metadata)
OTU = otu_table(mphlan_species_counts, taxa_are_rows = TRUE)
mphlan_species_decontam <- phyloseq(OTU, sampleData)

##visualizing library size ##this matters somewhat less because we're using prevalence based filtering
df <- as.data.frame(sample_data(mphlan_species_decontam)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(mphlan_species_decontam) #generate a library size column
df <- df[order(df$LibrarySize),] #order by library size
df$Index <- seq(nrow(df)) #set an x value 1 per row, in library size order to plot easier
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + 
  geom_point() + scale_y_continuous(limits=c(0,50000000)) ##graph it
#looks meh. never done this with a mphlan_species table though

#first, create a logical vector identifying negatives:
sample_data(mphlan_species_decontam)$is.neg <- sample_data(mphlan_species_decontam)$Sample_or_Control == "Control Sample"

#now run the contaminant identifier
isContam.prev05 <- isContaminant(mphlan_species_decontam, method ='prevalence', neg='is.neg', threshold=0.5)
head(isContam.prev05)
table(isContam.prev05$contaminant)

##36 mphlan_speciess identified as contaminants. could it also be cross-contam? 

##make a phyloseq object of presence/absence in neg cotrols and true samples
mphlan_species_decontam.pa <- transform_sample_counts(mphlan_species_decontam, function(abund) 1*(abund>0)) ##transforming by abundance >0
mphlan_species_decontam.pa.neg <- prune_samples(sample_data(mphlan_species_decontam.pa)$Sample_or_Control == "Control Sample", mphlan_species_decontam.pa) #grab negs
mphlan_species_decontam.pa.pos <- prune_samples(sample_data(mphlan_species_decontam.pa)$Sample_or_Control == "True Sample", mphlan_species_decontam.pa) #grab pos'

# Make data.frame and plot of prevalence in positive and negative samples
df.mphlan_species_decontam.pa <- data.frame(mphlan_species_decontam.pos=taxa_sums(mphlan_species_decontam.pa.pos), mphlan_species_decontam.neg=taxa_sums(mphlan_species_decontam.pa.neg),
                                            contaminant=isContam.prev05$contaminant)
ggplot(data=df.mphlan_species_decontam.pa, aes(x=mphlan_species_decontam.neg, y=mphlan_species_decontam.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

##this is a pretty good looking plot, actually

##export csv otu tables: decontammed table (1), and table of contaminants (2)

#(1)

urine_shotgun_mphlan_species_abund.decontammed <- prune_taxa(!isContam.prev05$contaminant, mphlan_species_decontam)

#convert to non-phyloseq artifact

mphlan_species_decontam.decontammed.table <- as(otu_table(urine_shotgun_mphlan_species_abund.decontammed, taxa_are_rows = TRUE), "matrix")
write.csv(mphlan_species_decontam.decontammed.table, 
          file = '//Users/zachlewis/Desktop/Projects/Urine_Shotgun/urine_shotgun_mphlan/urine_shotgun_mphlan_species_decontammed.csv')


#(2)
mphlan_species_decontam.contaminants <- prune_taxa(isContam.prev05$contaminant, mphlan_species_decontam)
mphlan_species_decontam.contaminants.table <- as(otu_table(mphlan_species_decontam.contaminants, taxa_are_rows=TRUE), "matrix")
write.csv(mphlan_species_decontam.contaminants.table, 
          file = "/Users/zachlewis/Desktop/Projects/Urine_Shotgun/urine_shotgun_mphlan/urine_shotgun_mphlan_species_contaminants.csv")

###manually removed unflagged taxa if they are part of the blacklist AND show up in negatives. 
##manually removing: 
s__Bradyrhizobium_diazoefficiens
s__Bradyrhizobium_viridifuturi
s__Sphingomonas_echinoides
Acinetobacter ursingii
Herbaspirillum huttiense
s__Xanthomonas_massiliensis

##adding back known zymo taxa: 
s__Limosilactobacillus_fermentum
s__Bifidobacterium_adolescentis

