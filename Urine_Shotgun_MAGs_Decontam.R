##decontam on the MAG abundance table for uricomplete. 
library(phyloseq)
library(tidyverse)
library(readxl)
library(decontam)

##import the mag table
setwd('/Users/zachlewis/Desktop/Projects/Urine_Shotgun/Stats_Files')
mag_abund <- read_csv("all_bin_abundance.csv") %>%
  column_to_rownames("Bin")

OTU = otu_table(mag_abund, taxa_are_rows = TRUE)

#import metadata 
metadata <- read_csv("Urine_Shotgun_Metadata.V3.csv") %>%
  filter(Sample != "NA") %>%
  column_to_rownames("Sample")
sampleData = sample_data(metadata)

mag_decontam <- phyloseq(OTU, sampleData)

#not normalizing because this is prevalence-based filtering. 

##visualizing library size; this isn't really necessary for prevalence-based filtering.
df <- as.data.frame(sample_data(mag_decontam)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(mag_decontam) #generate a library size column
df <- df[order(df$LibrarySize),] #order by library size
df$Index <- seq(nrow(df)) #set an x value 1 per row, in library size order to plot easier
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + 
  geom_point() + scale_y_continuous(limits=c(0,100000)) ##graph it
#looks meh. never done this with a mag table though

#first, create a logical vector identifying negatives:
sample_data(mag_decontam)$is.neg <- sample_data(mag_decontam)$Sample_or_Control == "Control Sample"

#now run the contaminant identifier
isContam.prev05 <- isContaminant(mag_decontam, method ='prevalence', neg='is.neg', threshold=0.5)
head(isContam.prev05)
table(isContam.prev05$contaminant)

##11/27 mags identified as contaminants. could it also be cross-contam?

##make a phyloseq object of presence/absence in neg cotrols and true samples
mag_decontam.pa <- transform_sample_counts(mag_decontam, function(abund) 1*(abund>0)) ##transforming by abundance >0
mag_decontam.pa.neg <- prune_samples(sample_data(mag_decontam.pa)$Sample_or_Control == "Control Sample", mag_decontam.pa) #grab negs
mag_decontam.pa.pos <- prune_samples(sample_data(mag_decontam.pa)$Sample_or_Control == "True Sample", mag_decontam.pa) #grab pos'

# Make data.frame and plot of prevalence in positive and negative samples
df.mag_decontam.pa <- data.frame(mag_decontam.pos=taxa_sums(mag_decontam.pa.pos), mag_decontam.neg=taxa_sums(mag_decontam.pa.neg),
                                 contaminant=isContam.prev05$contaminant)
ggplot(data=df.mag_decontam.pa, aes(x=mag_decontam.neg, y=mag_decontam.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

##there are clear contaminants based solely on prevalence but one very prevalent non-flagged
##contaminant appears

##export csv otu tables: decontammed table (1), and table of contaminants (2)

#(1)

urine_shotgun_mag_abund.decontammed <- prune_taxa(!isContam.prev05$contaminant, mag_decontam)

#convert to non-phyloseq artifact

mag_decontam.decontammed.table <- as(otu_table(urine_shotgun_mag_abund.decontammed, taxa_are_rows = TRUE), "matrix")
write.csv(mag_decontam.decontammed.table, 
          file = '/Users/zachlewis/Desktop/Projects/Urine_Shotgun/Stats_Files/urine_shotgun_mag_decontammed.csv')


#(2)
mag_decontam.contaminants <- prune_taxa(isContam.prev05$contaminant, mag_decontam)
mag_decontam.contaminants.table <- as(otu_table(mag_decontam.contaminants, taxa_are_rows=TRUE), "matrix")
write.csv(mag_decontam.contaminants.table, 
          file = "/Users/zachlewis/Desktop/Projects/Urine_Shotgun/Stats_Files/urine_shotgun_mag_contaminants.csv")


