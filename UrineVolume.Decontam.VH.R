#HowLow CUP samples

Urine.data <- featuretable.filtered.table2.from_biom.forR

MAT <- as.matrix(Urine.data)
# Calling is.matrix() function
is.matrix(MAT)

library(phyloseq); packageVersion("phyloseq")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

browseVignettes("decontam")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")

BiocManager::install("phyloseq")
browseVignettes("phyloseq")
library(phyloseq); packageVersion("phyloseq")

matrix_for_otu <- MAT
class(matrix_for_otu) <- "numeric"

library("phyloseq")

OTU = otu_table(MAT, taxa_are_rows = FALSE)
head(otu_table(OTU))

sampledata = sample_data(decontam.metadata.urine)
sampledata

---------------------------
  
  ps <- phyloseq(OTU, sampledata)
ps

sample_names(sampledata)

get_variable(ps, "Sample_or_Control")

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point() + scale_y_continuous(limits=c(0,45000)) 

sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

# 1569 taxa are not contam, 5 taxa are contam including, 103 129 144 199 940

head(which(contamdf.prev$contaminant))
write.csv(contamdf.prev, file = "export_csv")

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

# 1536 taxa are not contam, 38 taxa are contam 

write.csv(contamdf.prev05, file = "export.csv")

#remove taxa
ps.noncontam3 <- prune_taxa(!contamdf.prev05$contaminant, ps)
ps.noncontam3

# Extract abundance matrix from the phyloseq object (updated feature table)
OTU2 = as(otu_table(ps.noncontam3), "matrix")
# Coerce to data.frame
OTUdf = as.data.frame(OTU2)

write.csv(OTUdf, file = "export.csv")

?isContaminant


----------------------------------------------------------------------------------------------------
 