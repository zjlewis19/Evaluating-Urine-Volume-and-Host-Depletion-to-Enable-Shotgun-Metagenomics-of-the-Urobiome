###Compiling all code from HLCYG in preparation for manuscript submission, and rerunning some analyses
###to increase rigor. 

##this code was run on the Ohio Supercomputer. Host Depletion samples were sequenced alongside another project, 
##and were denoised and clustered together. 

cd /users/PAS1331/lewis3359/Gator-HL.Dilutions
#activate conda and qiime
module activate miniconda3
source activate qiime2-2023.2

qiime dada2 denoise-paired \
--i-demultiplexed-seqs Gator-HL-demux.qza \
--p-trim-left-f 5 \
--p-trim-left-r 5 \
--p-trunc-len-f 250 \
--p-trun-len-r 231 \
--o-table table.qza \
--o-representative-sequences rep-seqs.qza \
--o-denoising-stats denoising-stats.qza

qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file March2023.Combined.Mapping.V1.txt

##switched to local machine. 
##qiime version now 2024.2

##filtering to urine project samples
qiime feature-table filter-samples \
--i-table table.qza \
--m-metadata-file HL.Dilution.Metadata.V2.txt \
--o-filtered-table HL-dilution-table.qza

qiime feature-table summarize \
--i-table HL-dilution-table.qza \
--m-sample-metadata-file HL.Dilution.Metadata.V2.txt \
--o-visualization HL-dilution-table.qzv

##sent to decontam, and re-imported to qiime. 
cd /Users/zachlewis/Desktop/Projects/HLCYG/HLCYG.Reanalysis
qiime tools import \
  --input-path HL.decontammed.table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format \
  --output-path HL-decontammed-table.qza

##filter sequences to only urine sequences, to reduce alignment burden
qiime feature-table filter-seqs \
--i-data rep-seqs.qza \
--i-table HL-decontammed-table.qza \
--o-filtered-data HLCYG-seqs.qza

##generate a tree
qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences HLCYG-seqs.qza \
    --o-alignment alinged-rep-seqs.qza \
    --o-masked-alignment masked-aligned-rep-seqs.qza \
    --o-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree.qza

##generate a taxonomy
qiime feature-classifier classify-sklearn \
  --i-classifier SILVA138.V4.99.qza \
  --i-reads HLCYG-seqs.qza \
  --o-classification HLCYG-taxonomy.qza

qiime metadata tabulate \
--m-input-file HLCYG-taxonomy.qza \
--o-visualization HLCYG-taxonomy.qzv

##filter mito, eukarya, unassigned sequences.
qiime taxa filter-table \
--i-table HL-decontammed-table.qza \
--i-taxonomy HLCYG-taxonomy.qza \
--p-exclude Unassigned,Eukaryota,Mitochondria \
--o-filtered-table HLCYG-taxa-filtered-table.qza

qiime feature-table summarize \
--i-table HLCYG-taxa-filtered-table.qza \
--m-sample-metadata-file HL.Dilution.Metadata.V2.txt \
--o-visualization HL-taxa-filtered-table.qzv

##rarefy to see optimal sampling depth
qiime diversity alpha-rarefaction \
    --i-table HLCYG-taxa-filtered-table.qza \
    --p-max-depth 10000 \
    --m-metadata-file HL.Dilution.Metadata.V2.txt \
    --o-visualization HLCYG-taxa-filtered-rarefaction-plot-10k.qzv
##1836 is an appropriate depth for including all samples except negatives and nebnext which are
##getting excluded anyway. 

##diversity testing. first filter out nebnext, negatives, positives, which are checked separately
qiime feature-table filter-samples \
--i-table HLCYG-taxa-filtered-table.qza \
--m-metadata-file HL.Dilution.Metadata.V2.txt \
--p-where '[Sample_type_negative]="Urine"' \
--o-filtered-table HLCYG-urine.qza

qiime feature-table filter-samples \
--i-table HLCYG-urine.qza \
--m-metadata-file HL.Dilution.Metadata.V2.txt \
--p-where 'NOT [KitID]="N"' \
--o-filtered-table HLCYG-urine-no-neb-table.qza

qiime feature-table summarize \
--i-table HLCYG-urine-no-neb-table.qza \
--m-sample-metadata-file HL.Dilution.Metadata.V2.txt \
--o-visualization HLCYG-urine-no-neb-table.qzv

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table HLCYG-urine-no-neb-table.qza \
  --p-sampling-depth 1836 \
  --m-metadata-file HL.Dilution.Metadata.V2.txt \
  --output-dir HLCYG-1836-core-metrics


##beta diversity testing - weighted unifrac, unweighted unifrac, bray curtis, jaccard, by dog and kit
##kit
qiime diversity beta-group-significance \
  --i-distance-matrix HLCYG-1836-core-metrics/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file HL.Dilution.Metadata.V2.txt \
  --m-metadata-column KitID \
  --o-visualization  HLCYG-1836-core-metrics/unweighted-unifrac-kit-id-significance.qzv \
  --p-pairwise
qiime diversity beta-group-significance \
  --i-distance-matrix HLCYG-1836-core-metrics/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file HL.Dilution.Metadata.V2.txt \
  --m-metadata-column KitID \
  --o-visualization  HLCYG-1836-core-metrics/weighted-unifrac-kit-id-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix HLCYG-1836-core-metrics/bray_curtis_distance_matrix.qza \
  --m-metadata-file HL.Dilution.Metadata.V2.txt \
  --m-metadata-column KitID \
  --o-visualization  HLCYG-1836-core-metrics/bray-curtis-kit-id-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix HLCYG-1836-core-metrics/jaccard_distance_matrix.qza \
  --m-metadata-file HL.Dilution.Metadata.V2.txt \
  --m-metadata-column KitID \
  --o-visualization  HLCYG-1836-core-metrics/jaccard-kit-id-significance.qzv \
  --p-pairwise

##dog
qiime diversity beta-group-significance \
  --i-distance-matrix HLCYG-1836-core-metrics/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file HL.Dilution.Metadata.V2.txt \
  --m-metadata-column Dog \
  --o-visualization  HLCYG-1836-core-metrics/unweighted-unifrac-dog-significance.qzv \
  --p-pairwise
qiime diversity beta-group-significance \
  --i-distance-matrix HLCYG-1836-core-metrics/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file HL.Dilution.Metadata.V2.txt \
  --m-metadata-column Dog \
  --o-visualization  HLCYG-1836-core-metrics/weighted-unifrac-dog-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix HLCYG-1836-core-metrics/bray_curtis_distance_matrix.qza \
  --m-metadata-file HL.Dilution.Metadata.V2.txt \
  --m-metadata-column Dog \
  --o-visualization  HLCYG-1836-core-metrics/bray-curtis-dog-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix HLCYG-1836-core-metrics/jaccard_distance_matrix.qza \
  --m-metadata-file HL.Dilution.Metadata.V2.txt \
  --m-metadata-column Dog \
  --o-visualization  HLCYG-1836-core-metrics/jaccard-dog-significance.qzv \
  --p-pairwise

##these are for exporting raw counts to prism, not for the actual tests, 
##because qiime doesn't support friedman  
qiime diversity alpha-group-significance \
  --i-alpha-diversity HLCYG-1836-core-metrics/faith_pd_vector.qza \
  --m-metadata-file HL.Dilution.Metadata.V2.txt \
  --o-visualization HLCYG-1836-core-metrics/faith-pd-group-significance.qzv
qiime diversity alpha-group-significance \
  --i-alpha-diversity HLCYG-1836-core-metrics/observed_features_vector.qza \
  --m-metadata-file HL.Dilution.Metadata.V2.txt \
  --o-visualization HLCYG-1836-core-metrics/observed-features-group-significance.qzv
  qiime diversity alpha-group-significance \
  --i-alpha-diversity HLCYG-1836-core-metrics/shannon_vector.qza \
  --m-metadata-file HL.Dilution.Metadata.V2.txt \
  --o-visualization HLCYG-1836-core-metrics/shannon-group-significance.qzv


####code for generating barplots across kits
##filtering to at least 0.5% in 10% of samples for better plots
qiime feature-table filter-features-conditionally \
--i-table HLCYG-urine-no-neb-table.qza \
--p-abundance 0.005 \
--p-prevalence 0.1 \
--o-filtered-table abundance-filtered-table-hlcyg.qza

qiime taxa barplot \
--i-table abundance-filtered-table-hlcyg.qza \
--i-taxonomy HLCYG-taxonomy.qza \
--o-visualization hlcyg-abundance-filtered-barplots.qzv \
--m-metadata-file HL.Dilution.Metadata.V2.txt 

####summing by kit (i.e. all dogs summed together within each kit, to look at kit variability)
qiime feature-table group \
  --i-table abundance-filtered-table-hlcyg.qza \
  --m-metadata-file HL.Dilution.Metadata.V2.txt  \
  --m-metadata-column KitID \
  --p-axis sample \
  --p-mode sum \
  --o-grouped-table hlcyg-kit-grouped-table.qza

##make a barplot - don't pass a metadata file
qiime taxa barplot \
--i-table hlcyg-kit-grouped-table.qza \
--i-taxonomy HLCYG-taxonomy.qza \
--o-visualization hlcyg-kit-grouped-barplot.qzv 