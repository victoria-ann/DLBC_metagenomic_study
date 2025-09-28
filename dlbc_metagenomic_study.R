#This code block is used to determine kurtosis for the dataset
install.packages("moments")
library(moments)
shannon <- read.table("/Users/toripalecek/Documents/BMS690/alpha-diversity.tsv", header=TRUE, sep="\t")
kurtosis(shannon$shannon_entropy)
# kurtosis value of 2.17

#Create a box plot for Shannon Diversity
metadata <- read.delim("/Users/toripalecek/Documents/BMS690/blym_metadata.tsv")

colnames(shannon)[1] <- "SampleID"

shannon_df <- left_join(shannon, metadata, by = c(`SampleID` = "sampleid")) 

library(ggplot2)

ggplot(shannon_df, aes(x = group, y = shannon_entropy, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.9) +
  labs(title = "Shannon Diversity",
       x = NULL, y = "Shannon Diversity Index") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set1") +
  annotate("text", x = 1.0, y = max(shannon_df$shannon) + 0.1, label = " Kruskal-Wallis p = 0.25", size = 5)

#Create a box plot for Chao1 Diversity

chao1<- read.table("/Users/toripalecek/Documents/BMS690/chao1.tsv", header=TRUE, sep="\t")

colnames(chao1)[1] <- "SampleID"

chao1_df <- left_join(chao1, metadata, by = c(`SampleID` = "sampleid")) 

kruskal.test(chao1 ~ group, data = chao1_df)

ggplot(chao1_df, aes(x = group, y = chao1, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.9) +
  labs(title = "Chao1 Richness",
       x = NULL, y = "Chao1 Richness") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set1") +
  annotate("text", x = 1.0, y = max(chao1_df$chao1) + 0.1, label = " Kruskal-Wallis p = 0.14", size = 5)

#Create a box plot for Simpson Diversity
simpson<- read.table("/Users/toripalecek/Documents/BMS690/simpson.tsv", header=TRUE, sep="\t")

colnames(simpson)[1] <- "SampleID"

colnames(simpson)[2] <- "simpson"

simpson_df <- left_join(simpson, metadata, by = c(`SampleID` = "sampleid")) 

kruskal.test(simpson ~ group, data = simpson_df)

ggplot(simpson_df, aes(x = group, y = simpson, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.9) +
  labs(title = "Simpson Diversity",
       x = NULL, y = "Simpson Diversity Index") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set1") +
  annotate("text", x = 1.0, y = max(simpson$simpson) + 0.1, label = " Kruskal-Wallis p = 0.25", size = 5)

#Create a box plot for Observed Diversity

observed<- read.table("/Users/toripalecek/Documents/BMS690/observed.tsv", header=TRUE, sep="\t")

colnames(observed)[1] <- "SampleID"

colnames(observed)[2] <- "observed"

observed_df <- left_join(observed, metadata, by = c(`SampleID` = "sampleid")) 

kruskal.test(observed ~ group, data = observed_df)

ggplot(observed_df, aes(x = group, y = observed, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.9) +
  labs(title = "Observed Diversity",
       x = NULL, y = "Observed Diversity Index") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set1") +
  annotate("text", x = 1.0, y = max(observed$observed) + 0.1, label = " Kruskal-Wallis p = 0.14", size = 5)

#Create a box plot for Faith PD Diversity
faith<- read.table("/Users/toripalecek/Documents/BMS690/faith.tsv", header=TRUE, sep="\t")

colnames(faith)[1] <- "SampleID"

colnames(faith)[2] <- "faith"

faith_df <- left_join(faith, metadata, by = c(`SampleID` = "sampleid")) 

kruskal.test(faith ~ group, data = faith_df)

ggplot(faith_df, aes(x = group, y = faith, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.9) +
  labs(title = "Faith's Phylogenetic Diversity",
       x = NULL, y = "Faith's PD") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set1") +
  annotate("text", x = 1.0, y = max(faith$faith) + 0.1, label = " Kruskal-Wallis p = 0.25", size = 5)

#Create a PCoA Plot for Bray-Curtis
library(tidyverse)
install.packages("vegan")
library(vegan)
install.packages("ape")
library(ape)
library(ggplot2)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
library(phyloseq)

bray_curtis<- read.table("/Users/toripalecek/Documents/BMS690/bray_curtis.tsv", header=TRUE, sep="\t")

bray <- as.data.frame(bray_curtis)
rownames(bray) <- bray[[1]]
bray <- bray[,-1]
bray_dist <- as.dist(bray)

bray_pcoa <- ape::pcoa(bray_dist)

bray_df <- as.data.frame(bray_pcoa$vectors[,1:2])
bray_df$sampleid <- rownames(bray_df)
bray_df <- left_join(bray_df, metadata, by = "sampleid")

metadata$group <- as.factor(metadata$group)

permanova_bray <- adonis2(bray_dist ~ group, data = metadata)

bray_pval <- permanova_bray$`Pr(>F)`[1]

ggplot(bray_df, aes(x = Axis.1, y = Axis.2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCoA - Bray-Curtis",
    x = paste0("PCoA 1 (", round(bray_pcoa$values$Relative_eig[1]*100, 1), "%)"),
    y = paste0("PCoA 2 (", round(bray_pcoa$values$Relative_eig[2]*100, 1), "%)")
  ) + annotate("text", x = Inf, y = Inf, label = paste0("PERMANOVA p = ", signif(bray_pval, 3)),
               hjust = 1.1, vjust = 1.5, size = 3, fontface = "italic") + theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set1")


#Create a PCoA Plot for Unweighted Unifrac
unweighted_unifrac<- read.table("/Users/toripalecek/Documents/BMS690/unweighted_unifrac.tsv", header=TRUE, sep="\t")

unifrac <- as.data.frame(unweighted_unifrac)
rownames(unifrac) <- unifrac[[1]]
unifrac <- unifrac[,-1]
unifrac_dist <- as.dist(unifrac)

unifrac_pcoa <- ape::pcoa(unifrac_dist)

unifrac_df <- as.data.frame(unifrac_pcoa$vectors[,1:2])
unifrac_df$sampleid <- rownames(unifrac_df)
unifrac_df <- left_join(unifrac_df, metadata, by = "sampleid")

permanova_unifrac <- adonis2(unifrac_dist ~ group, data = metadata)

unifrac_pval <- permanova_unifrac$`Pr(>F)`[1]



ggplot(unifrac_df, aes(x = Axis.1, y = Axis.2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCoA - Unweighted UniFrac",
    x = paste0("PCoA 1 (", round(unifrac_pcoa$values$Relative_eig[1]*100, 1), "%)"),
    y = paste0("PCoA 2 (", round(unifrac_pcoa$values$Relative_eig[2]*100, 1), "%)")
  ) + annotate("text", x = Inf, y = Inf, label = paste0("PERMANOVA p = ", signif(unifrac_pval, 3)),
               hjust = 1.1, vjust = 1.5, size = 3, fontface = "italic") +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set1")


#This code block is used to determine how many of each taxonomic level there are
library(dplyr)
library(tidyr)
install.packages("tidyverse")
library(tidyverse)
install.packages("readr")
library(readr)

taxonomy <- read_tsv("/Users/toripalecek/Documents/BMS690/taxonomy.tsv")

taxonomy_clean <- taxonomy %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right") %>%
  mutate(across(Kingdom:Species, ~ gsub("^[a-z]__", "", .)))

# Count unique values at each taxonomic level
sapply(taxonomy_clean[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")], 
       function(x) length(unique(na.omit(x))))

# Import feature table 
feature_table_raw <- read.delim(
  "/Users/toripalecek/Documents/BMS690/feature-table.tsv",
  skip = 1,          # Skip the metadata line
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
head(feature_table_raw)
colnames(feature_table_raw)

# Set ASV IDs as row names
rownames(feature_table_raw) <- feature_table_raw$`#OTU ID`

# Remove the `#OTU ID` column
feature_table_raw$`#OTU ID` <- NULL

#Peform PLS-DA analysis
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("mixOmics")

library(mixOmics)

metadata <- read.delim("/Users/toripalecek/Documents/BMS690/blym_metadata.tsv")

feature_table_t <- t(feature_table_raw[, -1])

colnames(feature_table_t) <- feature_table_raw$`#OTU ID`

rownames(feature_table_t) <- colnames(feature_table_raw)[-1]

all(rownames(feature_table_t) %in% metadata$sampleid)

metadata_matched <- metadata[match(rownames(feature_table_t), metadata$sampleid), ]

group_factor <- as.factor(metadata_matched$group)

# Check for near-zero or zero variance columns
nzv <- apply(feature_table_t, 2, function(col) var(col, na.rm = TRUE) != 0)

# Filter the features
feature_table_filtered <- feature_table_t[, nzv]

# Now run PLS-DA
plsda_result <- plsda(X = feature_table_filtered, Y = group_factor, ncomp = 2)

perf_actual <- perf(plsda_result, validation = "Mfold", folds = 5, progressBar = FALSE)

observed_ber <- mean(perf_actual$error.rate$BER["comp1", ])

set.seed(123)  # for reproducibility
n_permutations <- 100
permuted_ber <- numeric(n_permutations)

# If you're starting from a phyloseq object:
X <- as.data.frame(t(otu_table(physeq)))  # transpose to make samples as rows
Y <- as.factor(sample_data(physeq)$group)  # or whatever grouping variable you’re using

# Make sure rownames match!
X <- X[rownames(sample_data(physeq)), ]

for (i in 1:n_permutations) {
  Y_perm <- sample(Y)}  # Shuffle the group labels


plsda_perm <- plsda(X, Y_perm, ncomp = 2)

perf_perm <- perf(plsda_perm, validation = "Mfold", folds = 5, progressBar = FALSE)

if (!is.null(perf_perm$error.rate$BER)) {
  ber_perm <- mean(perf_perm$error.rate$BER["comp1", ])  # Get BER for the first component
  permuted_ber[i] <- ber_perm
} else {
  permuted_ber[i] <- NA  # Handle if BER is not available
}

permuted_ber <- permuted_ber[!is.na(permuted_ber)]

# Compute observed BER (for the actual model)
observed_ber <- mean(perf_actual$error.rate$BER["comp1", ])

# Calculate p-value: how many permuted BERs were <= observed BER
p_value <- mean(permuted_ber <= observed_ber)

# Plot the permutation distribution
hist(permuted_ber, breaks = 20, main = "Permutation Test for PLS-DA",
     xlab = "Balanced Error Rate (BER)", col = "skyblue", border = "white")
abline(v = observed_ber, col = "red", lwd = 2)
legend("topright", legend = paste("Observed BER =", round(observed_ber, 3),
                                  "\nP-value =", round(p_value, 3)), bty = "n")

# Output the p-value
print(paste("Permutation test p-value:", round(p_value, 3)))

# Example variance explained — update with actual values from your model
comp1_var <- round(plsda_result$explained_variance$X[1] * 100, 1)
comp2_var <- round(plsda_result$explained_variance$X[2] * 100, 1)

plotIndiv(plsda_result,
          comp = c(1, 2),
          group = group_factor,
          legend = TRUE,
          ellipse = TRUE,
          star = TRUE,
          title = "PLS-DA Plot",
          xlab = paste0("PC1 (", comp1_var, "%)"),
          ylab = paste0("PC2 (", comp2_var, "%)"),
          cex = 2,  # Adjust font size of sample names (smaller for better fit)
          lwd = 0.02)   # Adjust line width for better visibility


#Create a phyloseq object
library(biomformat)
library(phyloseq)

feature_table <- read.table("/Users/toripalecek/Documents/BMS690/feature-table.tsv", sep = "\t", row.names = 1)
taxonomy <- read.table("/Users/toripalecek/Documents/BMS690/taxonomy.tsv", sep = "\t", header = TRUE, row.names = 1)
feature_table <- feature_table_raw

all.equal(rownames(feature_table), rownames(taxonomy))  # Should return TRUE now

otu_table <- otu_table(as.matrix(feature_table), taxa_are_rows = TRUE)

taxonomy <- tax_table(as.matrix(taxonomy))

sample_data <- sample_data(metadata)

sample_names(otu_table)

sample_names(sample_data)

head(metadata$sampleid)

rownames(metadata) <- sample_names(otu_table)

sample_data <- sample_data(metadata)

physeq <- phyloseq(otu_table, taxonomy, sample_data)

#Split Taxon in phyloseq object into its taxonomic rankings
tax_df <- as.data.frame(tax_table(physeq)) %>%
  rownames_to_column("FeatureID")

tax_split <- tax_df %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") %>%
  column_to_rownames("FeatureID")

tax_table(physeq) <- tax_table(as.matrix(tax_split))

rank_names(physeq)

#Make a relative abundance heat map
library(phyloseq)
library(ggplot2)
library(reshape2)

phyloseq_phylum <- tax_glom(physeq, "Phylum")

phyloseq_phylum_rel <- transform_sample_counts(phyloseq_phylum, function(x) x / sum(x))

phylo_df <- psmelt(phyloseq_phylum_rel)

phylo_df$group <- sample_data(phyloseq_phylum)[[ "group"]]

ggplot(phylo_df, aes(x = group, y = Phylum, fill = Abundance)) +
  geom_tile() +
  facet_grid(. ~ group, space = "free", scales = "free") + 
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "Phylum Relative Abundance by Sample Group", 
       x = "Sample Group", 
       y = "Phylum") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() + 
  theme(legend.position = "bottom") 


phyloseq_species <- tax_glom(physeq, "Species")

phyloseq_species_rel <- transform_sample_counts(phyloseq_species, function(x) x / sum(x))

species_df <- psmelt(phyloseq_species_rel)

species_df$group <- sample_data(phyloseq_species)[[ "group"]]

ggplot(species_df, aes(x = group, y = Species, fill = Abundance)) +
  geom_tile() +
  facet_grid(. ~ group, space = "free", scales = "free") + 
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "Species Relative Abundance by Sample Group", 
       x = "Sample Group", 
       y = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() + 
  theme(legend.position = "bottom") 


phyloseq_genus <- tax_glom(physeq, "Genus")

phyloseq_genus_rel <- transform_sample_counts(phyloseq_genus, function(x) x / sum(x))

genus_df <- psmelt(phyloseq_genus_rel)

genus_df$group <- sample_data(phyloseq_genus)[[ "group"]]

ggplot(genus_df, aes(x = group, y = Genus, fill = Abundance)) +
  geom_tile() +
  facet_grid(. ~ group, space = "free", scales = "free") + 
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "Genus Relative Abundance by Sample Group", 
       x = "Sample Group", 
       y = "Genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() + 
  theme(legend.position = "bottom") 


phyloseq_family <- tax_glom(physeq, "Family")

phyloseq_family_rel <- transform_sample_counts(phyloseq_family, function(x) x / sum(x))

family_df <- psmelt(phyloseq_family_rel)

family_df$group <- sample_data(phyloseq_family)[[ "group"]]

ggplot(family_df, aes(x = group, y = Family, fill = Abundance)) +
  geom_tile() +
  facet_grid(. ~ group, space = "free", scales = "free") + 
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "Family Relative Abundance by Sample Group", 
       x = "Sample Group", 
       y = "Family") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() + 
  theme(legend.position = "bottom") 


#Targeted Relative Abundance Box Plots

library(phyloseq)
library(dplyr)
library(ggplot2)

# Select the taxonomic rank you want (e.g., "Genus")
physeq_genus <- tax_glom(physeq, taxrank = "Genus")

# Extract the OTU table and taxonomy data
otu_table_genus <- otu_table(physeq_genus)
taxonomy_genus <- tax_table(physeq_genus)

# Convert to data frame for easy manipulation
otu_df <- as.data.frame(otu_table_genus)
taxonomy_df <- as.data.frame(taxonomy_genus)

# Add taxonomic information to the OTU data frame
otu_df$Taxa <- rownames(taxonomy_df)
otu_df <- merge(otu_df, taxonomy_df, by.x = "Taxa", by.y = "row.names")

# Check what taxonomic values are available in the 'Genus' column
unique(otu_df$Genus)

# Remove leading and trailing spaces in the Genus column
otu_df$Genus <- trimws(otu_df$Genus)

# Verify the cleaned Genus values
unique(otu_df$Genus)

# Filter the data for specific genera
otu_df_filtered <- otu_df %>%
  filter(!is.na(Genus)) %>%
  filter(Genus %in% c("g__Allisonella", "g__Enterococcus", "g__Escherichia-Shigella", 
                      "g__Lachnospira", "g__Prevotella_2", "g__Roseburia", "g__Veillonella"))

# Ensure that the 'Genus' column is treated as character
otu_df_filtered <- otu_df_filtered %>%
  mutate(Genus = as.character(Genus)) %>%
  # Convert only relative abundance columns to numeric, excluding the categorical ones
  mutate(across(!Genus & !contains(c("Kingdom", "Phylum", "Class", "Order", "Family", "Species", "Confidence")), as.numeric))

# Now pivot_longer
# Ensure that only numeric columns are pivoted, excluding taxonomic information
otu_df_long <- otu_df_filtered %>%
  pivot_longer(cols = where(is.numeric),  # Pivot only numeric columns
               names_to = "SampleID",  # The column with sample names
               values_to = "RelativeAbundance")


# Add metadata (assuming 'physeq' has a 'sample_data' dataframe with 'Group' column)
# Rename the 'SampleID' column in otu_df_long to match 'sampleid' in metadata
otu_df_long <- otu_df_long %>%
  rename(sampleid = SampleID)  # Renaming 'SampleID' to 'sampleid'

# Add metadata (assuming 'physeq' has a 'sample_data' dataframe with 'sampleid' column)
metadata <- sample_data(physeq_genus)  # Replace with the actual name if different
otu_df_long <- otu_df_long %>%
  left_join(metadata, by = "sampleid")  # Join by 'sampleid'

# Create box plots by test group and Genus
ggplot(otu_df_long, aes(x = Genus, y = RelativeAbundance, fill = group)) +  # Replace 'Group' with your test group column
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Boxplot of Targeted Genus Relative Abundance by Test Group",
       x = "Genus",
       y = "Relative Abundance") +
  scale_fill_manual(values = c("red", "blue"))  # Adjust color palette if necessary



# Select the taxonomic rank you want (e.g., "Species")
physeq_species <- tax_glom(physeq, taxrank = "Species")

# Extract the OTU table and taxonomy data
otu_table_species <- otu_table(physeq_species)
taxonomy_species <- tax_table(physeq_species)

# Convert to data frame for easy manipulation
otu_df <- as.data.frame(otu_table_species)
taxonomy_df <- as.data.frame(taxonomy_species)

# Add taxonomic information to the OTU data frame
otu_df$Taxa <- rownames(taxonomy_df)
otu_df <- merge(otu_df, taxonomy_df, by.x = "Taxa", by.y = "row.names")

# Check what taxonomic values are available in the 'Species' column
unique(otu_df$Species)

# Remove leading and trailing spaces in the Species column
otu_df$Species <- trimws(otu_df$Species)

# Verify the cleaned Species values
unique(otu_df$Species)

# Filter the data for specific genera
otu_df_filtered <- otu_df %>%
  filter(!is.na(Species)) %>%
  filter(Species %in% c("s__Bacteroides_fragilis", "s__Clostridium_butyricum", "s__Escherichia_coli","s__Lactococcus_garvieae"))

# Ensure that the 'Species' column is treated as character
otu_df_filtered <- otu_df_filtered %>%
  mutate(Species = as.character(Species)) %>%
  # Convert only relative abundance columns to numeric, excluding the categorical ones
  mutate(across(!Species & !contains(c("Kingdom", "Phylum", "Class", "Order", "Family", "Species", "Confidence")), as.numeric))

# Now pivot_longer
# Ensure that only numeric columns are pivoted, excluding taxonomic information
otu_df_long <- otu_df_filtered %>%
  pivot_longer(cols = where(is.numeric),  # Pivot only numeric columns
               names_to = "SampleID",  # The column with sample names
               values_to = "RelativeAbundance")


# Add metadata (assuming 'physeq' has a 'sample_data' dataframe with 'Group' column)
# Rename the 'SampleID' column in otu_df_long to match 'sampleid' in metadata
otu_df_long <- otu_df_long %>%
  rename(sampleid = SampleID)  # Renaming 'SampleID' to 'sampleid'

# Add metadata (assuming 'physeq' has a 'sample_data' dataframe with 'sampleid' column)
metadata <- sample_data(physeq_species)  # Replace with the actual name if different
otu_df_long <- otu_df_long %>%
  left_join(metadata, by = "sampleid")  # Join by 'sampleid'

# Create box plots by test group and Species
ggplot(otu_df_long, aes(x = Species, y = RelativeAbundance, fill = group)) +  # Replace 'Group' with your test group column
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Boxplot of Targeted Species Relative Abundance by Test Group",
       x = "Species",
       y = "Relative Abundance") +
  scale_fill_manual(values = c("red", "blue"))  # Adjust color palette if necessary





# Select the taxonomic rank you want (e.g., "Phylum")
physeq_phylum <- tax_glom(physeq, taxrank = "Phylum")

# Extract the OTU table and taxonomy data
otu_table_phylum <- otu_table(physeq_phylum)
taxonomy_phylum <- tax_table(physeq_phylum)

# Convert to data frame for easy manipulation
otu_df <- as.data.frame(otu_table_phylum)
taxonomy_df <- as.data.frame(taxonomy_phylum)

# Add taxonomic information to the OTU data frame
otu_df$Taxa <- rownames(taxonomy_df)
otu_df <- merge(otu_df, taxonomy_df, by.x = "Taxa", by.y = "row.names")

# Check what taxonomic values are available in the 'Phylum' column
unique(otu_df$Phylum)

# Remove leading and trailing spaces in thePhylums column
otu_df$Phylum <- trimws(otu_df$Phylum)

# Verify the cleaned Phylum values
unique(otu_df$Phylum)

# Filter the data for specific genera
otu_df_filtered <- otu_df %>%
  filter(!is.na(Phylum)) %>%
  filter(Phylum %in% c("p__Proteobacteria"))

# Ensure that the 'Phylum' column is treated as character
otu_df_filtered <- otu_df_filtered %>%
  mutate(Species = as.character(Phylum)) %>%
  # Convert only relative abundance columns to numeric, excluding the categorical ones
  mutate(across(!Phylum & !contains(c("Kingdom", "Phylum", "Class", "Order", "Family", "Species", "Confidence")), as.numeric))

# Now pivot_longer
# Ensure that only numeric columns are pivoted, excluding taxonomic information
otu_df_long <- otu_df_filtered %>%
  pivot_longer(cols = where(is.numeric),  # Pivot only numeric columns
               names_to = "SampleID",  # The column with sample names
               values_to = "RelativeAbundance")


# Add metadata (assuming 'physeq' has a 'sample_data' dataframe with 'Group' column)
# Rename the 'SampleID' column in otu_df_long to match 'sampleid' in metadata
otu_df_long <- otu_df_long %>%
  rename(sampleid = SampleID)  # Renaming 'SampleID' to 'sampleid'

# Add metadata (assuming 'physeq' has a 'sample_data' dataframe with 'sampleid' column)
metadata <- sample_data(physeq_phylum)  # Replace with the actual name if different
otu_df_long <- otu_df_long %>%
  left_join(metadata, by = "sampleid")  # Join by 'sampleid'

# Create box plots by test group and Phylum
ggplot(otu_df_long, aes(x = Phylum, y = RelativeAbundance, fill = group)) +  # Replace 'Group' with your test group column
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Boxplot of Targeted Phylum Relative Abundance by Test Group",
       x = "Phylum",
       y = "Relative Abundance") +
  scale_fill_manual(values = c("red", "blue"))  # Adjust color palette if necessary



# Select the taxonomic rank you want (e.g., "Family")
physeq_family <- tax_glom(physeq, taxrank = "Family")

# Extract the OTU table and taxonomy data
otu_table_family <- otu_table(physeq_family)
taxonomy_family <- tax_table(physeq_family)

# Convert to data frame for easy manipulation
otu_df <- as.data.frame(otu_table_family)
taxonomy_df <- as.data.frame(taxonomy_family)

# Add taxonomic information to the OTU data frame
otu_df$Taxa <- rownames(taxonomy_df)
otu_df <- merge(otu_df, taxonomy_df, by.x = "Taxa", by.y = "row.names")

# Check what taxonomic values are available in the 'Phylum' column
unique(otu_df$Family)

# Remove leading and trailing spaces in thePhylums column
otu_df$Family <- trimws(otu_df$Family)

# Verify the cleaned Family values
unique(otu_df$Family)

# Filter the data for specific genera
otu_df_filtered <- otu_df %>%
  filter(!is.na(Family)) %>%
  filter(Family %in% c("f__Enterobacteriaceae"))

# Ensure that the 'Family' column is treated as character
otu_df_filtered <- otu_df_filtered %>%
  mutate(Family = as.character(Family)) %>%
  # Convert only relative abundance columns to numeric, excluding the categorical ones
  mutate(across(!Family & !contains(c("Kingdom", "Phylum", "Class", "Order", "Family", "Species", "Confidence")), as.numeric))

# Now pivot_longer
# Ensure that only numeric columns are pivoted, excluding taxonomic information
otu_df_long <- otu_df_filtered %>%
  pivot_longer(cols = where(is.numeric),  # Pivot only numeric columns
               names_to = "SampleID",  # The column with sample names
               values_to = "RelativeAbundance")


# Add metadata (assuming 'physeq' has a 'sample_data' dataframe with 'Group' column)
# Rename the 'SampleID' column in otu_df_long to match 'sampleid' in metadata
otu_df_long <- otu_df_long %>%
  rename(sampleid = SampleID)  # Renaming 'SampleID' to 'sampleid'

# Add metadata (assuming 'physeq' has a 'sample_data' dataframe with 'sampleid' column)
metadata <- sample_data(physeq_family)  # Replace with the actual name if different
otu_df_long <- otu_df_long %>%
  left_join(metadata, by = "sampleid")  # Join by 'sampleid'

# Create box plots by test group and Family
ggplot(otu_df_long, aes(x = Family, y = RelativeAbundance, fill = group)) +  # Replace 'Group' with your test group column
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Boxplot of Targeted Family Relative Abundance by Test Group",
       x = "Family",
       y = "Relative Abundance") +
  scale_fill_manual(values = c("red", "blue"))  # Adjust color palette if necessary


# Select the taxonomic rank you want (e.g., "Order")
physeq_order <- tax_glom(physeq, taxrank = "Order")

# Extract the OTU table and taxonomy data
otu_table_order <- otu_table(physeq_order)
taxonomy_order <- tax_table(physeq_order)

# Convert to data frame for easy manipulation
otu_df <- as.data.frame(otu_table_order)
taxonomy_df <- as.data.frame(taxonomy_order)

# Add taxonomic information to the OTU data frame
otu_df$Taxa <- rownames(taxonomy_df)
otu_df <- merge(otu_df, taxonomy_df, by.x = "Taxa", by.y = "row.names")

# Check what taxonomic values are available in the 'Phylum' column
unique(otu_df$Order)

# Remove leading and trailing spaces in thePhylums column
otu_df$Order <- trimws(otu_df$Order)

# Verify the cleaned Order values
unique(otu_df$Order)

# Filter the data for specific genera
otu_df_filtered <- otu_df %>%
  filter(!is.na(Order)) %>%
  filter(Order %in% c("o__Enterobacterales"))

# Ensure that the 'Order' column is treated as character
otu_df_filtered <- otu_df_filtered %>%
  mutate(Family = as.character(Order)) %>%
  # Convert only relative abundance columns to numeric, excluding the categorical ones
  mutate(across(!Order & !contains(c("Kingdom", "Phylum", "Class", "Order", "Family", "Species", "Confidence")), as.numeric))

# Now pivot_longer
# Ensure that only numeric columns are pivoted, excluding taxonomic information
otu_df_long <- otu_df_filtered %>%
  pivot_longer(cols = where(is.numeric),  # Pivot only numeric columns
               names_to = "SampleID",  # The column with sample names
               values_to = "RelativeAbundance")


# Add metadata (assuming 'physeq' has a 'sample_data' dataframe with 'Group' column)
# Rename the 'SampleID' column in otu_df_long to match 'sampleid' in metadata
otu_df_long <- otu_df_long %>%
  rename(sampleid = SampleID)  # Renaming 'SampleID' to 'sampleid'

# Add metadata (assuming 'physeq' has a 'sample_data' dataframe with 'sampleid' column)
metadata <- sample_data(physeq_order)  # Replace with the actual name if different
otu_df_long <- otu_df_long %>%
  left_join(metadata, by = "sampleid")  # Join by 'sampleid'

# Create box plots by test group and Order
ggplot(otu_df_long, aes(x = Order, y = RelativeAbundance, fill = group)) +  # Replace 'Group' with your test group column
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Boxplot of Targeted Order Relative Abundance by Test Group",
       x = "Order",
       y = "Relative Abundance") +
  scale_fill_manual(values = c("red", "blue"))  # Adjust color palette if necessary



# Select the taxonomic rank you want (e.g., "Class")
physeq_class <- tax_glom(physeq, taxrank = "Class")

# Extract the OTU table and taxonomy data
otu_table_class <- otu_table(physeq_class)
taxonomy_class <- tax_table(physeq_class)

# Convert to data frame for easy manipulation
otu_df <- as.data.frame(otu_table_class)
taxonomy_df <- as.data.frame(taxonomy_class)

# Add taxonomic information to the OTU data frame
otu_df$Taxa <- rownames(taxonomy_df)
otu_df <- merge(otu_df, taxonomy_df, by.x = "Taxa", by.y = "row.names")

# Check what taxonomic values are available in the 'Class' column
unique(otu_df$Class)

# Remove leading and trailing spaces in the Class column
otu_df$Class <- trimws(otu_df$Class)

# Verify the cleanedClassr values
unique(otu_df$Class)

# Filter the data for specific genera
otu_df_filtered <- otu_df %>%
  filter(!is.na(Class)) %>%
  filter(Class %in% c("c__Gammaproteobacteria"))

# Ensure that the 'Class' column is treated as character
otu_df_filtered <- otu_df_filtered %>%
  mutate(Family = as.character(Class)) %>%
  # Convert only relative abundance columns to numeric, excluding the categorical ones
  mutate(across(!Class & !contains(c("Kingdom", "Phylum", "Class", "Order", "Family", "Species", "Confidence")), as.numeric))

# Now pivot_longer
# Ensure that only numeric columns are pivoted, excluding taxonomic information
otu_df_long <- otu_df_filtered %>%
  pivot_longer(cols = where(is.numeric),  # Pivot only numeric columns
               names_to = "SampleID",  # The column with sample names
               values_to = "RelativeAbundance")


# Add metadata (assuming 'physeq' has a 'sample_data' dataframe with 'Group' column)
# Rename the 'SampleID' column in otu_df_long to match 'sampleid' in metadata
otu_df_long <- otu_df_long %>%
  rename(sampleid = SampleID)  # Renaming 'SampleID' to 'sampleid'

# Add metadata (assuming 'physeq' has a 'sample_data' dataframe with 'sampleid' column)
metadata <- sample_data(physeq_class)  # Replace with the actual name if different
otu_df_long <- otu_df_long %>%
  left_join(metadata, by = "sampleid")  # Join by 'sampleid'

# Create box plots by test group and Class
ggplot(otu_df_long, aes(x = Class, y = RelativeAbundance, fill = group)) +  # Replace 'Group' with your test group column
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Boxplot of Targeted Class Relative Abundance by Test Group",
       x = "Class",
       y = "Relative Abundance") +
  scale_fill_manual(values = c("red", "blue"))  # Adjust color palette if necessary

# Manual kruskal-wallis test
kruskal_results <- genus_df %>%
  group_by(Genus) %>%
  summarise(p = kruskal.test(Abundance ~ group)$p.value) %>%
  mutate(p_adj = p.adjust(p, method = "fdr"))


# ANCOM Anaysis wtih Volcano Plot
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")

library(ANCOMBC)
library(phyloseq)

?ancombc

library(ANCOMBC)
install.packages("microbiome")


library(phyloseq)
library(ANCOMBC)

# Extract the count table (OTU/ASV table)
counts <- as.data.frame(otu_table(physeq_species))
if (taxa_are_rows(physeq_species)) {
  counts <- counts
} else {
  counts <- t(counts)
}

# Extract the sample metadata
meta <- as.data.frame(sample_data(physeq_species))

# Extract taxonomy (optional)
taxa <- as.data.frame(tax_table(physeq_species))

str(meta)

meta_df <- as.data.frame(meta)
class(meta_df)
is.data.frame(meta_df)

meta_df$group <- as.factor(meta_df$group)

meta_df <- as.data.frame(as(sample_data(physeq), "data.frame"))
rownames(meta_df) <- meta_df$sampleid

out <- ancombc2(
  data = counts,
  meta_data = meta_df,
  group = "group",
  fix_formula = "group",
  tax_level = NULL,
  prv_cut = 0.10,
  lib_cut = 1000,
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  global = FALSE,
  n_cl = 1,
  verbose = TRUE
)

names(out)

res <- out$res
names(res)

sig_taxa <- rownames(res$diff_abn)[res$diff_abn$group == TRUE]
sig_taxa

sig_results <- data.frame(
  Taxon = rownames(res$lfc),
  LFC = res$lfc$group,
  Qval = res$q_val$group,
  Significant = res$diff_abn$group
)

sig_results <- sig_results[sig_results$Significant == TRUE, ]
head(sig_results)

library(ggplot2)

volcano_df <- data.frame(
  Taxon = rownames(res$lfc),
  LFC = res$lfc$group,
  Qval = res$q_val$group,
  Significant = res$diff_abn$group
)

ggplot(volcano_df, aes(x = LFC, y = -log10(Qval), color = Significant)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "ANCOM-BC2 Volcano Plot",
       x = "Log Fold Change",
       y = "-log10(FDR-adjusted p-value)") +
  theme_minimal()


