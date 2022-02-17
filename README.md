# Plant genotype shapes the bacterial microbiome of fruits, leaves, and soil in olive plants 

### Antonino Malacrinò, Saveria Mosca, Maria Giulia Li Destri Nicosia, Giovanni E. Agosteo, Leonardo Schena


## Abstract 

The plant microbiome plays an important role in plant biology, ecology, and evolution. While recent technological developments enabled the characterization of plant-associated microbiota, we still know little about the impact of different biotic and abiotic factors on the diversity and structures of these microbial communities. Here, we characterized the structure of bacterial microbiomes of fruits, leaves, and soil collected from two olive genotypes (Sinopolese and Ottobratica), testing the hypothesis that plant genotype would impact each compartment with a different magnitude. Results show that plant genotype differently influenced the diversity, structure, composition, and co-occurence network at each compartment (fruits, leaves, soil), with a stronger effect on fruits compared to leaves and soil. Thus, plant genotype seems to be an important factor in shaping the structure of plant microbiomes in our system, and can be further explored to gain functional insights leading to improvements in plant productivity, nutrition, and defenses.


### Load packages

```R
library("dplyr")
library("phyloseq") 
library("vegan")
library("data.table") 
library("car")
library("lme4")
library("parallel")
library("tidyr")
library("stringr")
library("stringi")
library("ape")
library("ggplot2")
library("MuMIn")
library("emmeans")
library("ggpubr")
```

### Clean data

```R
ps.16S <- subset_taxa(ps.16S, Order !="Chloroplast")
ps.16S <- subset_taxa(ps.16S, Family !="Mitochondria")
```

### Split dataset by compartment

```R
ps.fruits <- phyloseq::subset_samples(ps.16S, Sample_type == "Fruits")
ps.leaves <- phyloseq::subset_samples(ps.16S, Sample_type == "Leaves")
ps.soil <- phyloseq::subset_samples(ps.16S, Sample_type == "Soil")
```

### PERMANOVA -- model

```R
sampledf <- data.frame(sample_data(ps.16S))
dist.mat <- phyloseq::distance(ps.16S, method = "unifrac")
perm <- how(nperm = 999)
set.seed(100)
pmv <- adonis2(dist.mat ~ Sample_type * Variety, data = sampledf, permutations = perm)
pmv
```

### PERMANOVA -- posthoc

```R
pairwise.perm.manova(dist.mat, paste0(sampledf$Sample_type, "_", sampledf$Variety), nperm = 999, progress = TRUE, p.method = "fdr", F = T, R2 = T)
```




























