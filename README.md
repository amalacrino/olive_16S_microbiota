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
library("NetCoMi")
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

### PERMANOVA -- fruits (model)

```R
sampledf <- data.frame(sample_data(ps.fruits))
dist.mat <- phyloseq::distance(ps.fruits, method = "unifrac")
perm <- how(nperm = 999)
set.seed(100)
pmv <- adonis2(dist.mat ~  Variety, data = sampledf, permutations = perm)
pmv
```

### PERMANOVA -- leaves (model)

```R
sampledf <- data.frame(sample_data(ps.leaves))
dist.mat <- phyloseq::distance(ps.leaves, method = "unifrac")
perm <- how(nperm = 999)
set.seed(100)
pmv <- adonis2(dist.mat ~  Variety, data = sampledf, permutations = perm)
pmv
```

### PERMANOVA -- soil (model)

```R
sampledf <- data.frame(sample_data(ps.soil))
dist.mat <- phyloseq::distance(ps.soil, method = "unifrac")
perm <- how(nperm = 999)
set.seed(100)
pmv <- adonis2(dist.mat ~  Variety, data = sampledf, permutations = perm)
pmv
```

### Calculate diversity indexes

```R
div <- microbiome::alpha(ps.16S, index = "all")
otus <- as.data.frame((otu_table(ps.16S)))
tree <- phy_tree(ps.16S)
div.pd <- pd(otus, tree, include.root = FALSE)
div.16s <- cbind(sample_data(ps.16S), div)
div.16s <- cbind(div.16s, div.pd)
```

### Phylogenetic diversity -- model

```R
model <- lm(PD ~ Sample_type * Variety, data = div.16s)
Anova(model)
```

### Phylogenetic diversity -- posthoc

```R
m1 <- emmeans(model, "Variety", by = "Sample_type")
pairs(m1)
```

### Differentially abundant ASVs -- function
```R
library("BiocParallel")
register(MulticoreParam(8))

diff.otus.fun <- function(x){
  diagdds = phyloseq_to_deseq2(x, ~ 1)
  ts = counts(diagdds)
  geoMeans = apply(ts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  diagdds = estimateSizeFactors(diagdds, geoMeans=geoMeans)
  diagdds = estimateDispersions(diagdds)
  diagdds$group <- factor(paste0(diagdds$Variety))
  design(diagdds) <- ~ group
  dds <-DESeq(diagdds, betaPrior=FALSE, parallel = T)
  tax.table <- as.data.frame(tax_table(x))
  c1 <- results(dds, parallel = T)
  c1 <- as.data.frame(c1)
  c1 <- setDT(c1, keep.rownames = TRUE)[]
  c1 <- c1[,c("rn", "log2FoldChange", "padj")]
  DiffOTUs<- cbind(tax.table, c1)
  DiffOTUs <- DiffOTUs[which(DiffOTUs$padj < 0.05),]
  return(DiffOTUs)
}
```

### Differentially abundant ASVs -- fruits
```R
diff.fruits <- diff.otus.fun(ps.fruits)
diff.fruits %>% filter(log2FoldChange > 0) %>% count(Genus, sort = TRUE)
diff.fruits %>% filter(log2FoldChange < 0) %>% count(Genus, sort = TRUE)
```

### Differentially abundant ASVs -- leaves
```R
diff.leaves <- diff.otus.fun(ps.leaves)
diff.leaves %>% filter(log2FoldChange > 0) %>% count(Genus, sort = TRUE)
diff.leaves %>% filter(log2FoldChange < 0) %>% count(Genus, sort = TRUE)
```

### Differentially abundant ASVs -- soil
```R
diff.soil <- diff.otus.fun(ps.soil)
diff.soil %>% filter(log2FoldChange > 0) %>% count(Genus, sort = TRUE)
diff.soil %>% filter(log2FoldChange < 0) %>% count(Genus, sort = TRUE)
```

### Impact of plant genotype on overall microbial community shift -- generate dataset

```R
diff.otus.fun2 <- function(x, compartment){
  diagdds = phyloseq_to_deseq2(x, ~ 1)
  ts = counts(diagdds)
  geoMeans = apply(ts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  diagdds = estimateSizeFactors(diagdds, geoMeans=geoMeans)
  diagdds = estimateDispersions(diagdds)
  diagdds$group <- factor(paste0(diagdds$Variety))
  design(diagdds) <- ~ group
  dds <-DESeq(diagdds, betaPrior=FALSE, parallel = T)
  tax.table <- as.data.frame(tax_table(x))
  c1 <- results(dds, parallel = T)
  c1 <- as.data.frame(c1)
  c1 <- setDT(c1, keep.rownames = TRUE)[]
  c1 <- c1[,c("rn", "log2FoldChange", "padj")]
  DiffOTUs<- cbind(tax.table, c1)
  DiffOTUs$log2FoldChange <- abs(DiffOTUs$log2FoldChange)
  DiffOTUs$compartment <- paste0(compartment)
  return(DiffOTUs)
}

a <- diff.otus.fun2(ps.fruits, "fruits")
b <- diff.otus.fun2(ps.leaves, "leaves")
c <- diff.otus.fun2(ps.soil, "soil")

dsq <- rbind(a, b, c)
```

### Impact of plant genotype on overall microbial community shift -- model

```R
model.dsq <- lmer(log2FoldChange ~ compartment * (1|rn), data = dsq)
Anova(model.dsq)
```

### Network analysis -- functions
```{r, warning=FALSE, message=FALSE, echo=FALSE}
ncores = 8
taxa <- as.data.frame(tax_table(ps.16S))
taxa <- setDT(taxa, keep.rownames = TRUE)[]
taxa <- taxa[,c(1,5:7)]

ps.fruits.sin <- phyloseq::subset_samples(ps.16S, Sample_type == "Fruits" & Variety == "Sinopolese")
ps.fruits.ott <- phyloseq::subset_samples(ps.16S, Sample_type == "Fruits" & Variety == "Ottobratica")
ps.leaves.sin <- phyloseq::subset_samples(ps.16S, Sample_type == "Leaves" & Variety == "Sinopolese")
ps.leaves.ott <- phyloseq::subset_samples(ps.16S, Sample_type == "Leaves" & Variety == "Ottobratica")
ps.soil.sin <- phyloseq::subset_samples(ps.16S, Sample_type == "Soil" & Variety == "Sinopolese")
ps.soil.ott <- phyloseq::subset_samples(ps.16S, Sample_type == "Soil" & Variety == "Ottobratica")

netTest <- function(x, g1, g2){
        comp_net <- netCompare(x, permTest = FALSE, verbose = FALSE, cores = ncores)
        a <- summary(comp_net, 
                groupNames = c(g1, g2),
                showCentr = c("degree", "between", "closeness"), 
                numbNodes = 5)
        return(a)
}

```

### Network analysis -- fruits
```{r, warning=FALSE, message=FALSE, echo=FALSE}
psObj1 <- ps.fruits.sin
psObj2 <- ps.fruits.ott
grNet1 <- "Sinopolese"
grNet2 <- "Ottobratica"
netTest(network, grNet1, grNet2)
```

### Network analysis -- leaves
```{r, warning=FALSE, message=FALSE, echo=FALSE}
psObj1 <- ps.leaves.sin
psObj2 <- ps.leaves.ott
grNet1 <- "Sinopolese"
grNet2 <- "Ottobratica"
netTest(network, grNet1, grNet2)
```

### Network analysis -- soil
```{r, warning=FALSE, message=FALSE, echo=FALSE}
psObj1 <- ps.soil.sin
psObj2 <- ps.soil.ott
grNet1 <- "Sinopolese"
grNet2 <- "Ottobratica"
netTest(network, grNet1, grNet2)
```
