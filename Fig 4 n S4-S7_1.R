## 1) Needed packages

# Graph tools
library("ggplot2")
library("RColorBrewer")
library("ggpubr")
library("dplyr")
library("patchwork") # ggplot modifiers
library(ggprism) # Prism-inspired graphics theme for ggplot
library(cowplot) # more plot modifier tools
library(ggh4x) # facet plotter

# Microbiome-specific tools
library("ShortRead") # Manage and explore DNA libraries (step 3)
library("dada2") # Infere sample from sequences (steps 4-6)
library("microbiome") #
library("phyloseq") #
library("DECIPHER") # Alternative algorithm to assign taxonomy to ASV table (step 6) and necessary for philogenetic tree (step 8) **** Check if this library does not conflict with steps 1-5
library("phangorn") # Philogenetic tree algorithms
library("ape") # Necessary to use some phagorn functions
library("devtools")
library("ranacapa") # Necessary for rarefraction of data
# library("d3heatmap") # Plot interactive D3 heatmaps
library("htmlwidgets") # Export D3 plots to standalone HTML file
library("FactoMineR") # PCA fuction
library("factoextra") # PCA visualization tools
library("ggbiplot") # Biploting
library("ade4") # 
library("vegan") # PERMANOVA
library("DESeq2") #
library("reshape2") #
library(GUniFrac) # ZicoSeq
library("ComplexUpset") # Upset plot
library("MicrobiotaProcess") # Microbiome analysis toolbox. Used just for the get_upset function
library(microshades) # Plot abundance with microshades of colors

## 2) Importing the phyloseq file

# setwd("/Users/davidcajasmunoz/Library/CloudStorage/GoogleDrive-dadavid.cajas@gmail.com/Mi unidad/Academia/Postgrado/PUCV/Tesis/Resultados/Biofertilización/2/5-ADN/R/Seq")
setwd("G:/My Drive/Academia/Postgrado/PUCV/Tesis/Resultados/Biofertilización/2/5-ADN/R/Seq")

ps<-readRDS("phyloseq132.rds")
ps
print("ps file imported successfuly")

# We can inspect the components of the phyloseq file
# View(otu_table(ps)@.Data)
# meta(ps)
# View(tax_table(ps)@.Data)
# ntaxa(ps)
# nsamples(ps)
# psmelt(ps)

# We can add or modify specific parts of the phyloseq file
# merge_phyloseq()

## 3) Cleaning and preparing phyloseq file

# 3.1) Eliminate non-bacterial taxa
rank_names(ps)
get_taxa_unique(ps, "Kingdom")
ps_bac <- subset_taxa(ps, Kingdom == "Bacteria")
get_taxa_unique(ps_bac, "Kingdom")
print("Non-bacterial taxa deleted")

# 3.2) Eliminate extremely low abundance (< 100 counts) ASVs

ps_bac2 <- prune_taxa(taxa_sums(ps_bac) > 100, ps_bac)

ntaxa(ps_bac)
ntaxa(ps_bac2)

# 3.2) Add read count to metadata
readcounts <- readcount(ps_bac2)
sample_data(ps_bac2)$Readcount_no_rarefaction <- readcounts
#View(meta(ps_bac2))
print("Metadata now shows Read count")

# 3.3) Add sample names to metadata
sample_data(ps_bac2)$Sample_names<-as.factor(sample_names(ps_bac2))
print("Metadata now shows sample names")

# 3.4) Change column types of variables to factor

sample_variables(ps_bac2)
str(meta(ps_bac2))

sample_data(ps_bac2)$Sample_names<-as.factor(sample_names(ps_bac2))
sample_data(ps_bac2)$Treatment<-as.factor(sample_data(ps_bac2)$Treatment)
sample_data(ps_bac2)$Replica<-as.factor(sample_data(ps_bac2)$Replica)
sample_data(ps_bac2)$Fertilization<-as.factor(sample_data(ps_bac2)$Fertilization)
sample_data(ps_bac2)$Fertilization <- factor(sample_data(ps_bac2)$Fertilization, levels = c("No","Hoagland","AMH3-8","No_plant")) #Reordering levels in Fertilization variable

str(meta(ps_bac2))
print("Metadata type fixed")

# 3.5) Study the effect of sample Readcount over community

# 3.5.1) Rarefaction curve analysis

t <- proc.time() # Timer
rarefacurve <- ggrare(ps_bac2, step = 500, color = "Copper", label = "Sample_names", se = TRUE)
proc.time()-t
rarefacurve+ggtitle("Rarefaction Curve")+xlim(c(0,40000))+
  theme(plot.title=element_text(size=13, 
                                family="American Typewriter",
                                color="black",
                                hjust=0.5,
                                lineheight=1.2))+
  scale_x_continuous(breaks=seq(0, 40000, 5000))
print(rarefacurve)

# 3.5.2) Scatterplot diversity index vs reads

adiv_original <- microbiome::alpha(ps_bac2, index = "all")

adiv_original$ReadsPerSample <- sample_sums(ps_bac2)
p1 <- ggscatter(adiv_original, "ReadsPerSample","observed", 
                add = "loess") + stat_cor(method = "pearson") + ggtitle(label = "Observed diversity")

p2 <- ggscatter(adiv_original,"ReadsPerSample", "diversity_shannon", 
                add = "loess") + stat_cor(method = "pearson") + ggtitle(label = "Shannon index")

p3 <- ggscatter(adiv_original, "ReadsPerSample","diversity_fisher", 
                add = "loess") + stat_cor(method = "pearson") + ggtitle(label = "Fisher index")
ggarrange(p1,p2,p3, ncol=3, nrow = 1)

# 3.6) Rarefaction of data

# To ensure reproducibility, it is best to set seed for random processes
set.seed(1)

# Summary of readcounts per sample
summary(sample_sums(ps_bac2))

ps_rar <- rarefy_even_depth(ps_bac2, sample.size = 20000)
# ps_rar_2k <- rarefy_even_depth(ps_bac, sample.size = 2000)

# Save cleaned phyloseq object on file
saveRDS(ps_rar, "phyloseq_clean.rds")
#saveRDS(ps_rar_2k, "phyloseq_clean2k.rds")
ps_rar <- readRDS("phyloseq_clean.rds")
#ps_rar2 <- readRDS("phyloseq_clean2k.rds")
print(paste("Data rarefacted at",readcount(ps_rar)[1],"reads per sample",sep = " "))

# 3.6.2) Update read count to metadata
readcounts <- readcount(ps_rar)
sample_data(ps_rar)$Readcount_rarefacted <- readcounts
# View(meta(ps_rar))

# 3.7) Sort subsets of the phyloseq file

sample_variables(ps_rar)

# Subset rhizospheric samples

# ps_rar_r <- rarefy_even_depth(ps_r, sample.size = 20000)
ps_rar_r <- subset_samples(ps_rar, Fertilization != "No_plant")

# Subset by Soil copper content

ps_nocu <- subset_samples(ps_rar, Copper == "FALSE")
ps_nocu <- prune_taxa(taxa_sums(ps_nocu) > 100, ps_nocu)

ps_yescu <- subset_samples(ps_rar, Copper == "TRUE")
ps_yescu <- prune_taxa(taxa_sums(ps_yescu) > 100, ps_yescu)
ps_yescu

# Aggregate replicas for abundance plot

melt_nocu <- ps_nocu %>%
  merge_samples("Fertilization") %>%
  # tax_glom("Genus") %>% # NOTE: If activate this method; aggregates ASVs by genus, thus losing all the "unknown" genera, even if other ranks are known.
  phyloseq::transform_sample_counts(function(x) { x/sum(x) }) %>%
  psmelt() %>%
  filter(Abundance > 0)
melt_nocu[melt_nocu["Sample"] == "No", "Sample"] = "Control"
# melt_nocu[melt_nocu["Sample"] == "No_plant", "Sample"] = "Bulk soil" # EN
melt_nocu[melt_nocu["Sample"] == "No_plant", "Sample"] = "Suelo NR" #ES
melt_nocu$Sample <- factor(melt_nocu$Sample, levels = c('Control',
                                                        'Hoagland',
                                                        'AMH3-8',
                                                        'Suelo NR'))

melt_yescu <- ps_yescu %>%
  merge_samples("Fertilization") %>%
  # tax_glom("Genus") %>% # NOTE: If activate this method; aggregates ASVs by genus, thus losing all the "unknown" genera, even if other ranks are known.
  phyloseq::transform_sample_counts(function(x) { x/sum(x) }) %>%
  psmelt() %>%
  filter(Abundance > 0)
melt_yescu[melt_yescu["Sample"] == "No", "Sample"] = "Control"
# melt_yescu[melt_yescu["Sample"] == "No_plant", "Sample"] = "Bulk soil" # EN
melt_yescu[melt_yescu["Sample"] == "No_plant", "Sample"] = "Suelo NR" #ES
melt_yescu$Sample <- factor(melt_yescu$Sample, levels = c('Control',
                                                        'Hoagland',
                                                        'AMH3-8',
                                                        'Suelo NR'))

#Bind back together No copper and Yes Copper melted dataframes
melt_full <- rbind(melt_nocu,melt_yescu)
melt_full[melt_full["Copper"] == "0", "Copper"] = "Sin Cobre"
melt_full[melt_full["Copper"] == "1", "Copper"] = "Con Cobre"
melt_full$Copper <- factor(melt_full$Copper, levels = c('Sin Cobre',
                                                          'Con Cobre'))

## 3.8) Subsets of microbiome sample

# 3.8.1) Check abundance distribution of taxa

# 3.8.1.1) By Phylum

phyl <- aggregate_taxa(ps_rar, level = "Phylum")
ntaxa(phyl)
get_taxa_unique(phyl,"Phylum")
tax_ab_phyl <- sort(taxa_sums(phyl), decreasing = TRUE)
tax_ab_phyl

phyl_nocu <- aggregate_taxa(ps_nocu, level = "Phylum")
ntaxa(phyl_nocu)
get_taxa_unique(phyl_nocu,"Phylum")
tax_nocu_phyl <- sort(taxa_sums(phyl_nocu), decreasing = TRUE)
tax_nocu_phyl

phyl_yescu <- aggregate_taxa(ps_yescu, level = "Phylum")
ntaxa(phyl_yescu)
get_taxa_unique(phyl_yescu,"Phylum")
tax_yescu_phyl <- sort(taxa_sums(phyl_yescu), decreasing = TRUE)
tax_yescu_phyl

# 3.8.1.2) By genus

gen <- aggregate_taxa(ps_rar, level = "Genus")
ntaxa(gen)
tax_ab_gen <- sort(taxa_sums(gen), decreasing = TRUE)
tax_ab_gen

gen_nocu <- aggregate_taxa(ps_nocu, level = "Genus")
ntaxa(gen_nocu)
tax_nocu_gen <- sort(taxa_sums(gen_nocu), decreasing = TRUE)
tax_nocu_gen

gen_yescu <- aggregate_taxa(ps_yescu, level = "Genus")
ntaxa(gen_yescu)
tax_yescu_gen <- sort(taxa_sums(gen_yescu), decreasing = TRUE)
tax_yescu_gen

# 3.8.2.1) Generate rare taxa (< 1% abundance) phyloseq objects

ps_rar_rare <- subset_taxa(ps_rar, abundances(ps_rar, "compositional") < 0.01)
ps_rar_dom <- subset_taxa(ps_rar, abundances(ps_rar, "compositional") > 0.01)

ps_nocu_rare <- subset_taxa(ps_nocu, abundances(ps_nocu, "compositional") < 0.01)
ps_nocu_dom <- subset_taxa(ps_nocu, abundances(ps_nocu, "compositional") > 0.01)

ps_yescu_rare <- subset_taxa(ps_yescu, abundances(ps_yescu, "compositional") < 0.01)
ps_yescu_dom <- subset_taxa(ps_yescu, abundances(ps_yescu, "compositional") > 0.01)

# 3.8.2.1) Generate new phyloseq file excluding 5 most abundant phyla

ps_rarep_nocu <- ps_nocu
for (taxarank in 1:5) {
  ps_rarep_nocu <- subset_taxa(ps_rarep_nocu, Phylum != top_taxa(phyl_nocu, 5)[taxarank])
}
ntaxa(aggregate_taxa(ps_rarep_nocu,level = "Phylum"))

ps_rarep_yescu <- ps_yescu
for (taxarank in 1:5) {
  ps_rarep_yescu <- subset_taxa(ps_rarep_yescu, Phylum != top_taxa(phyl_yescu, 5)[taxarank])
}
ntaxa(aggregate_taxa(ps_rarep_yescu,level = "Phylum"))

print(paste(ntaxa(phyl_yescu) - ntaxa(aggregate_taxa(ps_rarep_yescu,level = "Phylum")),"phylum have been deleted to get a new file containing", ntaxa(aggregate_taxa(ps_rarep_nocu,level = "Phylum")), "rare phyla in non-copper-containing samples and", ntaxa(aggregate_taxa(ps_rarep_yescu,level = "Phylum")), "rare phyla in copper-containing samples"))

# 3.8.2.2) Keep 4 most abundant Phyla

ps_domp_nocu <- subset_taxa(ps_nocu, Phylum == top_taxa(phyl_nocu, 4))
ntaxa(aggregate_taxa(ps_domp_nocu,level = "Phylum"))
get_taxa_unique(ps_domp_nocu,"Phylum")

print(paste(ntaxa(phyl_nocu) - ntaxa(aggregate_taxa(ps_domp_nocu,level = "Phylum")),"phyla have been deleted to get a new file containing", get_taxa_unique(ps_domp_nocu,"Phylum")[1],",",get_taxa_unique(ps_domp_nocu,"Phylum")[2],",",get_taxa_unique(ps_domp_nocu,"Phylum")[3],"and",get_taxa_unique(ps_domp_nocu,"Phylum")[4], "dominant phyla in non-copper-containing samples"))

ps_domp_yescu <- subset_taxa(ps_yescu, Phylum == top_taxa(phyl_yescu, 4))
ntaxa(aggregate_taxa(ps_domp_yescu,level = "Phylum"))
get_taxa_unique(ps_domp_yescu,"Phylum")

print(paste(ntaxa(phyl_yescu) - ntaxa(aggregate_taxa(ps_domp_yescu,level = "Phylum")),"phyla have been deleted to get a new file containing", get_taxa_unique(ps_domp_yescu,"Phylum")[1],",",get_taxa_unique(ps_domp_yescu,"Phylum")[2],",",get_taxa_unique(ps_domp_yescu,"Phylum")[3],"and",get_taxa_unique(ps_domp_yescu,"Phylum")[4], "dominant phyla in non-copper-containing samples"))

# 3.8.2.3) Generate new phyloseq file excluding 5 most abundant gena

ps_rareg_nocu <- ps_nocu
for (taxarank in 1:5) {
  ps_rareg_nocu <- subset_taxa(ps_rareg_nocu, Genus != top_taxa(gen_nocu, 5)[taxarank])
}
ntaxa(aggregate_taxa(ps_rareg_nocu,level = "Genus"))

ps_rareg_yescu <- ps_yescu
for (taxarank in 1:5) {
  ps_rareg_yescu <- subset_taxa(ps_rareg_yescu, Genus != top_taxa(gen_yescu, 5)[taxarank])
}
ntaxa(aggregate_taxa(ps_rareg_yescu,level = "Genus"))

print(paste(ntaxa(gen_yescu) - ntaxa(aggregate_taxa(ps_rareg_yescu,level = "Genus")),"gena have been deleted to get a new file containing", ntaxa(aggregate_taxa(ps_rareg_nocu,level = "Genus")), "rare gena in non-copper-containing samples and", ntaxa(aggregate_taxa(ps_rareg_yescu,level = "Genus")), "rare gena in copper-containing samples"))

# 3.8.2.4) Keep 4 most abundant gena

ps_domg_nocu <- subset_taxa(ps_nocu, Genus == top_taxa(gen_nocu, 4))
ntaxa(aggregate_taxa(ps_domg_nocu,level = "Genus"))
get_taxa_unique(ps_domg_nocu,"Genus")

print(paste(ntaxa(gen_nocu) - ntaxa(aggregate_taxa(ps_domg_nocu,level = "Genus")),"gena have been deleted to get a new file containing", get_taxa_unique(ps_domg_nocu,"Genus")[1],",",get_taxa_unique(ps_domg_nocu,"Genus")[2],"and",get_taxa_unique(ps_domg_nocu,"Genus")[3], "dominant gena in non-copper-containing samples"))

ps_domg_yescu <- subset_taxa(ps_yescu, Genus == top_taxa(gen_nocu, 4))
ntaxa(aggregate_taxa(ps_domg_yescu,level = "Genus"))
get_taxa_unique(ps_domg_yescu,"Genus")

print(paste(ntaxa(gen_yescu) - ntaxa(aggregate_taxa(ps_domg_yescu,level = "Genus")),"gena have been deleted to get a new file containing", get_taxa_unique(ps_domg_yescu,"Genus")[1],",",get_taxa_unique(ps_domg_yescu,"Genus")[2],"and",get_taxa_unique(ps_domg_yescu,"Genus")[3], "dominant gena in copper-containing samples"))

# Only Pseudomonas

ps_pseudo_nocu <- subset_taxa(ps_nocu, Genus == "Pseudomonas")
ntaxa(ps_pseudo_nocu)
ps_pseudo_yescu <- subset_taxa(ps_yescu, Genus == "Pseudomonas")
ntaxa(ps_pseudo_yescu)

# 3.9) Auxiliary objects and functions

# Removing specific elements on plots

pr <-  theme_prism(base_size = 14) + theme(plot.title = element_blank(), plot.subtitle=element_blank())
ny <- rremove("y.axis") + rremove("ylab") + rremove("y.ticks") + rremove("y.text")
nx <- rremove("x.text")

#Color palettes

# Treatments
tre_col <- c('No' = "#4C4C4C", 
             'Hoagland'= "#0F80FF", 
             'AMH3-8' = "#62C7FA", 
             'No_plant' = "#A1A3A8")

tr_lab <- c('No' = 'Control',
            'Hoagland'= 'Hoagland',
            'AMH3-8' = 'AMH3-8',
            'No_plant' = 'Suelo NR')

# Copper
cop_col <- c('FALSE' = "#028202",
             'TRUE' = "#DE7920")
cop_lab <-c('FALSE' = 'Sin cobre',
            'TRUE' = 'Con cobre')

# Phylum

phy_col <- c('Acidobacteria' = "#7D3560",
             'Actinobacteria' = "#9D654C",
             'Bacteroidetes' = "#616161",
             'Chloroflexi' = "#148F77",
             'Cyanobacteria' = "#989c46",
             'Deinococcus-Thermus' = "#E0D65F",
             'Entotheonellaeota' = "#A5A5A5",
             'Epsilonbacteraeota' = "#42147C",
             'Firmicutes' = "#4E7705",
             'Gemmatimonadetes' = "#F6A04D",
             'Nitrospirae' = "#d64f7a",
             'Patescibacteria' = "#b55eb2",
             'Proteobacteria' = "#D70000",
             'Verrucomicrobia' = "#0c6a99")
phy_col2 <- setNames(phy_col, paste0("p__", names(phy_col)))
phy_lab <- c('Acidobacteria' = "Acidobacteria",
             'Actinobacteria' = "Actinobacteria",
             'Bacteroidetes' = "Bacteroidetes",
             'Chloroflexi' = "Chloroflexi",
             'Cyanobacteria' = "Cyanobacteria",
             'Deinococcus-Thermus' = "Deinococcus-Thermus",
             'Entotheonellaeota' = "Entotheonellaeota",
             'Epsilonbacteraeota' = "Epsilonbacteraeota",
             'Firmicutes' = "Firmicutes",
             'Gemmatimonadetes' = "Gemmatimonadetes",
             'Nitrospirae' = "Nitrospirae",
             'Patescibacteria' = "Patescibacteria",
             'Proteobacteria' = "Proteobacteria",
             'Verrucomicrobia' = "Verrucomicrobia")
phy_lab2 <- setNames(phy_lab, paste0("p__", names(phy_lab)))
phy_lab_r <- c('Acidobacteria' = "Aci.",
               'Actinobacteria' = "Act.",
               'Bacteroidetes' = "Bac.",
               'Chloroflexi' = "Chl.",
               'Cyanobacteria' = "Cya.",
               'Deinococcus-Thermus' = "D-T.",
               'Entotheonellaeota' = "Ent.",
               'Epsilonbacteraeota' = "Eps.",
               'Firmicutes' = "Fir.",
               'Gemmatimonadetes' = "Gem.",
               'Nitrospirae' = "Nit.",
               'Patescibacteria' = "Pat.",
               'Proteobacteria' = "Prot.",
               'Verrucomicrobia' = "Ver.")

# Gradient palettes

green_pal <- colorRampPalette(c("#F8FBF7", "#027F01"))(9)

orange_pal <- colorRampPalette(c("#FCFAF6", "#DC6D0B"))(9)

## 4) Exploratory ecologic analysis

# 4.1) Diversity analysis

# 4.1.1) Alpha diversity analysis

# 4.1.2) Boxplot diversity index vs variable

compar <- list( c("No", "No_plant"), c("No", "Hoagland"), c("No", "AMH3-8"))
compar.r <- list(c("No", "Hoagland"), c("No", "AMH3-8"),c("AMH3-8", "Hoagland"))
compar.all <- list( c("No", "No_plant"), c("No", "Hoagland"), c("No", "AMH3-8"))

adivplot <-  plot_richness(ps_rar, x="Treatment", measures="Shannon")+geom_boxplot(aes(fill=meta(ps_rar)$Fertilization))+theme_bw()+
  stat_compare_means(ref.group = "all", label = "p.signif", bracket.size = 0.7, fontface = "bold", size = 4.5) + 
  scale_color_manual(
    values=c("#4C4C4C", 
             "#0F80FF", 
             "#62C7FA", 
             "#A1A3A8"),
    aesthetics = "fill",
    labels = c('Control','Hoagland','AMH3-8','Suelo NR'),
    name = "") + # labels
  scale_x_discrete(labels= c('Control','Hoagland','AMH3-8','Suelo NR')) + # labels on X axis ticks
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
        axis.text = element_text(face = "bold", size = 12), # Axis text format
        axis.title = element_text(face = "bold", size = 12), # Axis title format 
        plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
        strip.background = element_blank(), # Remove strip
        strip.text = element_blank(),
        legend.text = element_text(face = "bold", size = 12), 
        axis.text.x = element_blank(), # remove X labels
        ) +
  ylab("Diversidad de Shannon-Weaver") +
  xlab("") +
  ylim(c(3.5,4.5))

adiv_nocu <- plot_richness(ps_nocu, x="Fertilization", measures="Shannon")+geom_boxplot(aes(fill=meta(ps_nocu)$Fertilization))+theme_bw()+
  scale_color_manual(
    values=c("#4C4C4C", 
             "#0F80FF", 
             "#62C7FA", 
             "#A1A3A8"),
    aesthetics = "fill",
    labels = c('Control','Hoagland','AMH3-8','Bulk Soil'),
    name = "") + # labels
  scale_x_discrete(labels= c('Control','Hoagland','AMH3-8','Bulk Soil')) + # labels on X axis ticks
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
        axis.text = element_text(face = "bold", size = 14), # Axis text format
        axis.title = element_text(face = "bold", size = 18), # Axis title format 
        plot.title = element_blank(), # Plot title format
        strip.background = element_blank(), # Remove strip
        strip.text = element_blank(),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_blank() # remove X labels
        ) +
  stat_compare_means(comparisons = compar, label = "p.signif", bracket.size = 0.7, fontface = "bold", size = 4.5) +
  ylab("Shannon diversity") +
  xlab("No Copper") + 
  ylim(c(3.5,4.5)) +
  ggtitle("")

adiv_yescu <- plot_richness(ps_yescu, x="Fertilization", measures="Shannon")+geom_boxplot(aes(fill=meta(ps_yescu)$Fertilization))+theme_bw()+
  scale_color_manual(
    values=c("#4C4C4C", 
             "#0F80FF", 
             "#62C7FA", 
             "#A1A3A8"),
    aesthetics = "fill",
    labels = c('Control','Hoagland','AMH3-8','Bulk Soil'),
    name = "") + # labels
  scale_x_discrete(labels= c('Control','Hoagland','AMH3-8','Bulk Soil')) + # labels on X axis ticks
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
        axis.text = element_text(face = "bold", size = 14), # Axis text format
        axis.title = element_text(face = "bold", size = 18), # Axis title format 
        plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
        strip.background = element_blank(), # Remove strip
        strip.text = element_blank(),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_blank() # remove X labels
        ) +
  stat_compare_means(comparisons = compar, label = "p.signif", bracket.size = 0.7, fontface = "bold", size = 4.5) +
  ylab("Shannon Diversity") +
  xlab("Copper") + 
  ylim(c(3.5,4.5)) +
  ggtitle("")

adiv_rare <- plot_richness(ps_rar_rare, x="Fertilization", measures="Shannon")+xlab("")+geom_boxplot(aes(fill=meta(ps_rar_rare)$Fertilization))+theme_bw()+ggtitle("Global (ASV's raros)")+
  scale_color_manual(
    values=c("#4C4C4C", 
             "#0F80FF", 
             "#62C7FA", 
             "#A1A3A8"),
    aesthetics = "fill",
    labels = c('Control','Hoagland','AMH3-8','Suelo NR'),
    name = "") + # labels
  scale_x_discrete(labels= c('Control','Hoagland','AMH3-8','Suelo NR')) + # labels on X axis ticks
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
        axis.text = element_text(face = "bold", size = 12), # Axis text format
        axis.title = element_text(face = "bold", size = 12), # Axis title format 
        plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
        strip.background = element_blank(), # Remove strip
        strip.text = element_blank()) +
  # ylim(c(0,5))+
  stat_compare_means(comparisons = compar, label = "p.signif") +
  ylab("Diversidad de Shannon") +
  xlab("") + 
  ylim(c(3.5,5))

adiv_nocu_rare <- plot_richness(ps_nocu_rare, x="Fertilization", measures="Shannon")+xlab("")+geom_boxplot(aes(fill=meta(ps_nocu_rare)$Fertilization))+theme_bw()+ggtitle("Sin cobre (ASV's raros)")+
  scale_color_manual(
    values=c("#4C4C4C", 
             "#0F80FF", 
             "#62C7FA", 
             "#A1A3A8"),
    aesthetics = "fill",
    labels = c('Control','Hoagland','AMH3-8','Suelo NR'),
    name = "") + # labels
  scale_x_discrete(labels= c('Control','Hoagland','AMH3-8','Suelo NR')) + # labels on X axis ticks
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
        axis.text = element_text(face = "bold", size = 12), # Axis text format
        axis.title = element_text(face = "bold", size = 12), # Axis title format 
        plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
        strip.background = element_blank(), # Remove strip
        strip.text = element_blank()) +
  # ylim(c(0,5))+
  stat_compare_means(comparisons = compar, label = "p.signif") +
  ylab("Diversidad de Shannon") +
  xlab("") + 
  ylim(c(3.5,5))

adiv_yescu_rare <- plot_richness(ps_yescu_rare, x="Fertilization", measures="Shannon")+xlab("")+geom_boxplot(aes(fill=meta(ps_yescu_rare)$Fertilization))+theme_bw()+ggtitle("Con cobre (ASV's raros)")+
  scale_color_manual(
    values=c("#4C4C4C", 
             "#0F80FF", 
             "#62C7FA", 
             "#A1A3A8"),
    aesthetics = "fill",
    labels = c('Control','Hoagland','AMH3-8','Suelo NR'),
    name = "") + # labels
  scale_x_discrete(labels= c('Control','Hoagland','AMH3-8','Suelo NR')) + # labels on X axis ticks
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
        axis.text = element_text(face = "bold", size = 12), # Axis text format
        axis.title = element_text(face = "bold", size = 12), # Axis title format 
        plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
        strip.background = element_blank(), # Remove strip
        strip.text = element_blank()) +
  # ylim(c(0,5))+
  stat_compare_means(comparisons = compar, label = "p.signif") +
  ylab("Diversidad de Shannon") +
  xlab("") + 
  ylim(c(3.5,5))

# adiv_r_nocu <- plot_richness(ps_r_nocu, x="Fertilization", measures="Shannon")+geom_boxplot(aes(fill=meta(ps_r_nocu)$Fertilization))+theme_bw()+
#   scale_color_manual(
#     values=c("#4C4C4C", 
#              "#0F80FF", 
#              "#62C7FA", 
#              "#A1A3A8"),
#     aesthetics = "fill",
#     labels = c('Control','Hoagland','AMH3-8','Suelo NR'),
#     name = "") + # labels
#   scale_x_discrete(labels= c('Control','Hoagland','AMH3-8','Suelo NR')) + # labels on X axis ticks
#   theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
#         axis.text = element_text(face = "bold", size = 12), # Axis text format
#         axis.title = element_text(face = "bold", size = 12), # Axis title format 
#         plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
#         strip.background = element_blank(), # Remove strip
#         strip.text = element_blank(),
#         legend.text = element_text(face = "bold", size = 12),
#         axis.text.x = element_blank() # remove X labels
#   ) +
#   stat_compare_means(comparisons = compar.r, label = "p.signif", bracket.size = 0.7, fontface = "bold", size = 4.5) +
#   ylab("Diversidad de Shannon") +
#   xlab("Sin cobre") + 
#   ylim(c(3.5,4.5))
# 
# adiv_r_yescu <- plot_richness(ps_r_yescu, x="Fertilization", measures="Shannon")+geom_boxplot(aes(fill=meta(ps_r_yescu)$Fertilization))+theme_bw()+
#   scale_color_manual(
#     values=c("#4C4C4C", 
#              "#0F80FF", 
#              "#62C7FA", 
#              "#A1A3A8"),
#     aesthetics = "fill",
#     labels = c('Control','Hoagland','AMH3-8','Suelo NR'),
#     name = "") + # labels
#   scale_x_discrete(labels= c('Control','Hoagland','AMH3-8','Suelo NR')) + # labels on X axis ticks
#   theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
#         axis.text = element_text(face = "bold", size = 12), # Axis text format
#         axis.title = element_text(face = "bold", size = 12), # Axis title format 
#         plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
#         strip.background = element_blank(), # Remove strip
#         strip.text = element_blank(),
#         legend.text = element_text(face = "bold", size = 12),
#         axis.text.x = element_blank() # remove X labels
#   ) +
#   stat_compare_means(comparisons = compar.r, label = "p.signif", bracket.size = 0.7, fontface = "bold", size = 4.5) +
#   ylab("Diversidad de Shannon") +
#   xlab("Con cobre") + 
#   ylim(c(3.5,4.5))
# 
# ggarrange(adivplot,adiv_nocu + rremove("ylab") ,adiv_yescu+ rremove("ylab") ,adiv_rare,adiv_nocu_rare+ rremove("ylab") ,adiv_yescu_rare+ rremove("ylab") , legend ="none", ncol=3, nrow = 2)



# 4.1.2) Beta distance analysis

# See available analysis for beta distance measurement
unlist(distanceMethodList)

# Ordination analysis
# wunifrac , X r X , me , co, g, 19, rlb
# dpcoa , c , wb , r
set.seed(1)
dmet <- "r"
ordmet <- "PCoA"
ord <- ordinate(ps_rar, ordmet, dmet)
ord_nocu <- ordinate(ps_nocu, ordmet, dmet)
ord_yescu <- ordinate(ps_yescu, ordmet, dmet)

ord_nocu_rare <- ordinate(ps_nocu_rare, ordmet,  dmet)
ord_yescu_rare <- ordinate(ps_yescu_rare, ordmet, dmet)
#sample_variables(ps_rar)
# 4.1.2.1) Scatterplot

ord_plot <- plot_ordination(ps_rar, ord, type="samples", color = "Copper", title="Principle Coordinate Analysis (PCoA): ASV level") +
  geom_point(size = 5)+ 
  stat_ellipse(type = "t") +
  theme_bw() +
  scale_color_manual(name = "Suelo", values = c("#DE7920", 
                                               "#028202"), #paleta de colores para con cobre, sin cobre
                    labels = c('Sin cobre','Con cobre'))
ord_plot_nocu <- plot_ordination(ps_nocu, ord_nocu, type="samples", color = "Fertilization", title="Sin cobre") +
  geom_point(size = 4) + #Standard colors for AMH3-8, Control, Hoagland and Bulk soil
  stat_conf_ellipse(level = 0.92) +
  theme_bw() +
  scale_color_manual(values=c("#4C4C4C", 
                              "#0F80FF", 
                              "#62C7FA", 
                              "#A1A3A8"),
                     labels= c('Control','Hoagland','AMH3-8','Suelo NR')) + # labels
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
        axis.text = element_text(face = "bold", size = 12), # Axis text format
        axis.title = element_text(face = "bold", size = 12), # Axis title format 
        plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
        # legend.key.size = unit(0.5, "cm"), # Space between legend elements
        legend.text = element_text(size = 12), # Legend text format
        legend.title = element_blank()) + # Remove legend title 
        xlab("PCo 1 (31,9%)") + # Título eje X
        ylab("PCo 2 (22,5%)") # Título eje Y
ord_plot_yescu <- plot_ordination(ps_yescu, ord_yescu, type="samples", color = "Fertilization", title="Con cobre") +
  geom_point(size = 4) + 
  stat_conf_ellipse(level = 0.92) +
  theme_bw() +
  scale_color_manual(values=c("#4C4C4C", 
                              "#0F80FF", 
                              "#62C7FA", 
                              "#A1A3A8"),
                     labels= c('Control','Hoagland','AMH3-8','Suelo NR')) + # labels
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
        axis.text = element_text(face = "bold", size = 12), # Axis text format
        axis.title = element_text(face = "bold", size = 12), # Axis title format 
        plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
        # legend.key.size = unit(0.5, "cm"), # Space between legend elements
        legend.text = element_text(size = 12), # Legend text format
        legend.title = element_blank()) + # Remove legend title
        xlab("PCo 1 (39,4%)") + # Título eje X
        ylab("PCo 2 (27,2%)") # Título eje Y
ord_plot_nocu_rare <- plot_ordination(ps_nocu_rare, ord_nocu_rare, type="sample", color = "Fertilization", title="Sin cobre (ASVs <1% abundancia)") +
  geom_point(size = 5)+ 
  stat_conf_ellipse(level = 0.92) +
  theme_bw()
ord_plot_yescu_rare <- plot_ordination(ps_yescu_rare, ord_yescu_rare, type="sample", color = "Fertilization", title="Con cobre (ASVs <1% abundancia)") +
  geom_point(size = 5)+ 
  stat_conf_ellipse(level = 0.92) +
  theme_bw()

ggarrange(ord_plot_nocu,ord_plot_yescu, nrow=1, ncol=2, legend = "bottom", common.legend = TRUE, labels = "AUTO", font.label = list(size = 18))
## Export as 6 x 8 in pdf
# ggarrange(ord_plot, "" , ord_plot_nocu,ord_plot_yescu, ord_plot_nocu_rare, ord_plot_yescu_rare, nrow=3, ncol=2)
unlist(distanceMethodList)

# Separando por phylum y mostrando en 4 filas:
# ord_plot + facet_wrap(~Phylum, 4)

# 4.2) Microbiome structure

# Access abundances

# abs_ab <- abundances(ps_rar, "identity")
# 
# rel_ab <- abundances(ps_rar, "compositional")
# 
# norm_ab <- abundances(ps_rar, "Z")

# 4.3) Intersectional composition analysis

# 4.3.1) Upset plot

# Generate upset data frame

vnus <- get_upset(obj=ps_bac2, factorNames="Treatment")

# Generate metadata dataframe

upset_meta <- meta(ps_bac2)[as.character(seq(3,24,3)),1:4] 
upset_meta$Copper <- as.character(upset_meta$Copper)
upset_meta$Fertilization <- as.character(upset_meta$Fertilization)
upset_meta$set <- as.character(upset_meta$Treatment)
upset_meta[upset_meta["Fertilization"] == "No", "Fertilization"] = "Control"
upset_meta[upset_meta["Fertilization"] == "No_plant", "Fertilization"] = "Bulk soil"

# Metadata to Spanish

upset_meta_es <- upset_meta
colnames(upset_meta_es)[3:4] <- c("Cobre","Fertilización")
upset_meta_es[upset_meta_es["Fertilización"] == "Bulk soil", "Fertilización"] = "Suelo NR"

# Lists of intersections

ilist <- list(as.character(seq(1,8,1)), as.character(seq(1,4,1)), as.character(seq(5,8,1)), as.character(seq(1,3,1)), as.character(seq(5,7,1)))
nlist <- list("1", "2", "3")
ylist <- list("5", "6", "7")

# Merge taxa table into the Upset data frame

ustx <- cbind(vnus,rownames(vnus))
colnames(ustx)[9] <- "asv"

usty <- cbind(phyloseq::tax_table(ps_bac2),rownames(phyloseq::tax_table(ps_bac2)))
colnames(usty)[8] <- "asv"

ustax <- merge(ustx, usty, by = "asv", all.x = T)
rownames(ustax) <- ustax[,"asv"]
ustax <- ustax[,-1] # this is the final data frame object that will be used as input for the upset function

upsetbars <- upset(ustax
                   , intersect = sort(unique(as.vector(sample_data(ps_bac2)$Treatment)), decreasing = T)
                   , sort_sets = FALSE
                   , name = ""
                   , sort_intersections=FALSE
                   , intersections = ilist
                   , base_annotations = list('Intersection size'=intersection_size(mode='inclusive_intersection'
                                                                                    , mapping=aes(fill=Phylum)
                                                                                   , colour="black"
                                                                                   , size = 0.6
                                                                                   , text=list(vjust=-8, color="black")
                                                                                   , show.legend = T) + 
                                               scale_fill_manual(values=phy_col) +
                                               ylim(0,250) + 
                                               ylab('Conserved ASVs') + 
                                               theme_prism(base_size = 12) +
                                               theme(legend.position = "right",
                                                     legend.key.size = unit(1, "cm"), 
                                                     legend.text = element_text(size=12)) + 
                                               rremove("xlab") + 
                                               rremove("x.text")
                                             )
                   , set_sizes=upset_set_size() + ylab('Mean ASVs')
                   , matrix =  (intersection_matrix() + scale_y_discrete(position = "right"))
                   , stripes = scales::alpha(c(rep("#DE7920",4),rep("#028202",4)),0.3)
                   , labeller = ggplot2::as_labeller(c("1" = "Control"
                                                       , "2" = "Hoagland"
                                                       , "3" = "AMH3-8"
                                                       , "4" = "Bulk Soil"
                                                       , "5" = "Control"
                                                       , "6" = "Hoagland"
                                                       , "7" = "AMH3-8"
                                                       , "8" = "Bulk Soil"
                                                       )
                                                     )
                   , queries=list(
                                   upset_query(set='1', fill='#4C4C4C', color='#4C4C4C')
                                  , upset_query(set='2', fill='#0F80FF', color='#0F80FF')
                                  , upset_query(set='3', fill='#62C7FA', color='#62C7FA')
                                  , upset_query(set='4', fill='#A1A3A8', color='#A1A3A8')
                                  , upset_query(set='5', fill='#4C4C4C', color='#4C4C4C')
                                  , upset_query(set='6', fill='#0F80FF', color='#0F80FF')
                                  , upset_query(set='7', fill='#62C7FA', color='#62C7FA')
                                  , upset_query(set='8', fill='#A1A3A8', color='#A1A3A8')
                                  )
                   , wrap = T
                   ) + ggtitle('') + theme(plot.title = element_text(face="bold", hjust=0, size = 18))

upset_nocu <- upset(ustax
                   , intersect = sort(as.character(seq(1,3,1)), decreasing = T)
                   , sort_sets = FALSE
                   , name = element_blank()
                   , sort_intersections=FALSE
                   , intersections = nlist
                   , base_annotations = list('Intersection size'=intersection_size(mapping=aes(fill=Phylum)
                                                                                   , colour="black"
                                                                                   , size = 0.6
                                                                                   , text=list(vjust=-2.5, color="black")
                                                                                   , show.legend = F) + 
                                               scale_fill_manual(values=phy_col) +
                                               ylim(0,50) + 
                                               ylab('Unique ASVs\n Absolute')  + 
                                               theme_prism(base_size = 12) + 
                                               rremove("xlab") + 
                                               rremove("x.text")
                                             )
                   , set_sizes= F
                   , labeller = ggplot2::as_labeller(c("1" = "Control"
                                                       , "2" = "Hoagland"
                                                       , "3" = "AMH3-8"
                                                       , "4" = "Bulk Soil"
                                                       , "5" = "Control"
                                                       , "6" = "Hoagland"
                                                       , "7" = "AMH3-8"
                                                       , "8" = "Bulk Soil"
                                                       )
                                                     )
                   # , themes = theme_prism()
                   , stripes = scales::alpha("#028202",0.3)
                   , queries=list(
                     upset_query(set='1', fill='#4C4C4C', color='#4C4C4C')
                     , upset_query(set='2', fill='#0F80FF', color='#0F80FF')
                     , upset_query(set='3', fill='#62C7FA', color='#62C7FA')
                     , upset_query(set='4', fill='#A1A3A8', color='#A1A3A8')
                   )
                   , annotations =list(
                     'Unique ASVs\n Relative (%)'=list(
                       aes=aes(x=intersection, fill=Phylum),
                       geom=list(
                         geom_bar(stat='count', position='fill', na.rm=T, colour="black", size = 0.6, show.legend = F), # Color of the bars
                         geom_text(
                           aes(
                             label=after_stat(count), # Labels for the bars
                             color=ifelse(after_stat(count)>median(after_stat(count)), 'show', 'hide')
                           ),
                           stat='count',
                           position=position_fill(vjust = .5)
                         ),
                         scale_y_continuous(labels=scales::percent_format()),
                         scale_color_manual(values=c('show'='black', 'hide'='transparent'), guide='none')
                         , theme_prism(base_size = 12)
                         , rremove("xlab")
                         , rremove("x.text")
                         , rremove("x.ticks") 
                         , rremove("x.axis") 
                         , scale_fill_manual(values=phy_col)
                       )
                     )
                   )
                   , wrap = T
                   ) + ggtitle('No Copper') + ylab("Unique ASVs") + theme(plot.title = element_text(face="bold", hjust= 0.5, size = 18))

upset_yescu <- upset(ustax
                    , intersect = sort(as.character(seq(5,7,1)), decreasing = T)
                    , sort_sets = FALSE
                    , name = element_blank()
                    , sort_intersections=FALSE
                    , intersections = ylist
                    , base_annotations = list('Intersection size'=intersection_size(mapping=aes(fill=Phylum)
                                                                                    , colour="black"
                                                                                    , size = 0.6
                                                                                    , text=list(vjust=-2.5, color="black")
                                                                                    , show.legend = F) + 
                                                scale_fill_manual(values=phy_col) +
                                                ylim(0,50) + 
                                                ylab('Unique ASVs\n Absolute')  + 
                                                theme_prism(base_size = 12) + 
                                                rremove("xlab") + 
                                                rremove("x.text")
                    )
                    , set_sizes= F
                    , labeller = ggplot2::as_labeller(c("1" = "Control"
                                                        , "2" = "Hoagland"
                                                        , "3" = "AMH3-8"
                                                        , "4" = "Bulk Soil"
                                                        , "5" = "Control"
                                                        , "6" = "Hoagland"
                                                        , "7" = "AMH3-8"
                                                        , "8" = "Bulk Soil"
                                                        )
                                                      )
                    , stripes = scales::alpha("#DE7920",0.3)
                    , queries=list(
                      upset_query(set='5', fill='#4C4C4C', color='#4C4C4C')
                      , upset_query(set='6', fill='#0F80FF', color='#0F80FF')
                      , upset_query(set='7', fill='#62C7FA', color='#62C7FA')
                      , upset_query(set='8', fill='#A1A3A8', color='#A1A3A8')
                    )
                    , annotations =list(
                      'Unique ASVs\n Relative (%)'=list(
                        aes=aes(x=intersection, fill=Phylum),
                        geom=list(
                          geom_bar(stat='count', position='fill', na.rm=T, colour="black", size = 0.6, show.legend = F), # Color of the bars
                          geom_text(
                            aes(
                              label=after_stat(count), # Labels for the bars
                              color=ifelse(after_stat(count)>median(after_stat(count)), 'show', 'hide')
                            ),
                            stat='count',
                            position=position_fill(vjust = .5)
                          ),
                          scale_y_continuous(labels=scales::percent_format()),
                          scale_color_manual(values=c('show'='black', 'hide'='transparent'), guide='none')
                          , theme_prism(base_size = 12)
                          , rremove("xlab")
                          , rremove("x.text")
                          , rremove("x.ticks") 
                          , rremove("x.axis") 
                          , scale_fill_manual(values=phy_col)
                        )
                      )
                    )
                    , wrap = T
)+ ggtitle(label = 'Copper') + theme(plot.title = element_text(face="bold", hjust= 0.5, size = 18))


(upsetbars| (upset_nocu/upset_yescu))


# # 4.3.2) Venn diagram
# 
# # vnls <- get_vennlist(obj=ps_rar, factorNames="Treatment")
# 
# vnls_nocu <- get_vennlist(obj=ps_nocu, factorNames="Fertilization")
# 
# vnls_yescu <- get_vennlist(obj=ps_yescu, factorNames="Fertilization")
# 
# vnls_r_nocu <- get_vennlist(obj=ps_r_nocu, factorNames="Fertilization")
# 
# vnls_r_yescu <- get_vennlist(obj=ps_r_yescu, factorNames="Fertilization")
# 
# 
# # library(nVennR)
# # 
# # v_ <- plotVenn(vnls, v_, nCycles = 1e6, outFile= "venn.svg")
# # 
# # save(v_, file = "v_.Rdata")
# 
# venn_nocu <- venn.diagram(vnls_nocu,
#                           height=5,
#                           width=5, 
#                           filename=NULL, 
#                           fill=c("#4C4C4C", 
#                                  "#0F80FF", 
#                                  "#62C7FA", 
#                                  "#A1A3A8"), #Standard colors for Control, Hoagland, AMH3-8 and Bulk soil
#                           # # cat.col=c("#62C7FA", 
#                           #           # "#0F80FF", 
#                           #           # "#4C4C4C", 
#                           #           # "#A1A3A8"), #Standard colors for AMH3-8, Control, Hoagland and Bulk soil
#                           category.names = c("Control", "Hoagland","AMH3-8","Suelo NR"),
#                           print.mode = "raw",
#                           na = "remove",
#                           main = "Sin cobre",
#                           main.cex = 1.5,
#                           alpha = 0.8, 
#                           fontfamily = "Arial",
#                           main.fontfamily = "Arial",
#                           fontface = "bold",
#                           cat.fontfamily = "Arial",
#                           cex = 1.2,
#                           cat.cex = 1.3,
#                           cat.default.pos = "outer",
#                           cat.dist=0.21,
#                           margin = 0, 
#                           lwd = 3,
#                           lty ='dotted',
#                           imagetype = "tiff")
# venn_yescu <- venn.diagram(vnls_yescu,
#                            height=5,
#                            width=5, 
#                            filename=NULL, 
#                            fill=c("#4C4C4C", 
#                                   "#0F80FF", 
#                                   "#62C7FA", 
#                                   "#A1A3A8"), #Standard colors for Control, Hoagland, AMH3-8 and Bulk soil
#                            # # cat.col=c("#62C7FA", 
#                            #           # "#0F80FF", 
#                            #           # "#4C4C4C", 
#                            #           # "#A1A3A8"), #Standard colors for AMH3-8, Control, Hoagland and Bulk soil
#                            category.names = c("Control", "Hoagland","AMH3-8","Suelo NR"),
#                            print.mode = "raw",
#                            na = "remove",
#                            main = "Con cobre",
#                            main.cex = 1.5,
#                            alpha = 0.8, 
#                            fontfamily = "Arial",
#                            main.fontfamily = "Arial",
#                            fontface = "bold",
#                            cat.fontfamily = "Arial",
#                            cex = 1.2,
#                            cat.cex = 1.3,
#                            cat.default.pos = "outer",
#                            cat.dist=0.21,
#                            margin = 0, 
#                            lwd = 3,
#                            lty ='dotted',
#                            imagetype = "tiff")
# venn_r_nocu <- venn.diagram(vnls_r_nocu,
#                             height=5,
#                             width=5, 
#                             filename=NULL, 
#                             fill=c("#4C4C4C", 
#                                    "#0F80FF", 
#                                    "#62C7FA"), #Standard colors for Control, Hoagland and AMH3-8
#                             # # cat.col=c("#62C7FA", 
#                             #           # "#0F80FF", 
#                             #           # "#4C4C4C", 
#                             #           # "#A1A3A8"), #Standard colors for AMH3-8, Control, Hoagland and Bulk soil
#                             category.names = c("Control", "Hoagland","AMH3-8"),
#                             print.mode = "percent",
#                             na = "remove",
#                             main = "Sin cobre",
#                             main.cex = 1.5,
#                             alpha = 0.8, 
#                             fontfamily = "Arial",
#                             main.fontfamily = "Arial",
#                             fontface = "bold",
#                             cat.fontfamily = "Arial",
#                             cex = 1.2,
#                             cat.cex = 1.3,
#                             cat.default.pos = "outer",
#                             cat.dist=0.21,
#                             margin = 0, 
#                             lwd = 3,
#                             lty ='dotted',
#                             imagetype = "tiff")
# venn_r_yescu <- venn.diagram(vnls_r_yescu,
#                              height=5,
#                              width=5, 
#                              filename=NULL, 
#                              fill=c("#4C4C4C", 
#                                     "#0F80FF", 
#                                     "#62C7FA"), #Standard colors for Control, Hoagland and AMH3-8
#                              # # cat.col=c("#62C7FA", 
#                              #           # "#0F80FF", 
#                              #           # "#4C4C4C", 
#                              #           # "#A1A3A8"), #Standard colors for AMH3-8, Control, Hoagland and Bulk soil
#                              category.names = c("Control", "Hoagland","AMH3-8"),
#                              print.mode = "percent",
#                              sigdigs = 3,
#                              na = "remove",
#                              main = "Con cobre",
#                              main.cex = 1.5,
#                              alpha = 0.8, 
#                              fontfamily = "Arial",
#                              main.fontfamily = "Arial",
#                              fontface = "bold",
#                              cat.fontfamily = "Arial",
#                              cex = 1.2,
#                              cat.cex = 1.3,
#                              cat.default.pos = "outer",
#                              cat.dist=0.21,
#                              margin = 0, 
#                              lwd = 3,
#                              lty ='dotted',
#                              imagetype = "tiff")
# grid.arrange(venn_nocu, venn_yescu, venn_r_nocu, venn_r_yescu, ncol=2)
# 
# grid.arrange(venn_nocu, venn_yescu, ncol=2)

# 4.4) Microbiome taxonomic structure analysis

# 4.4.1) By at genera by phylum level

# Using 5 most abundant phyla as main colors
melt_full_dom <- create_color_dfs(melt_full,selected_groups = top_taxa(phyl, 5) 
                                  , cvd = TRUE)

# Replace grayshade palette
melt_full_dom$cdf[melt_full_dom$cdf["hex"] == "#616161", "hex"] = "#54533a"
melt_full_dom$cdf[melt_full_dom$cdf["hex"] == "#8B8B8B", "hex"] = "#8a896b"
melt_full_dom$cdf[melt_full_dom$cdf["hex"] == "#B7B7B7", "hex"] = "#b4b58f"
melt_full_dom$cdf[melt_full_dom$cdf["hex"] == "#D6D6D6", "hex"] = "#d4d49d"
melt_full_dom$cdf[melt_full_dom$cdf["hex"] == "#F5F5F5", "hex"] = "#ecf0a3"

# Replace blueshade palette by red
melt_full_dom$cdf[melt_full_dom$cdf["hex"] == "#098BD9", "hex"] = "#D70000"
melt_full_dom$cdf[melt_full_dom$cdf["hex"] == "#56B4E9", "hex"] = "#E93333"
melt_full_dom$cdf[melt_full_dom$cdf["hex"] == "#7DCCFF", "hex"] = "#FF5252"
melt_full_dom$cdf[melt_full_dom$cdf["hex"] == "#BCE1FF", "hex"] = "#FF7F7F"
melt_full_dom$cdf[melt_full_dom$cdf["hex"] == "#E7F4FF", "hex"] = "#FFA5A5"

# Generate a custom legend for this palette
phy_legend <-custom_legend(melt_full_dom$mdf, melt_full_dom$cdf, legend_key_size = 0.8, legend_text_size = 12)


# Using  5 "arbitrary" (significant differences) biomarker phyla as main colors
melt_full_diff <- create_color_dfs(melt_full,selected_groups = c("Actinobacteria", "Gemmatimonadetes", "Proteobacteria", "Firmicutes") , cvd = TRUE)

# Generate a custom legend for this palette
phy_diff_legend <-custom_legend(melt_full_diff$mdf, melt_full_diff$cdf)


# Plot with vertical bars

comp.genphy_full_v <- plot_microshades(melt_full_dom$mdf, melt_full_dom$cdf) + 
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0,0.55), 
    # trans='log2',
    expand = expansion(0)) +
  geom_hline(yintercept=0.01, linetype="dashed", color = "blue") + 
  theme_classic() +
  theme(legend.position = "none"
        , axis.text.x= element_text(face = "bold", 
                                    color = c('Control' = "#4C4C4C", 
                                               'Hoagland'= "#0F80FF", 
                                               'AMH3-8' = "#62C7FA", 
                                               'Suelo NR' = "#A1A3A8"), 
                                    angle = 45,
                                    size = 10,
                                    vjust = 0.85)
        , axis.text.y= element_text(face = "bold")
        , axis.title.y = element_text(face = "bold",
                                      size = 12)
        )  +
  scale_x_discrete(labels=c("Control" = "C",
                            "Hoagland" = "H",
                            "AMH3-8" = "A", 
                            "Suelo NR" = "NR")) +
  facet_grid2(vars(Copper), vars(factor(Phylum, levels=top_taxa(phyl,12)))
              , axes = "x"
              , remove_labels = "none"
              , labeller = labeller(.cols = phy_lab)
              ) + # Grid view of phylum by Copper
  # facet_nested_wrap(~Copper+Phylum, ncol = 22) + # Nested view
  theme(plot.margin = margin(6,20,6,6),
        strip.text.x = element_text(
          size = 10, face = "bold"
        ),
        strip.text.y = element_text(
          size = 12, face = "bold"
        )) +
  ylab("Abundancia relativa (%)") + xlab("")

ab_genphyl_v <- plot_grid(comp.genphy_full_v, phy_legend,  rel_widths = c(1, .1)) # best exported at 7 x 16 inch landscape

# Segmented

comp.genphy_full_v_1 <- plot_microshades(melt_full_dom$mdf, melt_full_dom$cdf) + 
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0,0.55), 
    # trans='log2',
    expand = expansion(0)) +
  geom_hline(yintercept=0.01, linetype="dashed", color = "#3147ad", size = 0.8, alpha = 0.7) + 
  theme_classic() +
  theme(legend.position = "none"
        , axis.text.x= element_text(face = "bold", 
                                    color = c('Control' = "#4C4C4C", 
                                              'Hoagland'= "#0F80FF", 
                                              'AMH3-8' = "#62C7FA", 
                                              'Suelo NR' = "#A1A3A8"), 
                                    angle = 45,
                                    size = 15,
                                    vjust = 0.85)
        , axis.text.y= element_text(size = 14, face = "bold")
        , axis.title.y = element_text(face = "bold",
                                      size = 16)
  )  +
  scale_x_discrete(labels=c("Control" = "C",
                            "Hoagland" = "H",
                            "AMH3-8" = "A", 
                            "Suelo NR" = "NR")) +
  facet_grid2(vars(Copper), vars(factor(Phylum, levels=top_taxa(phyl,5)))
              , axes = "x"
              , remove_labels = "none"
              , labeller = labeller(.cols = phy_lab)
  ) + # Grid view of phylum by Copper
  theme(plot.margin = margin(6,6,6,6),
        strip.text.x = element_text(
          size = 14, face = "bold"
        ),
        strip.text.y = element_text(
          size = 18, face = "bold"
        )
  ) +
  ylab("Abundancia relativa (%)") + xlab("")

comp.genphy_full_v_2 <- plot_microshades(subset(melt_full_dom$mdf, !(Phylum %in% top_taxa(phyl,5))), melt_full_dom$cdf) + 
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0,0.015), 
    expand = expansion(0)) +
  geom_hline(yintercept=0.01, linetype="dashed", color = "#3147ad", size = 0.8, alpha = 0.7) + 
  theme_classic() +
  theme(legend.position = "none"
        , axis.text.x= element_text(face = "bold", 
                                    color = c('Control' = "#4C4C4C", 
                                              'Hoagland'= "#0F80FF", 
                                              'AMH3-8' = "#62C7FA", 
                                              'Suelo NR' = "#A1A3A8"), 
                                    angle = 45,
                                    size = 15,
                                    vjust = 0.85)
        , axis.text.y= element_text(size = 14, face = "bold")
        , axis.title.y = element_text(face = "bold",
                                      size = 16)
  )  +
  scale_x_discrete(labels=c("Control" = "C",
                            "Hoagland" = "H",
                            "AMH3-8" = "A", 
                            "Suelo NR" = "NR")) +
  facet_grid2(vars(Copper), vars(factor(Phylum, levels=top_taxa(phyl,12)[6:12]))
              , axes = "x"
              , remove_labels = "none"
              , labeller = labeller(.cols = phy_lab)
  )+ 
  theme(plot.margin = margin(6,6,6,6),
        strip.text.x = element_text(
          size = 14, face = "bold"
        ),
        strip.text.y = element_text(
          size = 18, face = "bold"
        )
  ) +
  ylab("Abundancia relativa (%)") + xlab("")

phy_legend_2 <-custom_legend(melt_full_dom$mdf, melt_full_dom$cdf, legend_key_size = 1.1, legend_text_size = 18)

ab_genphyl_v_parts <- plot_grid(plot_grid(comp.genphy_full_v_1, comp.genphy_full_v_2,  nrow = 2), phy_legend_2,  rel_widths = c(1, .135))

# Divide in 2 stacked plots. The upper for dominant and the lower for rare phyla

comp.genphy_full_h <- plot_microshades(melt_full_dom$mdf, melt_full_dom$cdf) + 
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0,0.55), 
    # trans='log2',
    expand = expansion(0)) +
  geom_hline(yintercept=0.01, linetype="dashed", color = "#8B8B8B") + 
  theme_classic() +
  theme(legend.position = "none"
        , axis.text.x= element_text(face = "bold", 
                                    angle = 45,
                                    size = 10,
                                    vjust = 0.85)
        , axis.text.y= element_text(face = "bold", 
                                    color = c('Control' = "#4C4C4C", 
                                              'Hoagland'= "#0F80FF", 
                                              'AMH3-8' = "#62C7FA", 
                                              'Suelo NR' = "#A1A3A8"), 
                                    size = 10
                                    )
        , axis.title.x = element_text(face = "bold",
                                      size = 12)
  )  +
  scale_x_discrete(labels=c("Control" = "C",
                            "Hoagland" = "H",
                            "AMH3-8" = "A", 
                            "Suelo NR" = "NR")) +
  coord_flip() + 
  facet_grid2(vars(factor(Phylum, levels=top_taxa(phyl,12))), vars(Copper)
              , axes = "y"
              , remove_labels = "y"
              , labeller = labeller(.rows = phy_lab_r)
              )+ # Grid view of phylum by Copper
  theme(plot.margin = margin(6,20,6,6),
        strip.text.y = element_text(
          size = 10
        ), strip.text.x = element_text(
          size = 18, 
          face = "bold", 
          vjust = 2.2
        )
        , strip.background.x = element_blank()
        ) +
  ylab("Abundancia relativa (%)") + xlab("")
        

ab_genphyl_h <- plot_grid(comp.genphy_full_h, phy_legend,  rel_widths = c(1, .1)) # best exported at 8 x 16 inch landscape

# Divide in 2 stacked plots. The upper for dominant and the lower for rare phyla
# subset(melt_full_dom$mdf, (Phylum %in% top_taxa(phyl,5))) ## Subset dominant taxa
comp.genphy_full_h_1 <- plot_microshades(melt_full_dom$mdf, melt_full_dom$cdf) + 
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0,0.55), 
    # trans='log2',
    expand = expansion(0)) +
  geom_hline(yintercept=0.01, linetype="dashed", color = "#3147ad", size = 0.8, alpha = 0.7) + 
  theme_classic() +
  theme(legend.position = "none"
        , axis.text.x= element_text(face = "bold", 
                                    angle = 45,
                                    size = 10,
                                    vjust = 0.85)
        , axis.text.y= element_text(face = "bold", 
                                    color = c('Control' = "#4C4C4C", 
                                              'Hoagland'= "#0F80FF", 
                                              'AMH3-8' = "#62C7FA", 
                                              'Suelo NR' = "#A1A3A8"), 
                                    size = 10
        )
        , axis.title.x = element_text(face = "bold",
                                      size = 12)
  )  +
  scale_x_discrete(labels=c("Control" = "C",
                            "Hoagland" = "H",
                            "AMH3-8" = "A", 
                            "Suelo NR" = "BS")) +
  coord_flip() + 
  facet_grid2(vars(factor(Phylum, levels=top_taxa(phyl,5))), vars(Copper)
              , axes = "y"
              , remove_labels = "y"
              , labeller = labeller(.rows = phy_lab_r)
  )+ 
  theme(plot.margin = margin(6,6,0,23),
        strip.text.y = element_text(
          size = 11
        ) 
        # , strip.text.x = element_text(
        #   size = 17, 
        #   face = "bold", 
        #   vjust = 2.2
        # )
        , strip.text.x = element_blank() # no titles
        , strip.background.x = element_blank()
        , panel.spacing.x = unit(2, "lines")
  ) +
  ylab("") + xlab("")

comp.genphy_full_h_2 <- plot_microshades(subset(melt_full_dom$mdf, !(Phylum %in% top_taxa(phyl,5))), melt_full_dom$cdf) + 
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0,0.015), 
    # trans='log2',
    expand = expansion(0)) +
  geom_hline(yintercept=0.01, linetype="dashed", color = "#3147ad", size = 0.8, alpha = 0.7) + 
  theme_classic() +
  theme(legend.position = "none"
        , axis.text.x= element_text(face = "bold", 
                                    angle = 45,
                                    size = 10,
                                    vjust = 0.85)
        , axis.text.y= element_text(face = "bold", 
                                    color = c('Control' = "#4C4C4C", 
                                              'Hoagland'= "#0F80FF", 
                                              'AMH3-8' = "#62C7FA", 
                                              'Suelo NR' = "#A1A3A8"), 
                                    size = 10
        )
        , axis.title.x = element_text(face = "bold",
                                      size = 12)
  )  +
  scale_x_discrete(labels=c("Control" = "C",
                            "Hoagland" = "H",
                            "AMH3-8" = "A", 
                            "Suelo NR" = "BS")) +
  coord_flip() + 
  facet_grid2(vars(factor(Phylum, levels=top_taxa(phyl,12)[6:12])), vars(Copper)
              , axes = "y"
              , remove_labels = "y"
              , labeller = labeller(.rows = phy_lab_r)
  )+ 
  theme(plot.margin = margin(0,6,6,23),
        strip.text.y = element_text(
          size = 11
        ), strip.text.x = element_blank()
        , strip.background.x = element_blank()
        , panel.spacing.x = unit(2, "lines")
  ) +
  ylab("Relative abundance (%)") + xlab("")

ab_genphyl_h_parts <- plot_grid(plot_grid(comp.genphy_full_h_1, comp.genphy_full_h_2,  nrow = 2), phy_legend,  rel_widths = c(1, .1))

# 4.4.1) Microbiome structure graph at Phylum level

ps.phylum.rel <- ps_rar %>%
  aggregate_taxa(level = "Phylum") %>%  
  microbiome::transform(transform = "compositional")

p.comp.phylum.rel <- plot_composition(ps.phylum.rel,
                                           taxonomic.level = "Phylum",
                                           sample.sort = "Fertilization",
                                           x.label = "Fertilization",
                                           group_by = "Copper",
                                           average_by = "Treatment") +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", y = "Relative abundance",
       title = "Phylum microbiome structure",
       subtitle = "",
       caption = "") + 
  theme_bw()
  theme(axis.title= element_text(size=20, colour = "black"),
        axis.title.x=element_text(size=15,face = "bold",hjust = 0.5),
        axis.title.y=element_text(size=15,face = "bold",hjust = 0.5),
        axis.text.x=element_text(size=15, angle=90,hjust=0.5,vjust=0.2))
  
  ps.phylum.rel.rare <- ps_rar_rare %>%
    aggregate_taxa(level = "Phylum") %>%  
    microbiome::transform(transform = "compositional")
  
ps.phylum.rel.rare <- subset_taxa(ps.phylum.rel, abundances(ps.phylum.rel, "compositional") < 0.01)
  
  p.comp.phylum.rel <- plot_composition(ps.phylum.rel,
                                        taxonomic.level = "Phylum",
                                        sample.sort = "Fertilization",
                                        x.label = "Fertilization",
                                        group_by = "Copper",
                                        average_by = "Treatment") +
    guides(fill = guide_legend(ncol = 1)) +
    labs(x = "Samples", y = "Relative abundance",
         title = "Phylum microbiome structure",
         subtitle = "",
         caption = "") + 
    theme_bw()
  theme(axis.title= element_text(size=20, colour = "black"),
        axis.title.x=element_text(size=15,face = "bold",hjust = 0.5),
        axis.title.y=element_text(size=15,face = "bold",hjust = 0.5),
        axis.text.x=element_text(size=15, angle=90,hjust=0.5,vjust=0.2))


ps.phylum.rel.nocu <- ps_nocu %>%
  aggregate_taxa(level = "Phylum") %>%  
  microbiome::transform(transform = "compositional")

ps.phylum.rel.yescu <- ps_yescu %>%
  aggregate_taxa(level = "Phylum") %>%  
  microbiome::transform(transform = "compositional")

p.comp.phylum.rel.nocu <- plot_composition(ps.phylum.rel.nocu,
                               taxonomic.level = "Phylum",
                               sample.sort = "Fertilization",
                               x.label = "Fertilization",
                               average_by = "Fertilization") +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", y = "Relative abundance",
       title = "Phylum microbiome structure",
       subtitle = "No Copper",
       caption = "") + 
  
  theme(axis.title= element_text(size=20, colour = "black"),
        axis.title.x=element_text(size=15,face = "bold",hjust = 0.5),
        axis.title.y=element_text(size=15,face = "bold",hjust = 0.5),
        axis.text.x=element_text(size=15, angle=90,hjust=0.5,vjust=0.2))
p.comp.phylum.rel.yescu <- plot_composition(ps.phylum.rel.yescu,
                                            taxonomic.level = "Phylum",
                                            sample.sort = "Fertilization",
                                            x.label = "Fertilization",
                                            average_by = "Fertilization") +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", y = "Relative abundance",
       title = "Phylum microbiome structure",
       subtitle = "Yes Copper",
       caption = "") + 
  
  theme(axis.title= element_text(size=20, colour = "black"),
        axis.title.x=element_text(size=15,face = "bold",hjust = 0.5),
        axis.title.y=element_text(size=15,face = "bold",hjust = 0.5),
        axis.text.x=element_text(size=15, angle=90,hjust=0.5,vjust=0.2))

ps.dom.phylum.rel.nocu <- ps_nocu_dom %>%
  aggregate_taxa(level = "Phylum") %>%  
  microbiome::transform(transform = "compositional")

ps.dom.phylum.rel.yescu <- ps_yescu_dom %>%
  aggregate_taxa(level = "Phylum") %>%  
  microbiome::transform(transform = "compositional")

p.comp.dom.phylum.rel.nocu <- plot_composition(ps.dom.phylum.rel.nocu,
                                           taxonomic.level = "Phylum",
                                           sample.sort = "Fertilization",
                                           x.label = "Fertilization",
                                           average_by = "Fertilization") +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", y = "Relative abundance",
       title = "Dominant phyla microbiome structure",
       subtitle = "No Copper",
       caption = "") + 
  
  theme(axis.title= element_text(size=20, colour = "black"),
        axis.title.x=element_text(size=15,face = "bold",hjust = 0.5),
        axis.title.y=element_text(size=15,face = "bold",hjust = 0.5),
        axis.text.x=element_text(size=15, angle=90,hjust=0.5,vjust=0.2))
p.comp.dom.phylum.rel.yescu <- plot_composition(ps.dom.phylum.rel.yescu,
                                            taxonomic.level = "Phylum",
                                            sample.sort = "Fertilization",
                                            x.label = "Fertilization",
                                            average_by = "Fertilization") +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", y = "Relative abundance",
       title = "Dominant phyla microbiome structure",
       subtitle = "Yes Copper",
       caption = "") + 
  
  theme(axis.title= element_text(size=20, colour = "black"),
        axis.title.x=element_text(size=15,face = "bold",hjust = 0.5),
        axis.title.y=element_text(size=15,face = "bold",hjust = 0.5),
        axis.text.x=element_text(size=15, angle=90,hjust=0.5,vjust=0.2))

ps.rare.phylum.rel.nocu <- ps_nocu_rare %>%
  aggregate_taxa(level = "Phylum") %>%  
  microbiome::transform(transform = "compositional")

ps.rare.phylum.rel.yescu <- ps_yescu_rare %>%
  aggregate_taxa(level = "Phylum") %>%  
  microbiome::transform(transform = "compositional")

p.comp.rare.phylum.rel.nocu <- plot_composition(ps.rare.phylum.rel.nocu,
                                                taxonomic.level = "Phylum",
                                                sample.sort = "Fertilization",
                                                x.label = "Fertilization",
                                                average_by = "Fertilization") +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", y = "Relative abundance",
       title = "Rare phyla microbiome structure",
       subtitle = "No Copper",
       caption = "") + 
  
  theme(axis.title= element_text(size=20, colour = "black"),
        axis.title.x=element_text(size=15,face = "bold",hjust = 0.5),
        axis.title.y=element_text(size=15,face = "bold",hjust = 0.5),
        axis.text.x=element_text(size=15, angle=90,hjust=0.5,vjust=0.2))
p.comp.rare.phylum.rel.yescu <- plot_composition(ps.rare.phylum.rel.yescu,
                                                 taxonomic.level = "Phylum",
                                                 sample.sort = "Fertilization",
                                                 x.label = "Fertilization",
                                                 average_by = "Fertilization") +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", y = "Relative abundance",
       title = "Rare phyla microbiome structure",
       subtitle = "Yes Copper",
       caption = "") + 
  
  theme(axis.title= element_text(size=20, colour = "black"),
        axis.title.x=element_text(size=15,face = "bold",hjust = 0.5),
        axis.title.y=element_text(size=15,face = "bold",hjust = 0.5),
        axis.text.x=element_text(size=15, angle=90,hjust=0.5,vjust=0.2))

abundance_phylum_plot <- ggarrange(p.comp.dom.phylum.rel.nocu, p.comp.dom.phylum.rel.yescu,p.comp.rare.phylum.rel.nocu, p.comp.rare.phylum.rel.yescu, ncol=2, nrow=2)

#
# 4.4.2) Microbiome structure graph at Genus level
#

plot_composition(phyl,
                 taxonomic.level = "Phylum",
                 sample.sort = "Fertilization",
                 x.label = "Fertilization",
                 otu.sort = "abundance",
                 group_by = "Copper") +
  geom_col(position = position_dodge(width = 10)) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Sample", y = "Abundance",
       title = "Genus microbiome structure",
       subtitle = "",
       caption = "") + 
  
  theme(axis.title= element_text(size=20, colour = "black"),
        axis.title.x=element_text(size=15,face = "bold",hjust = 0.5),
        axis.title.y=element_text(size=15,face = "bold",hjust = 0.5),
        axis.text.x=element_text(size=15, angle=90,hjust=0.5,vjust=0.2))

ps.genus.rel.nocu <- ps_nocu %>%
  aggregate_taxa(level = "Genus") %>%  
  microbiome::transform(transform = "compositional")

ps.genus.rel.yescu <- ps_yescu %>%
  aggregate_taxa(level = "Genus") %>%  
  microbiome::transform(transform = "compositional")

p.comp.genus.rel.nocu <- plot_composition(ps.genus.rel.nocu,
                                    taxonomic.level = "Genus",
                                    sample.sort = "Sample_ID",
                                    x.label = "Fertilization") +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", y = "Relative abundance",
       title = "Genus microbiome structure",
       subtitle = "No Copper",
       caption = "") + 
  
  theme(axis.title= element_text(size=20, colour = "black"),
        axis.title.x=element_text(size=15,face = "bold",hjust = 0.5),
        axis.title.y=element_text(size=15,face = "bold",hjust = 0.5),
        axis.text.x=element_text(size=15, angle=90,hjust=0.5,vjust=0.2))
p.comp.genus.rel.yescu <- plot_composition(ps.genus.rel.yescu,
                                     taxonomic.level = "Genus",
                                     sample.sort = "Fertilization",
                                     x.label = "Fertilization") +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", y = "Relative abundance",
       title = "Genus microbiome structure",
       subtitle = "Yes Copper",
       caption = "") + 
  
  theme(axis.title= element_text(size=20, colour = "black"),
        axis.title.x=element_text(size=15,face = "bold",hjust = 0.5),
        axis.title.y=element_text(size=15,face = "bold",hjust = 0.5),
        axis.text.x=element_text(size=15, angle=90,hjust=0.5,vjust=0.2))

ps.rare.genus.rel.nocu <- ps_rareg_nocu %>%
  aggregate_taxa(level = "Genus") %>%  
  microbiome::transform(transform = "compositional")

ps.rare.genus.rel.yescu <- ps_rareg_yescu %>%
  aggregate_taxa(level = "Genus") %>%  
  microbiome::transform(transform = "compositional")

p.comp.rare.genus.rel.nocu <- plot_composition(ps.rare.genus.rel.nocu,
                                                taxonomic.level = "Genus",
                                                sample.sort = "Fertilization",
                                                x.label = "Fertilization") +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", y = "Relative abundance",
       title = "Rare genera microbiome structure",
       subtitle = "No Copper",
       caption = "") + 
  
  theme(axis.title= element_text(size=20, colour = "black"),
        axis.title.x=element_text(size=15,face = "bold",hjust = 0.5),
        axis.title.y=element_text(size=15,face = "bold",hjust = 0.5),
        axis.text.x=element_text(size=15, angle=90,hjust=0.5,vjust=0.2))
p.comp.rare.genus.rel.yescu <- plot_composition(ps.rare.genus.rel.yescu,
                                                 taxonomic.level = "Genus",
                                                 sample.sort = "Fertilization",
                                                 x.label = "Fertilization") +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", y = "Relative abundance",
       title = "Rare genera microbiome structure",
       subtitle = "Yes Copper",
       caption = "") + 
  
  theme(axis.title= element_text(size=20, colour = "black"),
        axis.title.x=element_text(size=15,face = "bold",hjust = 0.5),
        axis.title.y=element_text(size=15,face = "bold",hjust = 0.5),
        axis.text.x=element_text(size=15, angle=90,hjust=0.5,vjust=0.2))

abundance_genus_plot <- ggarrange(p.comp.genus.rel.nocu, p.comp.genus.rel.yescu, p.comp.rare.genus.rel.nocu, p.comp.rare.genus.rel.yescu, ncol=2, nrow = 2)

### Heatmap of communities

ps.phylum.rel <- ps_rar %>%
  aggregate_taxa(level = "Phylum") %>%  
  microbiome::transform(transform = "compositional")

hm_phylum <- plot_heatmap(ps.phylum.rel, taxa.label="Phylum")

ps.phylum.rel.nocu <- ps_nocu %>%
  aggregate_taxa(level = "Phylum") %>%  
  microbiome::transform(transform = "compositional")

hm_phylum_nocu <- plot_heatmap(ps.phylum.rel.nocu, sample.label = "Fertilization", sample.order = "Fertilization" ,  taxa.label="Phylum")

ps.phylum.rel.yescu <- ps_yescu %>%
  aggregate_taxa(level = "Phylum") %>%  
  microbiome::transform(transform = "compositional")

hm_phylum_yescu <- plot_heatmap(ps.phylum.rel.yescu, sample.label = "Fertilization", sample.order = "Fertilization" , taxa.label="Phylum")

ggarrange(hm_phylum_nocu, hm_phylum_yescu, ncol=2)


#
# Microbiome structure graph at Phylum level
#

ps.fam.rel <- ps_rar %>%
  aggregate_taxa(level = "Family") %>%  
  subset_taxa(Family != "Unknown_Family") %>%   
  subset_taxa(Family != "Unknown") %>%   
  microbiome::transform(transform = "compositional")

# plot_composition(ps.fam.rel,
#                  plot.type = "heatmap", 
#                  taxonomic.level = "Family",
#                   average_by = "Fertilization",
#                  x.label = "Fertilization", 
#                  otu.sort = "abundance") +
#   guides(fill = guide_legend(ncol = 1)) +
#   labs(x = "Samples", y = "Relative abundance",
#        title = "Rare genera microbiome structure",
#        subtitle = "No Copper",
#        caption = "")

# hm_fam <- plot_heatmap(ps.fam.rel, sample.label = "Treatment" , taxa.label="Family")

ps.fam.rel.nocu <- ps_nocu %>%
  aggregate_taxa(level = "Family") %>%  
  subset_taxa(Family != "Unknown_Family") %>%   
  subset_taxa(Family != "Unknown") %>%   
  microbiome::transform(transform = "compositional")

hm_fam_nocu <- plot_heatmap(ps.fam.rel.nocu
                            , taxa.label = "Family"
                            , sample.label = "Fertilization"
                            )

ps.fam.rel.yescu <- ps_yescu %>%
  aggregate_taxa(level = "Family") %>%  
  subset_taxa(Family != "Unknown") %>%   
  microbiome::transform(transform = "compositional")

hm_fam_yescu <- plot_heatmap(ps.fam.rel.nocu, taxa.label="Family")

####

ps.genus.rel.nocu <- ps_nocu %>%
  aggregate_taxa(level = "Genus") %>%  
  subset_taxa(Genus != "Unknown") %>%   
  microbiome::transform(transform = "compositional")
  

hm_genus_nocu <- plot_heatmap(ps.genus.rel.nocu, method = "PCoA", distance = "bray", sample.label = "Fertilization", sample.order = "Fertilization" ,  taxa.label="Genus", taxa.order = "Phylum")

ps.genus.rel.yescu <- ps_yescu %>%
  aggregate_taxa(level = "Genus") %>%  
  subset_taxa(Genus != "Unknown") %>% 
  microbiome::transform(transform = "compositional", scale = 100)


hm_genus_yescu <- plot_heatmap(ps.genus.rel.yescu, sample.label = "Fertilization", sample.order = "Fertilization" , taxa.label="Genus", taxa.order = "Phylum")

ggarrange(hm_genus_nocu, hm_genus_yescu, ncol=2)

ggarrange(hm_phylum_nocu, hm_phylum_yescu, hm_genus_nocu, hm_genus_yescu, ncol=2, nrow = 2)

# gen.nocu.hm <- d3heatmap(abundances(ps.genus.rel.nocu), 
#           hclustfun = hclust,
#           #scale = "row", # generates different scales based on the branches of the rows
#           Rowv = T, 
#           Colv = T,
#           colors = "RdYlGn",
#           k_row = 2, # Numero de grupos en las filas
#           k_col = 2, # Numero de grupos en las columnas
#           labCol = sample_data(ps.genus.rel.nocu)$Fertilization, 
#           digits = 10L
# )
# 
# d3heatmap(abundances(ps.genus.rel.yescu),
#           hclustfun = hclust,
#           #scale = "row", # generates different scales based on the branches of the rows
#           Rowv = T,
#           Colv = T,
#           #reorderfun = ,
#           colors = "YlOrRd",
#           k_row = 2, # Numero de grupos en las filas
#           k_col = 2, # Numero de grupos en las columnas
#           labCol = sample_data(ps.genus.rel.yescu)$Fertilization,
#           digits = 10L
# )
# 
# 
# saveWidget(gen.nocu.hm, "test.html")
# 
# d3heatmap(abundances(ps.fam.rel.nocu), 
#                          hclustfun = hclust,
#                          #scale = "row", # generates different scales based on the branches of the rows
#                          Rowv = T, 
#                          Colv = T,
#                          colors = "BuGn",
#                          k_row = 2, # Numero de grupos en las filas
#                          k_col = 2, # Numero de grupos en las columnas
#                          labCol = sample_data(ps.fam.rel.nocu)$Fertilization, 
#                          digits = 10L
# )

# 4.4.3) Taxonomic PCA

### Create new function from ggbiplot

ggbiplot2 <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                      obs.scale = 1 - scale, var.scale = scale, 
                      groups = NULL, ellipse = FALSE, ellipse.prob = 0.68, 
                      labels = NULL, labels.size = 3, alpha = 1, 
                      var.axes = TRUE, 
                      circle = FALSE, circle.prob = 0.69, 
                      varname.size = 3, varname.adjust = 1.5, 
                      varname.abbrev = FALSE, 
                      arrow.color = muted("red"), 
                      arrow.linetype = "solid",
                      arrow.alpha = 1, 
                      arrow.size = 1,
                      varname.angle = 0,
                      varname.color = 'darkred',
                      varname.border = NA, 
                      varname.fill = NA,
                      varname.ff = "plain",
                      ellipse.linetype = "solid",
                      ellipse.linewidth = 0.3,
                      ellipse.fill = 0,
                      intercept = FALSE)
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  
  stopifnot(length(choices) == 2)
  
  # Recover the SVD
  if(inherits(pcobj, 'prcomp')){
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$rotation
  } else if(inherits(pcobj, 'princomp')) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$loadings
  } else if(inherits(pcobj, 'PCA')) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN="/")
  } else if(inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  } else {
    stop('Expected a object of class prcomp, princomp, PCA, or lda')
  }
  
  # Scores
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))
  
  # Directions
  v <- sweep(v, 2, d^var.scale, FUN='*')
  df.v <- as.data.frame(v[, choices])
  
  names(df.u) <- c('xvar', 'yvar')
  names(df.v) <- names(df.u)
  
  if(pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  
  # Scale the radius of the correlation circle so that it corresponds to 
  # a data ellipse for the standardized PC scores
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  
  # Scale directions
  v.scale <- rowSums(v^2)
  df.v <- r * df.v / sqrt(max(v.scale))
  
  # Change the labels for the axes
  if(obs.scale == 0) {
    u.axis.labs <- paste('standardized PC', choices, sep='')
  } else {
    u.axis.labs <- paste('PC', choices, sep='')
  }
  
  # Append the proportion of explained variance to the axis labels
  u.axis.labs <- paste(u.axis.labs, 
                       sprintf('(%0.1f%% explained var.)', 
                               100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  
  # Score Labels
  if(!is.null(labels)) {
    df.u$labels <- labels
  }
  
  # Grouping variable
  if(!is.null(groups)) {
    df.u$groups <- groups
  }
  
  # Variable Names
  if(varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else {
    df.v$varname <- rownames(v)
  }
  
  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)
  
  # Base plot
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + 
    xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) + coord_equal()
  
  if(var.axes) {
    # Draw circle
    if(circle) 
    {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
      g <- g + geom_path(data = circle, color = muted('white'), 
                         size = 1/2, alpha = 1/3)
    }
    
    # Draw directions
    g <- g +
      geom_segment(data = df.v,
                   aes(x = 0, y = 0, xend = xvar*arrow.size, yend = yvar*arrow.size),
                   arrow = arrow(length = unit(1/2, 'picas')), 
                   color = arrow.color, linetype = arrow.linetype, alpha = arrow.alpha)
  }
  
  # Draw either labels or points
  if(!is.null(df.u$labels)) {
    if(!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    } else {
      g <- g + geom_text(aes(label = labels), size = labels.size)      
    }
  } else {
    if(!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    } else {
      g <- g + geom_point(alpha = alpha)      
    }
  }
  
  # Overlay a concentration ellipse if there are groups
  if(!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    
    ell <- ddply(df.u, 'groups', function(x) {
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'), 
                 groups = x$groups[1])
    })
    names(ell)[1:2] <- c('xvar', 'yvar')
    g <- g + geom_path(data = ell, aes(color = groups, group = groups),  
                       linetype = ellipse.linetype, linewidth = ellipse.linewidth)
  }
  
  # Label the variable axes
  if(var.axes) {
    g <- g + 
      geom_label(data = df.v, 
                 aes(label = varname, x = xvar*arrow.size, y = yvar*arrow.size, 
                     angle = varname.angle, hjust = hjust), 
                 color = varname.color, size = varname.size, label.size = varname.border, fill = varname.fill, fontface = varname.ff)
  }
  # Draw 0,0 intercept lines
  if(intercept){
    g <- g + 
      geom_hline(yintercept=0, linetype="dashed") + 
      geom_vline(xintercept=0, linetype="dashed")
  }
  # Change the name of the legend for groups
  # if(!is.null(groups)) {
  #   g <- g + scale_color_brewer(name = deparse(substitute(groups)), 
  #                               palette = 'Dark2')
  # }
  
  # TODO: Add a second set of axes
  
  return(g)
}

ps.phylum <- ps_rar %>%
  aggregate_taxa(level = "Phylum")

res.pca.phylum <- prcomp(x = t(abundances(ps.phylum)), center = TRUE, scale. = TRUE)

ps.phylum.nocu <- ps_nocu %>%
  aggregate_taxa(level = "Phylum") 
ps.phylum.nocu <- prune_taxa(taxa_sums(ps.phylum.nocu) > 0, ps.phylum.nocu)

res.pca.phylum.nocu <- prcomp(x = t(abundances(ps.phylum.nocu)), center = TRUE, scale. = TRUE)

ps.phylum.yescu <- ps_yescu %>%
  aggregate_taxa(level = "Phylum")
ps.phylum.yescu <- prune_taxa(taxa_sums(ps.phylum.yescu) > 0, ps.phylum.yescu)

res.pca.phylum.yescu <- prcomp(x = t(abundances(ps.phylum.yescu)), center = TRUE, scale. = TRUE)

# fviz_eig(res.pca.phylum)

fviz_pca_ind(res.pca.phylum, # Grafica las observaciones en los ejes PC1 y PC2.
             col.ind = "cos2", # Color según calidad de representación
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # pseudocolor de cálido a frío
             repel = FALSE     # No evitar sobrelape de texto
)

fviz_pca_var(res.pca.phylum, # Grafica las variables del df en los ejes PC1 y PC2.
             col.var = "contrib", # Color según contribución a los PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # pseudocolor de cálido a frío
             repel = TRUE     # Evitar sobrelape de texto
)

fviz_pca_biplot(res.pca.phylum, repel = TRUE, # Representa en un mismo gráfico las observaciones y las variables
                col.var = "#2E9FDF", # Un color para las variables
                geom.ind = "point", # Dibujar observaciones sin colocar etiquetas
                col.ind = "#696969"  # Un color para las observaciones
)


pca_biplot_nocu <- ggbiplot2(res.pca.phylum.nocu, 
          groups = meta(ps.phylum.nocu)$Fertilization, # Procesamiento en función de factor Treatment en dataframe, que incluye condicoines y tratamientos
          ellipse = TRUE, # Crear elipses
          ellipse.prob = 0.70, # Intervalo de confianza para la elipse
          ellipse.linewidth = 0.5, # Ancho de línea para la elipse
          obs.scale = 1, # Ajuste escala de puntos
          var.scale = 1, # De alguna forma modifica la escala del eje X. Ni pico idea cómo funciona
          arrow.size = 1, # Escala de las flechas de variables para facilitar visualización
          arrow.color = "#028202", # Color de flechas
          arrow.alpha = 1, # opacidad de flechas
          varname.color = "black", # Color de etiquetas de variables
          varname.adjust = 1.1, # Aleja un poco las etiquetas de las flechas
          varname.size = 3, # Tamaño de las etiquetas de las variables
          varname.ff = "bold", # etiquetas en negrita
          intercept = TRUE
) + theme_bw() +
  # stat_conf_ellipse(level = 0.95) +
  labs(title = "Sin cobre") +
  scale_color_manual(values=c("#4C4C4C", 
                              "#0F80FF", 
                              "#62C7FA", 
                              "#A1A3A8"),
                     labels= c('Control','Hoagland','AMH3-8','Suelo NR')) + # labels
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
        axis.text = element_text(face = "bold", size = 12), # Axis text format
        axis.title = element_text(face = "bold", size = 12), # Axis title format 
        plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
        legend.text = element_text(size = 12), # Legend text format
        legend.title = element_blank()) + # Remove legend title
  xlab("PC 1 (28,7%)") + # Título eje X
  ylab("PC 2 (22,7%)") # Título eje Y

pca_biplot_yescu <- ggbiplot2(res.pca.phylum.yescu, 
          groups = meta(ps.phylum.yescu)$Fertilization, # Procesamiento en función de factor Treatment en dataframe, que incluye condicoines y tratamientos
          ellipse = TRUE, # Crear elipses
          ellipse.prob = 0.70, # Intervalo de confianza para la elipse
          ellipse.linewidth = 0.5, # Ancho de línea para la elipse
          obs.scale = 1, # Ajuste escala de puntos
          var.scale = 1, # De alguna forma modifica la escala del eje X. Ni pico idea cómo funciona
          arrow.size = 1, # Escala de las flechas de variables para facilitar visualización
          arrow.color = "#DE7920", # Color de flechas
          arrow.alpha = 1, # opacidad de flechas
          varname.color = "black", # Color de etiquetas de variables
          varname.adjust = 1.1, # Aleja un poco las etiquetas de las flechas
          varname.size = 3, # Tamaño de las etiquetas de las variables
          varname.ff = "bold", # etiquetas en negrita
          intercept = TRUE
) + theme_bw() +
  labs(title = "Con cobre") +
  scale_color_manual(values=c("#4C4C4C", 
                              "#0F80FF", 
                              "#62C7FA", 
                              "#A1A3A8"),
                     labels= c('Control','Hoagland','AMH3-8','Suelo NR')) + # labels
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
        axis.text = element_text(face = "bold", size = 12), # Axis text format
        axis.title = element_text(face = "bold", size = 12), # Axis title format 
        plot.title = element_text( face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
        legend.text = element_text(size = 12), # Legend text format
        legend.title = element_blank()) + # Remove legend title
  xlab("PC 1 (28,5%)") + # Título eje X
  ylab("PC 2 (21,8%)") # Título eje Y

ggarrange(pca_biplot_nocu,pca_biplot_yescu, nrow=1, ncol=2, legend = "bottom", common.legend = TRUE, labels = "AUTO", font.label = list(size = 18), align = "hv")

# 4.4.4) Core microbiome analysis

## 4.4.4.1) Phylum level

ps.phylum.rel.r.nocu <- ps_r_nocu %>%
  aggregate_taxa(level = "Phylum") %>%  
  microbiome::transform(transform = "compositional")

ps.phylum.rel.r.yescu <- ps_r_yescu %>%
  aggregate_taxa(level = "Phylum") %>%  
  microbiome::transform(transform = "compositional")

ps.phylum.rel.r <- ps_rar_r %>%
  aggregate_taxa(level = "Phylum") %>%  
  microbiome::transform(transform = "compositional")

coreplot_nocu <- plot_core(ps.phylum.rel.r.nocu,
          plot.type = "heatmap",
          colours=green_pal,
          prevalences = seq(0.1, 1, 0.2), 
          detections = round(2^(-10:-1),3)) + 
  labs(title = "Sin cobre") +
  theme(
        axis.text = element_text(size = 8), # Axis text format
        axis.title = element_text(face = "bold", size = 10), # Axis title format 
        plot.title = element_text( face = "bold", size = rel(1.6), margin = margin(b = 5), hjust = 0.5), # Plot title format
        legend.text = element_text(size = 12), # Legend text format
        legend.title = element_blank()) +# Remove legend title
  xlab("Límite de detección")
  # ylab("Phylum")

coreplot_yescu <- plot_core(ps.phylum.rel.r.yescu,
                           plot.type = "heatmap",
                           colours=orange_pal,
                           prevalences = seq(0.1, 1, 0.2), 
                           detections = round(2^(-10:-1),3)) + 
  labs(title = "Con cobre") +
  theme(
    axis.text = element_text(size = 8), # Axis text format
    axis.title = element_text(face = "bold", size = 10), # Axis title format 
    plot.title = element_text( face = "bold", size = rel(1.6), margin = margin(b = 5), hjust = 0.5), # Plot title format
    legend.text = element_text(size = 12), # Legend text format
    legend.title = element_blank()) + # Remove legend title
  xlab("Límite de detección")

ggarrange(coreplot_nocu,coreplot_yescu, nrow=1, ncol=2, labels = "AUTO", font.label = list(size = 18), align = "hv")

coreplot_phy <- plot_core(ps.phylum.rel.r,
                           plot.type = "heatmap",
                           colours=brewer.pal(9, "BuGn"),
                           prevalences = seq(0.1, 1, 0.2), 
                           detections = round(2^(-10:-1),3)) + 
  geom_hline(yintercept=11.5, linetype="dashed", color = "#DE7920", size = .7) + 
  geom_hline(yintercept=9.5, linetype="dashed", color = "#DE7920", size = .7) + 
  geom_vline(xintercept=4.5, linetype="dashed", color = "#DE7920", size = .7) + 
  geom_hline(yintercept=6.5, linetype="dashed", color = "#DE7920", size = .7) +
  geom_vline(xintercept=1.5, linetype="dashed", color = "#DE7920", size = .7) +
  # geom_vline(xintercept=2.5, linetype="dashed", color = "blue") + 
  theme(
    axis.text = element_text(size = 8), # Axis text format
    axis.title = element_text(face = "bold", size = 10), # Axis title format 
    plot.title = element_text( face = "bold", size = rel(1.6), margin = margin(b = 5), hjust = 0.5), # Plot title format
    legend.text = element_text(size = 12), # Legend text format
    legend.title = element_blank()) +# Remove legend title
  xlab("Detection limit")
  # ylab("Phylum")

## 4.4.4.2) Genus level

ps.genus.rel.r.nocu <- ps_r_nocu %>%
  aggregate_taxa(level = "Genus") %>%  
  subset_taxa(Genus != "Unknown") %>%   
  microbiome::transform(transform = "compositional")

ps.genus.rel.r.yescu <- ps_r_yescu %>%
  aggregate_taxa(level = "Genus") %>%  
  subset_taxa(Genus != "Unknown") %>%   
  microbiome::transform(transform = "compositional")

ps.genus.rel.r <- ps_rar_r %>%
  aggregate_taxa(level = "Genus") %>%  
  subset_taxa(Genus != "Unknown") %>%   
  microbiome::transform(transform = "compositional")

coreplot_genus_nocu <- plot_core(ps.genus.rel.r.nocu,
                           plot.type = "heatmap",
                           colours=green_pal,
                           prevalences = seq(0.1, 1, 0.2), 
                           detections = round(2^(-10:-1),3)) + 
  # labs(title = "Sin cobre") +
  theme(
    axis.text = element_text(size = 8), # Axis text format
    axis.title = element_text(face = "bold", size = 10), # Axis title format 
    plot.title = element_blank(), # Plot title format
    legend.text = element_text(size = 12), # Legend text format
    legend.title = element_blank()) +# Remove legend title
  xlab("Detection limit")
  # ylab("Genus")

coreplot_genus_yescu <- plot_core(ps.genus.rel.r.yescu,
                            plot.type = "heatmap",
                            colours=orange_pal,
                            prevalences = seq(0.1, 1, 0.2), 
                            detections = round(2^(-10:-1),3)) +
  # labs(title = "Con cobre") +
  theme(
    axis.text = element_text(size = 8), # Axis text format
    axis.title = element_text(face = "bold", size = 10), # Axis title format 
    plot.title = element_blank(), # Plot title format
    legend.text = element_text(size = 12), # Legend text format
    legend.title = element_blank()) +# Remove legend title
  xlab("Detection limit")
  
ggarrange(coreplot_genus_nocu,coreplot_genus_yescu, nrow=1, ncol=2, labels = "AUTO", font.label = list(size = 18), align = "hv")

coreplot_gen <- plot_core(ps.genus.rel.r,
                          plot.type = "heatmap",
                          colours=brewer.pal(9, "BuGn"),
                          prevalences = seq(0.1, 1, 0.2), 
                          detections = round(2^(-10:-1),3)) + 
  geom_hline(yintercept=137.5, linetype="dashed", color = "#DE7920", size = .7) + 
  geom_hline(yintercept=133.5, linetype="dashed", color = "#DE7920", size = .7) + 
  geom_vline(xintercept=4.5, linetype="dashed", color = "#DE7920", size = .7) + 
  geom_hline(yintercept=124.5, linetype="dashed", color = "#DE7920", size = .7) + 
  geom_vline(xintercept=1.5, linetype="dashed", color = "#DE7920", size = .7) +
  theme(
    axis.text = element_text(size = 8), # Axis text format
    axis.title = element_text(face = "bold", size = 10), # Axis title format 
    plot.title = element_text( face = "bold", size = rel(1.6), margin = margin(b = 5), hjust = 0.5), # Plot title format
    legend.text = element_text(size = 12), # Legend text format
    legend.title = element_blank()) +# Remove legend title
  xlab("Detection limit")
  # ylab("Género")


# 4.2.5) Differential abundances analysis

###MicrobiotaProcess
set.seed(1)

ps_rar_r_AC_nocu <- rarefy_even_depth(ps_r_AC_nocu, sample.size = 20000)
ps_rar_r_HC_nocu <- rarefy_even_depth(ps_r_HC_nocu, sample.size = 20000)
ps_rar_r_AC_yescu <- rarefy_even_depth(ps_r_AC_yescu, sample.size = 20000)
ps_rar_r_HC_yescu <- rarefy_even_depth(ps_r_HC_yescu, sample.size = 20000)

taxa_table <- phyloseq::tax_table(ps_rar_r_AC_nocu)@.Data

biom <- diff_analysis(obj = ps_rar_r, 
                              # sampleda = sample_data(ps_rar_r_AC_nocu),
                              classgroup = "Fertilization",
                              # taxda = taxa_table,
                              # alltax = T,
                              standard_method = "hellinger",
                              mlfun = "lda",
                              clmin = 2,
                              filtermod = "pvalue",
                              firstcomfun = "kruskal_test",
                              firstalpha = 0.05,
                              strictmod = TRUE,
                              secondcomfun = "wilcox.test",
                              normalization = 1e+07,
                              ldascore = 1,
                              subclmin = 2,
                              subclwilc = TRUE,
                              secondalpha = 0.05,
                              bootnums = 100
)

biomarker <- ggdiffclade(
  obj=biom, 
  alpha=0.3, 
  linewd=0.3,
  skpointsize=1.5, 
  layout="radial",
  taxlevel=6, # Show the full name of up until the 6th highest LDA value
  cladetext = 3.3,
  bg.tree.color = "black", 
  bg.point.color = "black",
  removeUnknown=TRUE,
  reduce=TRUE # This argument is to remove the branch of unknown taxonomy.
) +
  scale_fill_manual(
    values=tre_col, 
    labels=tr_lab,
    guide = "none"
  ) +
  theme(
    panel.background=element_blank(),
    panel.margin= unit(c(0, 0, 0, 0), "lines"), 
    legend.position="right", 
    plot.margin=margin(0,0,0,0),
    legend.spacing.y=unit(0.2, "cm"), 
    legend.title=element_text(size=11),
    legend.text=element_text(size=10), 
    legend.box.spacing=unit(0.1,"cm")
  )

biomarker_phy <-  biomarker  +
  geom_hilight(
    data = td_filter(nodeClass == "Phylum"),
    mapping = aes(node = node, fill = label),
    alpha = 0.3
  ) + 
  scale_fill_manual(values = c(phy_col2,tre_col), labels = phy_lab2, guide = "none")

###

diff_ac_nocu <- diff_analysis(obj = ps_rar_r_AC_nocu, 
                      # sampleda = sample_data(ps_rar_r_AC_nocu),
                      classgroup = "Fertilization",
                      # taxda = taxa_table,
                      # alltax = T,
                      standard_method = "hellinger",
                       mlfun = "lda",
                       clmin = 2,
                       filtermod = "pvalue",
                       firstcomfun = "glm",
                       firstalpha = 0.05,
                       strictmod = TRUE,
                       secondcomfun = "wilcox_test",
                      normalization = 1e+06,
                      ldascore = 3,
                       subclmin = 2,
                       subclwilc = TRUE,
                       secondalpha = 0.05,
                      bootnums = 100
                      )

diff_hc_nocu <- diff_analysis(obj = ps_r_HC_nocu, 
                              # sampleda = sample_data(ps_nocu), 
                              classgroup = "Fertilization",
                              mlfun = "lda",
                              clmin = 2,
                              filtermod = "pvalue",
                              firstcomfun = "glm",
                              firstalpha = 0.05,
                              strictmod = TRUE,
                              secondcomfun = "wilcox_test",
                              normalization = 1e+06,
                              ldascore = 3,
                              subclmin = 2,
                              subclwilc = TRUE,
                              secondalpha = 0.05,
                              bootnums = 100)

diffplot_ac_nocu <- ggdiffbox(obj=diff_ac_nocu, box_notch=FALSE, 
                     colorlist=c("#4C4C4C", 
                                 "#62C7FA"), l_xlabtext="Abundancia relativa (%)")
diffplot_hc_nocu <- ggdiffbox(obj=diff_hc_nocu, box_notch=FALSE, 
                              colorlist=c("#4C4C4C", 
                                          "#0F80FF"), l_xlabtext="Abundancia relativa (%)")
ggarrange(diffplot_ac_nocu, diffplot_hc_nocu, ncol = 2)
# View(phyloseq::tax_table(ps_diff_ac_nc))
diff_ac_yescu <- diff_analysis(obj = ps_rar_r_AC_yescu, 
                               # sampleda = sample_data(ps_r_AC_yescu), 
                               classgroup = "Fertilization",
                               mlfun = "lda",
                               clmin = 2,
                               ratio = 1,
                               filtermod = "pvalue",
                               firstcomfun = "glm",
                               firstalpha = 0.05,
                               strictmod = TRUE,
                               secondcomfun = "wilcox_test",
                               normalization = 1e+07,
                               ldascore = 3,
                               subclmin = 2,
                               subclwilc = TRUE,
                               secondalpha = 0.15,
                               bootnum = 100)

diff_hc_yescu <- diff_analysis(obj = ps_r_HC_yescu, 
                               # sampleda = sample_data(ps_r_HC_yescu), 
                               classgroup = "Fertilization",
                               mlfun = "lda",
                               clmin = 2,
                               ratio = 1,
                               filtermod = "pvalue",
                               firstcomfun = "glm",
                               firstalpha = 0.05,
                               strictmod = TRUE,
                               secondcomfun = "wilcox_test",
                               normalization = 1e+06,
                               ldascore = 3,
                               subclmin = 2,
                               subclwilc = TRUE,
                               secondalpha = 0.15,
                               bootnum = 100)

diffplot_ac_yescu <- ggdiffbox(obj=diff_ac_yescu, box_notch=FALSE, 
                              colorlist=c("#4C4C4C", 
                                          "#62C7FA"), l_xlabtext="Abundancia relativa (%)")
diffplot_hc_yescu <- ggdiffbox(obj=diff_hc_yescu, box_notch=FALSE, 
                              colorlist=c("#4C4C4C", 
                                          "#0F80FF"), l_xlabtext="Abundancia relativa (%)")
ggarrange(diffplot_ac_yescu, diffplot_hc_yescu, ncol = 2)

ggarrange(diffplot_ac_nocu, diffplot_ac_yescu, diffplot_hc_nocu, diffplot_hc_yescu, nrow = 2, ncol = 2)

# Correlation taxonomy


library(ggnewscale)
library(WGCNA)
library(reshape2)

biomark_ac_nocu <- ggdiffclade(
  obj=diff_ac_nocu,
  alpha=0.3, 
  linewd=0.3,
  skpointsize=1.5, 
  layout="inward_circular",
  taxlevel=7,
  cladetext = 0,
  bg.tree.color = "black", 
  bg.point.color = "black",
  setColors=FALSE,
  xlim=16
) +
  scale_fill_manual(values=tre_col, 
                    labels=tr_lab,
                    guide=guide_legend(keywidth=0.5,
                                       keyheight=0.5,
                                       order=3,
                                       override.aes=list(alpha=1))
  ) +
  scale_size_continuous(range=c(1, 3),
                        guide=guide_legend(keywidth=0.5,keyheight=0.5,order=4,
                                           override.aes=list(shape=21))) +
  scale_colour_manual(values=rep("white", 100),guide="none")

genustab <- get_taxadf(ps_rar_r_AC_nocu, taxlevel=7)
genustab <- data.frame(t(otu_table(genustab)), check.names=FALSE)
genustab <- data.frame(apply(genustab, 2, function(x)x/sum(x)), check.names=FALSE)

cortest <- WGCNA::corAndPvalue(genustab, method="spearman", alternative="two.sided")
cortest$cor[upper.tri(cortest$cor, diag = TRUE)] <- NA
cortest$p[upper.tri(cortest$p, diag = TRUE)] <- NA
cortab1 <- na.omit(melt(t(cortest$cor))) %>% rename(from=Var1,to=Var2,cor=value)
corptab1 <- na.omit(melt(t(cortest$p))) %>% rename(pvalue=value)
cortab1$fdr <- p.adjust(corptab1$pvalue, method="fdr")

cortab1 <- cortab1 %>% mutate(correlation=case_when(cor>0 ~ "positive",cor < 0 ~ "negative",TRUE ~ "No"))
cortab2 <- cortab1 %>% filter(fdr <= 0.05) %>% filter(cor <= -0.5 | cor >= 0.8)
cortab2gem <- cortab2[grepl("Gemma", cortab2$from) | grepl("Gemma", cortab2$to), ]

biom_corr_ac_nocu <- biomark_ac_nocu +
  new_scale_color() +
  new_scale("size") +
  geom_tiplab(size=1.5, hjust=1) +
  geom_taxalink(
    data=cortab2gem,
    mapping=aes(taxa1=from,
                taxa2=to,
                colour=correlation,
                size=abs(cor)),
    alpha=0.4,
    ncp=10,
    hratio=1,
    offset=0.31
  ) +
  scale_size_continuous(range = c(0.2, 1),
                        guide=guide_legend(keywidth=1, keyheight=1,
                                           order=1, override.aes=list(alpha=1))
  ) +
  scale_colour_manual(values=c("negative" = "chocolate2", "positive" = "#009E73"),
                      guide=guide_legend(keywidth=1, keyheight=1,
                                         order=2, override.aes=list(alpha=1, size=1))) + 
  theme(legend.text = element_text(size = 12)
        , legend.title = element_text(size = 14)
  )

biomark_ac_yescu <- ggdiffclade(
  obj=diff_ac_yescu,
  alpha=0.3, 
  linewd=0.3,
  skpointsize=1.5, 
  layout="inward_circular",
  taxlevel=7,
  cladetext = 0,
  bg.tree.color = "black", 
  bg.point.color = "black",
  setColors=FALSE,
  xlim=16
) +
  scale_fill_manual(values=tre_col, 
                    labels=tr_lab,
                    guide=guide_legend(keywidth=0.5,
                                       keyheight=0.5,
                                       order=3,
                                       override.aes=list(alpha=1))
  ) +
  scale_size_continuous(range=c(1, 3),
                        guide=guide_legend(keywidth=0.5,keyheight=0.5,order=4,
                                           override.aes=list(shape=21))) +
  scale_colour_manual(values=rep("white", 100),guide="none")

genustab2 <- get_taxadf(ps_rar_r_AC_yescu, taxlevel=7)
genustab2 <- data.frame(t(otu_table(genustab2)), check.names=FALSE)
genustab2 <- data.frame(apply(genustab2, 2, function(x)x/sum(x)), check.names=FALSE)

cortest2 <- WGCNA::corAndPvalue(genustab2, method="spearman", alternative="two.sided")
cortest2$cor[upper.tri(cortest2$cor, diag = TRUE)] <- NA
cortest2$p[upper.tri(cortest2$p, diag = TRUE)] <- NA
cortab3 <- na.omit(melt(t(cortest2$cor))) %>% rename(from=Var1,to=Var2,cor=value)
corptab2 <- na.omit(melt(t(cortest2$p))) %>% rename(pvalue=value)
cortab3$fdr <- p.adjust(corptab2$pvalue, method="fdr")

cortab3 <- cortab3 %>% mutate(correlation=case_when(cor>0 ~ "positive",cor < 0 ~ "negative",TRUE ~ "No"))
cortab4 <- cortab3 %>% filter(fdr <= 0.05) %>% filter(cor <= -0.5 | cor >= 0.8)
cortab4gem <- cortab4[grepl("Gemma", cortab4$from) | grepl("Gemma", cortab4$to), ]

biom_corr_ac_yescu <- biomark_ac_yescu +
  new_scale_color() +
  new_scale("size") +
  geom_tiplab(size=1.5, hjust=1) +
  geom_taxalink(
    data=cortab4gem,
    mapping=aes(taxa1=from,
                taxa2=to,
                colour=correlation,
                size=abs(cor)),
    alpha=0.4,
    ncp=10,
    hratio=1,
    offset=0.31
  ) +
  scale_size_continuous(range = c(0.2, 1),
                        guide=guide_legend(keywidth=1, keyheight=1,
                                           order=1, override.aes=list(alpha=1))
  ) +
  scale_colour_manual(values=c("negative" = "chocolate2", "positive" = "#009E73"),
                      guide=guide_legend(keywidth=1, keyheight=1,
                                         order=2, override.aes=list(alpha=1, size=1))) + 
  theme(legend.text = element_text(size = 12)
        , legend.title = element_text(size = 14)
        )

ggarrange(biom_corr_ac_nocu, biom_corr_ac_yescu, common.legend = T, legend = "bottom", ncol = 2)

###

# sample_data(ps.genus.nocu)$Fertilization
# otu_table(ps_rar)
# sam_data(ps_rar)
# tax_table(ps_rar)
# write.csv(t(otu_table(ps_rar)), "otu_table.csv", row.names=TRUE)
# write.csv(meta(ps_rar), "metadata_exp.csv", row.names=TRUE)
# write.csv(tax_table(ps_rar), "tax_table.csv", row.names=TRUE)
# 
# library("animalcules")
# run_animalcules()
# 
