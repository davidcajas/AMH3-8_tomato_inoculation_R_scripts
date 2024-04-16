## 1) Needed packages

# Graph tools
library("ggplot2")
library("RColorBrewer")
library("ggpubr")
library("dplyr")
library("ggalluvial") # Plotting flowbart plots
library("ggh4x") # extension tools for ggplot2
library("ggbiplot") # Biploting
library("FactoMineR") # PCA fuction
library("factoextra") # PCA visualization tools
# library("d3heatmap") # Plot interactive D3 heatmaps
library("htmlwidgets") # Export D3 plots to standalone HTML file
library("VennDiagram") # Venn diagram plotting tool
library("grid") # Grob graphics tool
library("gridExtra") # Grob extension tools
library(ggtree) # Plot taxa trees
library(ggtreeExtra) # Plot taxa trees
library(tidytree) # Tools for taxa tree objects
library(ggstar)
library(forcats)
library("ggside")
library(ggprism)  # Prism-like plots

# Microbiome-specific tools
library("ShortRead") # Manage and explore DNA libraries (step 3)
library("dada2") # Infere sample from sequences (steps 4-6)
library("microbiome") #
library("phyloseq") #
library("DECIPHER") # Alternative algorithm to assign taxonomy to ASV table (step 6) and necessary for philogenetic tree (step 8) **** Check if this library does not conflict with steps 1-5
library("phangorn") # Philogenetic tree algorithms
library("ape") # Necessary to use some phagorn functions
library("seqinr") # Write fasta files
library("devtools")
library("MicrobiotaProcess") # Microbiome analysis toolbox
library("ranacapa") # Necessary for rarefraction of data
library("ade4") # 
library("vegan") # PERMANOVA
library("corrr") # Beta distance calculation
library("coin") # Kruskal-Wallis test
library("ggplotify") # as.ggplot function to convert aplot to ggplot graphs

## 2) Importing the phyloseq file

setwd("/Users/davidcajasmunoz/Library/CloudStorage/GoogleDrive-dadavid.cajas@gmail.com/Mi unidad/Academia/Postgrado/PUCV/Tesis/Resultados/Biofertilización/2/5-ADN/R/Seq")

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

# 3.2) Add read count to metadata
readcounts <- readcount(ps_bac)
sample_data(ps_bac)$Readcount_no_rarefaction <- readcounts
#View(meta(ps_bac2))
print("Metadata now shows Read count")

# 3.3) Add sample names to metadata
sample_data(ps_bac)$Sample_names<-as.factor(sample_names(ps_bac))
print("Metadata now shows sample names")

# 3.4) Change column types of variables to factor

sample_variables(ps_bac)
str(meta(ps_bac))

sample_data(ps_bac)$Sample_names<-as.factor(sample_names(ps_bac))
sample_data(ps_bac)$Treatment<-as.factor(sample_data(ps_bac)$Treatment)
sample_data(ps_bac)$Replica<-as.factor(sample_data(ps_bac)$Replica)
sample_data(ps_bac)$Fertilization<-as.factor(sample_data(ps_bac)$Fertilization)
sample_data(ps_bac)$Fertilization <- factor(sample_data(ps_bac)$Fertilization, levels = c("No","Hoagland","AMH3-8","No_plant")) #Reordering levels in Fertilization variable

str(meta(ps_bac))
print("Metadata type fixed")

# 3.5) More cleaning and sample subsetting

# 3.5.1) Remove low readcount samples

ps_bac2 <- subset_samples(ps_bac, Readcount_no_rarefaction > 20000)

# 3.5.2) Subset by Soil copper content

ps_nocu <- subset_samples(ps_bac2, Copper == "FALSE")
ps_nocu
ps_yescu <- subset_samples(ps_bac2, Copper == "TRUE")
ps_yescu

# 3.2) Eliminate extremely low abundance (< 100 counts) ASVs

ps_bac2 <-  prune_taxa(taxa_sums(ps_bac2) > 100, ps_bac2)

ps_nocu <- prune_taxa(taxa_sums(ps_nocu) > 100, ps_nocu)

ps_yescu <-  prune_taxa(taxa_sums(ps_yescu) > 100, ps_yescu)

# 3.5.2) Filter out bulk soil samples

#For some reason, filtering in ps_nocu the bulk soil samples gives the following error:
# "Component sample names do not match.
# Try sample_names()"
# An alternative workaround was implemented instead by subsampling out each individual sample. Since the metadata was lost in the process, it was extracted from the ps_nocu object, filtered and reinserted on ps_r_nocu object.

ps_r <- subset_samples(ps_bac2, Fertilization != "No_plant")

ps_r_nocu <- subset_samples(ps_nocu, sample_names(ps_nocu) != "12") %>% 
  subset_samples(sample_names(ps_nocu) != "11") %>% 
  subset_samples(sample_names(ps_nocu) != "10")

ps_r_nocu <- prune_taxa(taxa_sums(ps_r_nocu) > 0, ps_r_nocu) # Delete taxa with zero count



metadatanocu <- as.data.frame(sample_data(ps_nocu))
metadatanocu2 <- metadatanocu[metadatanocu$Fertilization!="No_plant",]
rownames(metadatanocu2) <- seq(1,9,1)
sample_data(ps_r_nocu) <-metadatanocu2

ps_r_yescu <- subset_samples(ps_yescu, Fertilization != "No_plant")

ps_r_yescu <- prune_taxa(taxa_sums(ps_r_yescu) > 0, ps_r_yescu) # Delete taxa with zero count

ps_bs <- subset_samples(ps_bac2, Fertilization == "No_plant")
ps_bs <- prune_taxa(taxa_sums(ps_bs) > 0, ps_bs)  # Delete taxa with zero count

# 3.5.2) Paired Phyloseq files for comparison between fertilization treatments

ps_r_AC_nocu <- subset_samples(ps_r_nocu, Fertilization != "Hoagland")
ps_r_AC_nocu <- prune_taxa(taxa_sums(ps_r_AC_nocu) > 0, ps_r_AC_nocu) # Delete taxa with zero count

ps_r_HC_nocu <- subset_samples(ps_r_nocu, sample_names(ps_r_nocu) != "8") %>% 
  subset_samples(sample_names(ps_r_nocu) != "7") %>% 
  subset_samples(sample_names(ps_r_nocu) != "7")
ps_r_HC_nocu <- prune_taxa(taxa_sums(ps_r_HC_nocu) > 0, ps_r_HC_nocu) # Delete taxa with zero count

metadatanocu3 <- as.data.frame(sample_data(ps_r_nocu))
metadatanocu4 <- metadatanocu3[metadatanocu3$Fertilization!="AMH3-8",]
rownames(metadatanocu4) <- seq(1,6,1)
sample_data(ps_r_HC_nocu) <-metadatanocu4

ps_r_AC_yescu <- subset_samples(ps_r_yescu, Fertilization != "Hoagland")
ps_r_AC_yescu <- prune_taxa(taxa_sums(ps_r_AC_yescu) > 0, ps_r_AC_yescu) # Delete taxa with zero count
ps_r_HC_yescu <- subset_samples(ps_r_yescu, Fertilization != "AMH3-8")
ps_r_HC_yescu <- prune_taxa(taxa_sums(ps_r_HC_yescu) > 0, ps_r_HC_yescu) # Delete taxa with zero count

# 3.6) Change to MicrobiotaProcess object

mp <- as.mpse(ps_bac2)

mp.nc <- as.mpse(ps_nocu)

mp.r.nc <- as.mpse(ps_r_nocu)

mp.r.ac.nc <- as.mpse(ps_r_AC_nocu)

mp.r.hc.nc <- as.mpse(ps_r_HC_nocu)

mp.yc <- as.mpse(ps_yescu)

mp.r.yc <- as.mpse(ps_r_yescu)

mp.r.ac.yc <- as.mpse(ps_r_AC_yescu)

mp.r.hc.yc <- as.mpse(ps_r_HC_yescu)

mp.r <- as.mpse(ps_r)

mp.bs <- as.mpse(ps_bs)

# 3.7) Auxiliary objects and functions

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

# 3.8) Rarefaction curve analysis

# 3.8.1) analysis
set.seed(1)
mp_rar <- mp_rrarefy(mp) %>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance
    , chunks = 400
  )

mp_rar.nc <- mp_rrarefy(mp.nc) %>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance
    , chunks = 400
  )

mp_rar.r.nc <- mp_rrarefy(mp.r.nc) %>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance
    , chunks = 400
  )

mp_rar.r.ac.nc <- mp_rrarefy(mp.r.ac.nc) %>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance
    , chunks = 400
  )

mp_rar.r.hc.nc <- mp_rrarefy(mp.r.hc.nc) %>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance
    , chunks = 400
  )

mp_rar.yc <- mp_rrarefy(mp.yc) %>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance
    , chunks = 400
  )

mp_rar.r.yc <- mp_rrarefy(mp.r.yc) %>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance
    , chunks = 400
  )

mp_rar.r.ac.yc <- mp_rrarefy(mp.r.ac.yc) %>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance
    , chunks = 400
  )

mp_rar.r.hc.yc <- mp_rrarefy(mp.r.hc.yc) %>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance
    , chunks = 400
  )

mp_rar.r <- mp_rrarefy(mp.r) %>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance
    , chunks = 400
  )

mp_rar.bs <- mp_rrarefy(mp.bs) %>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance
    , chunks = 400
  )

# 3.8.2) plot

mp_plot_rarecurve(mp_rar, 
                  .rare = RareAbundanceRarecurve, 
                  .alpha = Observe
) /
mp_plot_rarecurve(mp_rar.nc, 
                           .rare = RareAbundanceRarecurve, 
                           .alpha = Observe
) +
mp_plot_rarecurve(mp_rar.yc, 
                           .rare = RareAbundanceRarecurve, 
                           .alpha = Observe
) +
mp_plot_rarecurve(mp_rar.r.nc, 
                     .rare = RareAbundanceRarecurve, 
                     .alpha = Observe
) +
  mp_plot_rarecurve(mp_rar.r.yc, 
                    .rare = RareAbundanceRarecurve, 
                    .alpha = Observe
)


# 4) Exploratory analysis 

# 4.1) Alpha diversity analysis

# 4.1.1) analysis

mp_a <- mp_cal_alpha(mp_rar
             , .abundance=RareAbundance
             )

mp_a.nc <- mp_cal_alpha(mp_rar.nc
                     , .abundance=RareAbundance
)

mp_a.r.nc <- mp_cal_alpha(mp_rar.r.nc
                        , .abundance=RareAbundance
)

mp_a.yc <- mp_cal_alpha(mp_rar.yc
                     , .abundance=RareAbundance
)


mp_a.r.yc <- mp_cal_alpha(mp_rar.r.yc
                        , .abundance=RareAbundance
)

mp_a.r <- mp_cal_alpha(mp_rar.r
                        , .abundance=RareAbundance
)

mp_a.bs <- mp_cal_alpha(mp_rar.bs
             , .abundance=RareAbundance
)

# 4.1.2) plot

mp_plot_alpha(mp_a
  , .group = Treatment
  , .alpha = c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
) /
mp_plot_alpha(mp_a.nc
              , .group = Fertilization
              , .alpha = c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
) /
mp_plot_alpha(mp_a.yc
              , .group = Fertilization
              , .alpha = c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
) /
mp_plot_alpha(mp_a.r.nc
                 , .group = Fertilization
                 , .alpha = c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
) /
  mp_plot_alpha(mp_a.r.yc
                , .group = Fertilization
                , .alpha = c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
) 



mp_plot_alpha(mp_a.nc
              , .group = Fertilization
              , .alpha = c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
) /
  mp_plot_alpha(mp_a.yc
                , .group = Fertilization
                , .alpha = c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
  )

# 4.2) Beta distance analysis

# 4.2.1) Standarization and distance calculation

unlist(distanceMethodList)
# Ordination analysis
# wunifrac , X r X , me , co, g, 19, rlb
# dpcoa , c , wb , r
dmet <- "r"

mp_ab <- mp_decostand(mp_a
                     , .abundance = Abundance
                     ) %>% 
  mp_cal_dist(.abundance = hellinger # Default standarization method is Hellinger
              , distmethod = dmet
  ) %>% 
  mp_adonis(.abundance=hellinger
            , .formula=~Fertilization
            , distmethod=dmet
            , permutations=9999
            , action="add") %>% 
  mp_cal_pcoa(.abundance=hellinger
              , distmethod=dmet)


mp_ab.nc <- mp_decostand(mp_a.nc
                     , .abundance = Abundance
                     ) %>% 
  mp_cal_dist(.abundance = hellinger # Default standarization method is Hellinger
              , distmethod = dmet
  ) %>% 
  mp_adonis(.abundance=hellinger
            , .formula=~Fertilization
            , distmethod=dmet
            , permutations=9999
            , action="add") %>% 
  mp_cal_pcoa(.abundance=hellinger
              , distmethod=dmet)

mp_ab.r.nc <- mp_decostand(mp_a.r.nc
                         , .abundance = Abundance
) %>% 
  mp_cal_dist(.abundance = hellinger # Default standarization method is Hellinger
              , distmethod = dmet
  ) %>% 
  mp_adonis(.abundance=hellinger
            , .formula=~Fertilization
            , distmethod=dmet
            , permutations=9999
            , action="add") %>% 
  mp_cal_pcoa(.abundance=hellinger
              , distmethod=dmet)

mp_ab.yc <- mp_decostand(mp_a.yc
                     , .abundance = Abundance
                     )%>% 
  mp_cal_dist(.abundance = hellinger # Default standarization method is Hellinger
              , distmethod = dmet
  ) %>% 
  mp_adonis(.abundance=hellinger
            , .formula=~Fertilization
            , distmethod=dmet
            , permutations=9999
            , action="add") %>% 
  mp_cal_pcoa(.abundance=hellinger
              , distmethod=dmet)

mp_ab.r.yc <- mp_decostand(mp_a.r.yc
                         , .abundance = Abundance
)%>% 
  mp_cal_dist(.abundance = hellinger # Default standarization method is Hellinger
              , distmethod = dmet
  ) %>% 
  mp_adonis(.abundance=hellinger
            , .formula=~Fertilization
            , distmethod=dmet
            , permutations=9999
            , action="add") %>% 
  mp_cal_pcoa(.abundance=hellinger
              , distmethod=dmet)

# 4.2.2) Ordination plot

# centroids <- aggregate(mp_ab.nc$`PCo1 (20.14%)`, by = list(mp_ab.nc$Fertilization), FUN = mean)
# centroids$y <- aggregate(mp_ab.nc$`PCo2 (18.46%)`, by = list(mp_ab.nc$Fertilization), FUN = mean)$x
# mp_extract_internal_attr(mp_ab.nc, name=adonis)
# mp_extract_dist(mp_ab.nc, distmethod = dmet)
# distance(ps_nocu, method=dmet, type="samples")
# TukeyHSD(aov(distance(ps_nocu, method=dmet, type="samples") ~ meta(ps_nocu)$Fertilization, data = meta(ps_nocu)))
# TukeyHSD(aov(mp_extract_dist(mp_ab.nc, distmethod = dmet) ~ mp_ab.nc@colData$Fertilization, data = as.data.frame(mp_ab.nc@colData)))

pcoa_all <- mp_plot_ord(mp_ab
            ,.ord = pcoa
            ,.group = Copper
            ,.color = Copper
            ,.size = Shannon
            # ,.alpha = 0.8
            # ,.alpha = Shannon
            ,ellipse = T
            ,show.legend = FALSE # don't display the legend of stat_ellipse
) +
  geom_point() + 
  # geom_point(data=centroids,
  #            mapping=aes(x=centroids$x, y=centroids$y)) +
  scale_fill_manual(values=c('#028202', 
                                    '#DE7920'),
                    labels= c('No Copper','Copper')) + # labels
  scale_color_manual(values=c('#028202', 
                                     '#DE7920'),
                     labels= c('No Copper','Copper')) + # labels (must be the same in both fill and color for ggplot to show them as a single legend)
  guides(size = F) + #Don't show legend for size
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
        axis.text = element_text(face = "bold", size = 12), # Axis text format
        axis.title = element_text(face = "bold", size = 12), # Axis title format 
        plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
        # legend.key.size = unit(0.5, "cm"), # Space between legend elements
        legend.text = element_text(size = 14), # Legend text format
        legend.title = element_blank(), # Remove legend title 
        legend.position = "bottom", 
        ggside.axis.text = element_blank(),# Remove side boxplots' labels 
        ggside.axis.ticks = element_blank()) # Remove side boxplots' ticks 

pcoa_nc <- mp_plot_ord(mp_ab.nc
                       ,.ord = pcoa
                       ,.group = Fertilization
                       ,.color = Fertilization
                       ,.size = Shannon
                       # ,.alpha = 0.8
                       # ,.alpha = Shannon
                       ,ellipse = F
                       ,show.legend = FALSE # don't display the legend of stat_ellipse
) +
  geom_point() + 
  # geom_point(data=centroids,
  #            mapping=aes(x=centroids$x, y=centroids$y)) +
  stat_conf_ellipse(level = 0.95) +
  scale_fill_manual(values=c("#4C4C4C", 
                             "#0F80FF", 
                             "#62C7FA", 
                             "#A1A3A8"),
                    labels= c('Control','Hoagland','AMH3-8','Bulk soil')) + # labels
  scale_color_manual(values=c("#4C4C4C", 
                              "#0F80FF", 
                              "#62C7FA", 
                              "#A1A3A8"),
                     labels= c('Control','Hoagland','AMH3-8','Bulk soil')) + # labels (must be the same in both fill and color for ggplot to show them as a single legend)
  ggtitle("No Copper") +
  guides(size = F) + #Don't show legend for size
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
      axis.text = element_text(face = "bold", size = 12), # Axis text format
      axis.title = element_text(face = "bold", size = 12), # Axis title format 
      plot.title = element_text(face = "bold", size = rel(1.8), margin = margin(b = 15), hjust = 0.5), # Plot title format
      # legend.key.size = unit(0.5, "cm"), # Space between legend elements
      legend.text = element_text(size = 12), # Legend text format
      legend.title = element_blank(), # Remove legend title 
      ggside.axis.text = element_blank(),# Remove side boxplots' labels 
      ggside.axis.ticks = element_blank()) # Remove side boxplots' ticks 
  

pcoa_yc <- mp_plot_ord(mp_ab.yc
                       ,.ord = pcoa
                       ,.group = Fertilization
                       ,.color = Fertilization
                       ,.size = Shannon
                       # ,.alpha = 0.8
                       # ,.alpha = Shannon
                       ,ellipse = F
                       ,show.legend = FALSE # don't display the legend of stat_ellipse
) +
  geom_point() + 
  # geom_point(data=centroids,
  #            mapping=aes(x=centroids$x, y=centroids$y)) +
  stat_conf_ellipse(level = 0.95) +
  scale_fill_manual(values=c("#4C4C4C", 
                             "#0F80FF", 
                             "#62C7FA", 
                             "#A1A3A8"),
                    labels= c('Control','Hoagland','AMH3-8','Bulk soil')) + # labels
  scale_color_manual(values=c("#4C4C4C", 
                              "#0F80FF", 
                              "#62C7FA", 
                              "#A1A3A8"),
                     labels= c('Control','Hoagland','AMH3-8','Bulk soil')) + # labels (must be the same in both fill and color for ggplot to show them as a single legend)
  ggtitle("Copper") +
  guides(size = F) + #Don't show legend for size
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
        axis.text = element_text(face = "bold", size = 12), # Axis text format
        axis.title = element_text(face = "bold", size = 12), # Axis title format 
        plot.title = element_text(face = "bold", size = rel(1.8), margin = margin(b = 15), hjust = 0.5), # Plot title format
        # legend.key.size = unit(0.5, "cm"), # Space between legend elements
        legend.text = element_text(size = 12), # Legend text format
        legend.title = element_blank(), # Remove legend title 
        ggside.axis.text = element_blank(),# Remove side boxplots' labels 
        ggside.axis.ticks = element_blank()) # Remove side boxplots' ticks 

pcoa.r_nc <- mp_plot_ord(mp_ab.r.nc
                       ,.ord = pcoa
                       ,.group = Fertilization
                       ,.color = Fertilization
                       ,.size = Shannon
                       # ,.alpha = 0.8
                       # ,.alpha = Shannon
                       ,ellipse = F
                       ,show.legend = FALSE # don't display the legend of stat_ellipse
) +
  geom_point() + 
  # geom_point(data=centroids,
  #            mapping=aes(x=centroids$x, y=centroids$y)) +
  stat_conf_ellipse(level = 0.95) +
  scale_fill_manual(values=c("#4C4C4C", 
                             "#0F80FF", 
                             "#62C7FA", 
                             "#A1A3A8"),
                    labels= c('Control','Hoagland','AMH3-8','Suelo NR')) + # labels
  scale_color_manual(values=c("#4C4C4C", 
                              "#0F80FF", 
                              "#62C7FA", 
                              "#A1A3A8"),
                     labels= c('Control','Hoagland','AMH3-8','Suelo NR')) + # labels (must be the same in both fill and color for ggplot to show them as a single legend)
  ggtitle("Sin cobre") +
  guides(size = F) + #Don't show legend for size
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
        axis.text = element_text(face = "bold", size = 12), # Axis text format
        axis.title = element_text(face = "bold", size = 12), # Axis title format 
        plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
        # legend.key.size = unit(0.5, "cm"), # Space between legend elements
        legend.text = element_text(size = 12), # Legend text format
        legend.title = element_blank(), # Remove legend title 
        ggside.axis.text = element_blank(),# Remove side boxplots' labels 
        ggside.axis.ticks = element_blank()) # Remove side boxplots' ticks 


pcoa.r_yc <- mp_plot_ord(mp_ab.r.yc
                       ,.ord = pcoa
                       ,.group = Fertilization
                       ,.color = Fertilization
                       ,.size = Shannon
                       # ,.alpha = 0.8
                       # ,.alpha = Shannon
                       ,ellipse = F
                       ,show.legend = FALSE # don't display the legend of stat_ellipse
) +
  geom_point() + 
  # geom_point(data=centroids,
  #            mapping=aes(x=centroids$x, y=centroids$y)) +
  stat_conf_ellipse(level = 0.95) +
  scale_fill_manual(values=c("#4C4C4C", 
                             "#0F80FF", 
                             "#62C7FA", 
                             "#A1A3A8"),
                    labels= c('Control','Hoagland','AMH3-8','Suelo NR')) + # labels
  scale_color_manual(values=c("#4C4C4C", 
                              "#0F80FF", 
                              "#62C7FA", 
                              "#A1A3A8"),
                     labels= c('Control','Hoagland','AMH3-8','Suelo NR')) + # labels (must be the same in both fill and color for ggplot to show them as a single legend)
  ggtitle("Con cobre") +
  guides(size = F) + #Don't show legend for size
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
        axis.text = element_text(face = "bold", size = 12), # Axis text format
        axis.title = element_text(face = "bold", size = 12), # Axis title format 
        plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
        # legend.key.size = unit(0.5, "cm"), # Space between legend elements
        legend.text = element_text(size = 12), # Legend text format
        legend.title = element_blank(), # Remove legend title 
        ggside.axis.text = element_blank(),# Remove side boxplots' labels 
        ggside.axis.ticks = element_blank()) # Remove side boxplots' ticks 

ggarrange(pcoa_nc,pcoa_yc, nrow=1, ncol=2, legend = "bottom", common.legend = TRUE, labels = "AUTO", font.label = list(size = 18))
ggarrange(pcoa.r_nc,pcoa.r_yc, nrow=1, ncol=2, legend = "bottom", common.legend = TRUE, labels = "AUTO", font.label = list(size = 18))

# 4.2.3) Distance plot

mp_plot_dist(mp_ab, .distmethod = r)

mp_plot_dist(mp_ab.nc
             , .distmethod = r
             , .group = Fertilization
             , group.test = T
) +
  mp_plot_dist(mp_ab.yc
               , .distmethod = r
               , .group = Fertilization
               , group.test = T
  ) +
  mp_plot_dist(mp_ab.r.nc
               , .distmethod = r
               , .group = Fertilization
               , group.test = T
  ) +
  mp_plot_dist(mp_ab.r.yc
               , .distmethod = r
               , .group = Fertilization
               , group.test = T
  )

# 4.2.3) Hierarchical clustering

mp_abc <- mp_cal_clust(mp_ab, 
    .abundance = hellinger, 
    distmethod = "wunifrac",
    hclustmethod = "average", # (UPGAE)
    action = "add" # action is used to control which result will be returned
  )


clustering_map <- ggtree(mp_extract_internal_attr(mp_abc, name = 'SampleClust')) + #Tree itself
  geom_tippoint(aes(color=Fertilization)) + # Tip points colored by Fertilization treatment
  scale_color_manual(values=tre_col,
                     labels= tr_lab) + # Scale for AMH3-8, Hoagland, No, No_plant (alphabetical)
  scale_x_continuous(expand=c(0, 0.01)) + # X axis scale
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.position = "bottom") + 
  geom_strip('22', '21', barsize=1.2, color='#DE7920', 
             label="Copper", angle = 270, offset=.002, offset.text=.005, extend = 0.3, fontsize = 5, hjust = 0.5) + # Strip to mark Copper soils
  geom_strip('6', '5', barsize=1.2, color='#028202', 
             label="No Copper", angle = 270, offset=.002, offset.text=.005, extend = 0.3, fontsize = 5, hjust = 0.5) # Strip to mark No Copper soils
clustering_map

# 4.3) Abundance analysis

# 4.3.1) Calculate relative abundance by sample and by group

mp_aba <- mp_cal_abundance(mp_ab # for each samples
  , .abundance = RareAbundance
  ) %>% 
  mp_cal_abundance( # for each groups 
      .abundance = RareAbundance
    , .group = Treatment
  )

mp_aba.nc <- mp_cal_abundance(mp_ab.nc # for each samples
                           , .abundance = RareAbundance
                           ) %>% 
  mp_cal_abundance( # for each groups 
    .abundance = RareAbundance
    , .group = Fertilization
  )

mp_aba.r.nc <- mp_cal_abundance(mp_ab.r.nc # for each samples
                              , .abundance = RareAbundance
) %>% 
  mp_cal_abundance( # for each groups 
    .abundance = RareAbundance
    , .group = Fertilization
  )

mp_a.r.ac.nc <- mp_cal_abundance(mp_rar.r.ac.nc # for each samples
                                , .abundance = RareAbundance
) %>% 
  mp_cal_abundance( # for each groups 
    .abundance = RareAbundance
    , .group = Fertilization
  )

mp_a.r.hc.nc <- mp_cal_abundance(mp_rar.r.hc.nc # for each samples
                                 , .abundance = RareAbundance
) %>% 
  mp_cal_abundance( # for each groups 
    .abundance = RareAbundance
    , .group = Fertilization
  )

mp_aba.yc <- mp_cal_abundance(mp_ab.yc # for each samples
                           , .abundance = RareAbundance
                           ) %>% 
  mp_cal_abundance( # for each groups 
    .abundance = RareAbundance
    , .group = Fertilization
  )

mp_aba.r.yc <- mp_cal_abundance(mp_ab.r.yc # for each samples
                              , .abundance = RareAbundance
) %>% 
  mp_cal_abundance( # for each groups 
    .abundance = RareAbundance
    , .group = Fertilization
  )

mp_a.r.ac.yc <- mp_cal_abundance(mp_rar.r.ac.yc # for each samples
                                 , .abundance = RareAbundance
) %>% 
  mp_cal_abundance( # for each groups 
    .abundance = RareAbundance
    , .group = Fertilization
  )

mp_a.r.hc.yc <- mp_cal_abundance(mp_rar.r.hc.yc # for each samples
                                 , .abundance = RareAbundance
) %>% 
  mp_cal_abundance( # for each groups 
    .abundance = RareAbundance
    , .group = Fertilization
  )

mp_aab.r <- mp_cal_abundance(mp_a.r # for each samples
                              , .abundance = RareAbundance
) %>% 
  mp_cal_abundance( # for each groups 
    .abundance = RareAbundance
    , .group = Fertilization
  )

mp_aab.bs <- mp_cal_abundance(mp_a.bs # for each samples
                                 , .abundance = RareAbundance
) %>% 
  mp_cal_abundance( # for each groups 
    .abundance = RareAbundance
    , .group = Copper
  )


# 4.3.2) PCA Analysis for relative abundances

# mp_abp.nc <- mp_cal_pca(mp_aba
#                         , .abundance = RelRareAbundanceBySample 
#                         )
# 
# mp_plot_ord(mp_abp.nc
#             ,.ord = pca
#             ,.group = Fertilization
#             ,.color = Fertilization
#             ,.size = Observe
#             ,.alpha = Shannon
#             ,ellipse = F
#             ,show.legend = FALSE # don't display the legend of stat_ellipse
# ) +
#   stat_conf_ellipse(level = 0.95) +
#   scale_fill_manual(values=c("#62C7FA", 
#                              "#0F80FF", 
#                              "#4C4C4C", 
#                              "#A1A3A8")) +
#   scale_color_manual(values=c("#62C7FA", 
#                               "#0F80FF", 
#                               "#4C4C4C", 
#                               "#A1A3A8"),
#                      labels= c('AMH3-8','Hoagland','Control','Suelo NR')) + #labels
#   ggtitle("Sin cobre") +
  # theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
  #       axis.text = element_text(face = "bold", size = 12), # Axis text format
  #       axis.title = element_text(face = "bold", size = 12), # Axis title format
  #       plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
  #       # legend.key.size = unit(0.5, "cm"), # Space between legend elements
  #       legend.text = element_text(size = 12), # Legend text format
  #       legend.title = element_blank())# Remove legend title

# 4.3.2) Relative abundance flow bar plot by phylum

mp_plot_abundance(mp_rar
                  , RareAbundance
                  , .group = Copper
                  , taxa.class = Phylum
)

mp_plot_abundance(mp_aba.nc
                  , RareAbundance
                  , .group = Fertilization
                  , plot.group = TRUE
                  , taxa.class = Phylum
                  , topn = 6 # Plot he most abundant X taxa
                  , rmun = TRUE
                  , rm.zero = TRUE
)+ ggtitle(label = "Sin Cobre") +
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
        axis.text = element_text(face = "bold", size = 12), # Axis text format
        axis.title = element_text(face = "bold", size = 12), # Axis title format
        plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
        # legend.key.size = unit(0.5, "cm"), # Space between legend elements
        legend.text = element_text(size = 12), # Legend text format
        legend.title = element_blank()) + # Remove legend title
  ylab("Abundancia relativa (%)") +
  scale_x_discrete(labels = c("Control", "Hoagland", "AMH3-8", "Suelo NR")) +
mp_plot_abundance(mp_aba.yc
                  , RareAbundance
                  , .group = Fertilization
                  , plot.group = TRUE
                  , taxa.class = Phylum
                  , topn = 6 # Plot he most abundant X taxa
                  , rmun = TRUE
                  , rm.zero = TRUE
) + ggtitle(label = "Con Cobre") +
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
        axis.text = element_text(face = "bold", size = 12), # Axis text format
        axis.title = element_text(face = "bold", size = 12), # Axis title format
        plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
        # legend.key.size = unit(0.5, "cm"), # Space between legend elements
        legend.text = element_text(size = 12), # Legend text format
        legend.title = element_blank()) + # Remove legend title
  ylab("Abundancia relativa (%)") +
  scale_x_discrete(labels = c("Control", "Hoagland", "AMH3-8", "Suelo NR")) 


mp_plot_abundance(mp_aba.r.nc
                   , RareAbundance
                   , .group = Fertilization
                   , plot.group = TRUE
                   , taxa.class = Phylum
                   , rmun = TRUE
                   , rm.zero = TRUE
) +
  mp_plot_abundance(mp_aba.r.yc
                    , RareAbundance
                    , .group = Fertilization
                    , plot.group = TRUE
                    , taxa.class = Phylum
                    , rmun = TRUE
                    , rm.zero = TRUE
)



# 4.3.2) Relative abundance flow bar plot by genus

mp_plot_abundance(mp_aba.nc
                  , RareAbundance
                  , .group = Fertilization
                  , plot.group = TRUE
                  , taxa.class = Genus
                  , topn = 6 # Plot he most abundant X taxa
                  , rmun = TRUE
                  , rm.zero = TRUE
)  + theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
        axis.text = element_text(face = "bold", size = 12), # Axis text format
        axis.title = element_text(face = "bold", size = 12), # Axis title format
        plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
        # legend.key.size = unit(0.5, "cm"), # Space between legend elements
        legend.text = element_text(size = 12), # Legend text format
        legend.title = element_blank()) + # Remove legend title
  ylab("Abundancia relativa (%)") +
  scale_x_discrete(labels = tr_lab) +
  mp_plot_abundance(mp_aba.yc
                    , RareAbundance
                    , .group = Fertilization
                    , plot.group = TRUE
                    , taxa.class = Genus
                    , topn = 6 # Plot he most abundant X taxa
                    , rmun = TRUE
                    , rm.zero = TRUE
  )+ theme(panel.border = element_rect(color ="black", linewidth = 1), # Plot border thickness
           axis.text = element_text(face = "bold", size = 12), # Axis text format
           axis.title = element_text(face = "bold", size = 12), # Axis title format
           plot.title = element_text(face = "bold", size = rel(1.5), margin = margin(b = 15), hjust = 0.5), # Plot title format
           # legend.key.size = unit(0.5, "cm"), # Space between legend elements
           legend.text = element_text(size = 12), # Legend text format
           legend.title = element_blank()) + # Remove legend title
  ylab("Abundancia relativa (%)") +
  scale_x_discrete(labels = tr_lab) 

mp_plot_abundance(mp_aba.r.nc
                  , RareAbundance
                  , .group = Fertilization
                  , plot.group = TRUE
                  , taxa.class = Genus
                  , rmun = TRUE
                  , rm.zero = TRUE
) +
  mp_plot_abundance(mp_aba.r.yc
                    , RareAbundance
                    , .group = Fertilization
                    , plot.group = TRUE
                    , taxa.class = Genus
                    , rmun = TRUE
                    , rm.zero = TRUE
  )

# 4.3.3) Intersectional composition

mp_abau <- mp_cal_upset(mp_aba
                          , .group = Treatment
                          , .abundance = RareAbundance
)

mp_plot_upset(mp_abau
              , .group = Treatment
              , order_by = "degree"
              # , n_intersections = 6
              , intersections = list(as.character(seq(1,8,1)), as.character(seq(1,3,1)), "1", "2", "3", as.character(seq(1,4,1)), as.character(seq(5,8,1)), as.character(seq(5,7,1)), "5", "6", "7")
              , sets = as.character(seq(1,8,1))
              # , label_color= mp_extract_sample(mp_abav.nc)$Fertilization
              # , category.names = c("Control", "Suelo NR", "Hoagland","AMH3-8")
) 

mp_abav.nc <- mp_cal_venn(mp_aba.nc
            , .group = Fertilization
            , .abundance = RareAbundance
            )

mp_abav.yc <- mp_cal_venn(mp_aba.yc
                       , .group = Fertilization
                       , .abundance = RareAbundance
)

mp_abav.r.nc <- mp_cal_venn(mp_aba.r.nc
                          , .group = Fertilization
                          , .abundance = RareAbundance
)

mp_abav.r.yc <- mp_cal_venn(mp_aba.r.yc
                          , .group = Fertilization
                          , .abundance = RareAbundance
)

mp_plot_venn(mp_abav.nc
             , .group = Fertilization
             # , label_color= mp_extract_sample(mp_abav.nc)$Fertilization
             , category.names = c("Control", "Suelo NR", "Hoagland","AMH3-8")
) +
mp_plot_venn(mp_abav.yc
             , .group = Fertilization
             # , label_color= mp_extract_sample(mp_abav.yc)$Fertilization
             , category.names = c("Control", "Suelo NR", "Hoagland","AMH3-8")
) 
# + scale_fill_manual(values=c("#4C4C4C", 
                        # "#A1A3A8",
                        # "#0F80FF", 
                        # "#62C7FA" )) #Standard colors for Control, Hoagland, AMH3-8 and Bulk soil

mp_plot_venn(mp_abav.r.nc
             , .group = Fertilization
             # , label_color= mp_extract_sample(mp_abav.nc)$Fertilization
             , category.names = c("Control","Hoagland","AMH3-8")
) +
  mp_plot_venn(mp_abav.r.yc
               , .group = Fertilization
               # , label_color= mp_extract_sample(mp_abav.yc)$Fertilization
               , category.names = c("Control","Hoagland","AMH3-8")
  ) 


# 4.3.3) Differencial abundance analysis

mp_abaud <- mp_diff_analysis(mp_abau
                 , .abundance = RelRareAbundanceBySample
                 , .group = Fertilization
                 , tip.level = "OTU"
                 , first.test.method = "kruskal.test"
                 , first.test.alpha = 0.05
                 , filter.p = "pvalue"
                 , second.test.method = "wilcox.test"
                 , second.test.alpha = 0.05
                 , cl.min = 2
                 , seed = 1
                 , normalization = 1e+06
                 , ldascore = 2
                 , bootnums = 100
)

mp_abavd.nc <- mp_diff_analysis(mp_abav.nc
                            , .abundance = RelRareAbundanceByFertilization
                            , .group = Fertilization
                            , tip.level = "Genus"
                            , first.test.method = "glm"
                            , first.test.alpha = 0.2
                            , filter.p = "pvalue"
                            , second.test.method = "wilcox.test"
                            , second.test.alpha = 1
                            , cl.min = 2
                            , seed = 1
                            , normalization = 1e+06
                            , ldascore = 2
                            , bootnums = 100
)

mp_abad.r.nc <- mp_diff_analysis(mp_aba.r.nc
                               , .abundance = RelRareAbundanceByFertilization
                               , .group = Fertilization
                               , tip.level = "Genus"
                               , first.test.method = "glm"
                               , first.test.alpha = 0.1
                               , filter.p = "pvalue"
                               , second.test.method = "wilcox.test"
                               , second.test.alpha = 1
                               , cl.min = 2
                               , seed = 1
                               , normalization = 1e+06
                               , ldascore = 2
                               , bootnums = 100
)

mp_abad.r.ac.nc <- mp_diff_analysis(mp_a.r.ac.nc
                                 , .abundance = RelRareAbundanceByFertilization
                                 , .group = Fertilization
                                 , tip.level = "Genus"
                                 , first.test.method = "glm"
                                 , first.test.alpha = 0.05
                                 , filter.p = "pvalue"
                                 , second.test.method = "wilcox.test"
                                 , second.test.alpha = 0.23
                                 , cl.min = 2
                                 , seed = 1
                                 , normalization = 1e+06
                                 , ldascore = 2
                                 , bootnums = 100
)

mp_abad.r.hc.nc <- mp_diff_analysis(mp_a.r.hc.nc
                                    , .abundance = RelRareAbundanceByFertilization
                                    , .group = Fertilization
                                    , tip.level = "Genus"
                                    , first.test.method = "glm"
                                    , first.test.alpha = 0.05
                                    , filter.p = "pvalue"
                                    , second.test.method = "wilcox.test"
                                    , second.test.alpha = 0.23
                                    , cl.min = 2
                                    , seed = 1
                                    , normalization = 1e+06
                                    , ldascore = 2
                                    , bootnums = 100
)

# mp_abavd.yc <- mp_diff_analysis(mp_abav.yc
#                                , .abundance = RelRareAbundanceBySample
#                                , .group = Fertilization
#                                , tip.level = "OTU"
#                                , first.test.method = "kruskal.test"
#                                , first.test.alpha = 1
#                                , filter.p = "pvalue"
#                                , second.test.method = "wilcox.test"
#                                , second.test.alpha = 1
#                                , cl.min = 2
# )

mp_abad.r.yc <- mp_diff_analysis(mp_aba.r.yc
                                 , .abundance = RelRareAbundanceByFertilization
                                 , .group = Fertilization
                                 , tip.level = "Genus"
                                 , first.test.method = "glm"
                                 , first.test.alpha = 0.1
                                 , filter.p = "pvalue"
                                 , second.test.method = "wilcox.test"
                                 , second.test.alpha = 1
                                 , cl.min = 2
                                 , seed = 1
                                 , normalization = 1e+06
                                 , ldascore = 2
                                 , sample.prop.boot = 1
                                 , bootnums = 100
)

mp_abad.r.ac.yc <- mp_diff_analysis(mp_a.r.ac.yc
                                    , .abundance = RelRareAbundanceByFertilization
                                    , .group = Fertilization
                                    , tip.level = "Genus"
                                    , first.test.method = "glm"
                                    , first.test.alpha = 0.05
                                    , filter.p = "pvalue"
                                    , second.test.method = "wilcox.test"
                                    , second.test.alpha = 0.23
                                    , cl.min = 2
                                    , seed = 1
                                    , normalization = 1e+06
                                    , ldascore = 2
                                    , bootnums = 100
                                    , sample.prop.boot = 1
)

mp_abad.r.hc.yc <- mp_diff_analysis(mp_a.r.hc.yc
                                    , .abundance = RelRareAbundanceByFertilization
                                    , .group = Fertilization
                                    , tip.level = "Genus"
                                    , first.test.method = "glm"
                                    , first.test.alpha = 0.05
                                    , filter.p = "pvalue"
                                    , second.test.method = "wilcox.test"
                                    , second.test.alpha = 0.23
                                    , cl.min = 2
                                    , seed = 1
                                    , normalization = 1e+06
                                    , ldascore = 2
                                    , bootnums = 100
                                    , sample.prop.boot = 1
)

mp_aabd.r <- mp_diff_analysis(mp_aab.r
                             , .abundance = RelRareAbundanceByFertilization
                             , .group = Fertilization
                             , tip.level = "Genus"
                             , first.test.method = "kruskal_test"
                             , first.test.alpha = 0.1
                             , filter.p = "pvalue"
                             , second.test.method = "wilcox.test"
                             , second.test.alpha = 0.1
                             , cl.min = 2
                             , seed = 1
                             , normalization = 1e+07
                             , ldascore = 1
                             , bootnums = 100
)

# mp_plot_diff_boxplot(mp_abad.r.nc,
#                                       .group = Fertilization,
#                                       taxa.class = c(Genus) # select the taxonomy level to display
# ) %>%
#   set_diff_boxplot_color(
#     values = c("#62C7FA",
#                "#0F80FF",
#                "#4C4C4C"),
#     guide = guide_legend(title=NULL)
#   )

# 4.3.3.1) Boxplot for differential abundances

mp_plot_diff_boxplot(mp_abaud,
                     .group = Copper,
                     taxa.class = c(Genus) # select the taxonomy level to display
) %>%
  set_diff_boxplot_color(
    values = c("#028202","#DE7920"),
    guide = guide_legend(title=NULL)
  )

mp_plot_diff_boxplot(mp_aabd.r,
                     .group = Fertilization,
                     taxa.class = c(all) # select the taxonomy level to display
) %>%
  set_diff_boxplot_color(
    values = tre_col,
    guide = guide_legend(title=NULL)
  )

diffbox_ac_nc <- mp_plot_diff_boxplot(mp_abad.r.ac.nc,
  .group = Fertilization,
  taxa.class = c(all) # select the taxonomy level to display
) %>%
  set_diff_boxplot_color(
    values = c("#4C4C4C", 
               "#62C7FA"),
    guide = "none"
  )
  
diffbox_hc_nc <- mp_plot_diff_boxplot(mp_abad.r.hc.nc,
                       .group = Fertilization,
                       taxa.class = c(all) # select the taxonomy level to display
  ) %>%
  set_diff_boxplot_color(
    values = c("#4C4C4C", 
               "#0F80FF"),
    guide = "none"
  )

diffbox_ac_yc <- mp_plot_diff_boxplot(mp_abad.r.ac.yc,
                                      .group = Fertilization,
                                      taxa.class = c(all) # select the taxonomy level to display
) %>%
  set_diff_boxplot_color(
    values = c("#4C4C4C", 
               "#62C7FA"),
    guide = "none"
  )

diffbox_hc_yc <- mp_plot_diff_boxplot(mp_abad.r.hc.yc,
                                      .group = Fertilization,
                                      taxa.class = c(all) # select the taxonomy level to display
) %>%
  set_diff_boxplot_color(
    values = c("#4C4C4C", 
               "#0F80FF"),
    guide = "none"
  )

plot_list(diffbox_ac_nc, diffbox_hc_nc, diffbox_ac_yc, diffbox_hc_yc, byrow = F)

diffbox_r_nc <- mp_plot_diff_boxplot(mp_abad.r.nc,
                                      .group = Fertilization,
                                      group.abun = F,
                                      taxa.class = c(Genus,Phylum) # select the taxonomy level to display
) %>%
  set_diff_boxplot_color(
    values = c("#4C4C4C", 
               "#0F80FF",
               "#62C7FA"),
    guide = guide_legend(title=NULL)
  )

diffbox_r_yc <- mp_plot_diff_boxplot(mp_abad.r.yc,
                                      .group = Fertilization,
                                      group.abun = F,
                                      taxa.class = c(Genus,Phylum) # select the taxonomy level to display
) %>%
  set_diff_boxplot_color(
    values = c("#4C4C4C", 
               "#0F80FF",
               "#62C7FA"),
    guide = guide_legend(title=NULL)
  )

plot_list(diffbox_r_nc, diffbox_r_yc, byrow = T)

# 4.3.3.3) Radial tree plot for differential abundances


mp_aabd.r %>%
  mp_plot_diff_res(
    group.abun = TRUE,
    pwidth.abun=0.05,
    offset.effsize = 0.6
  ) +
  scale_fill_manual(values=tre_col) +
  scale_fill_manual(
    aesthetics = "fill_new", # The fill aes was renamed to "fill_new" for the abundance dotplot layer
    values = tre_col
  ) +
  scale_fill_manual(
    aesthetics = "fill_new_new", # The fill aes for hight light layer of tree was renamed to 'fill_new_new'
    values = phy_col
  )

biom_tree.ac.nc <- mp_abad.r.ac.nc %>%
  mp_plot_diff_res(
    group.abun = TRUE,
    pwidth.abun=0.05,
    offset.effsize = 0.6
  ) +
  scale_fill_manual(values=c("#62C7FA", 
                             "#4C4C4C")) +
  scale_fill_manual(
    aesthetics = "fill_new", # The fill aes was renamed to "fill_new" for the abundance dotplot layer
    values = c("#4C4C4C", 
               "#62C7FA")
  ) +
  scale_fill_manual(
    aesthetics = "fill_new_new", # The fill aes for hight light layer of tree was renamed to 'fill_new_new'
    values = c("#E41A1C", "#377EB8", "#4DAF4A",
               "#984EA3", "#FF7F00", "#FFFF33",
               "#A65628", "#F781BF", "#999999",
               "#A65628", "#F781BF", "#999999",
               "#A65628"
    )
  )

biom_tree.hc.nc <- mp_abad.r.hc.nc %>%
  mp_plot_diff_res(
    group.abun = TRUE,
    pwidth.abun=0.05,
    offset.effsize = 0.6
  ) +
  scale_fill_manual(values=c("#0F80FF", 
                             "#4C4C4C")) +
  scale_fill_manual(
    aesthetics = "fill_new", # The fill aes was renamed to "fill_new" for the abundance dotplot layer
    values = c("#4C4C4C", 
               "#0F80FF")
  ) +
  scale_fill_manual(
    aesthetics = "fill_new_new", # The fill aes for hight light layer of tree was renamed to 'fill_new_new'
    values = c("#E41A1C", "#377EB8", "#4DAF4A",
               "#984EA3", "#FF7F00", "#FFFF33",
               "#A65628", "#F781BF", "#999999",
               "#A65628", "#F781BF", "#999999",
               "#A65628"
    )
  )

biom_tree.ac.yc <- mp_abad.r.ac.yc %>%
  mp_plot_diff_res(
    group.abun = TRUE,
    pwidth.abun=0.01,
    offset.effsize = 0.6
  ) +
  scale_fill_manual(values=c("#62C7FA", 
                             "#4C4C4C")) +
  scale_fill_manual(
    aesthetics = "fill_new", # The fill aes was renamed to "fill_new" for the abundance dotplot layer
    values = c("#4C4C4C", 
               "#62C7FA")
  ) +
  scale_fill_manual(
    aesthetics = "fill_new_new", # The fill aes for hight light layer of tree was renamed to 'fill_new_new'
    values = c("#E41A1C", "#377EB8", "#4DAF4A",
               "#984EA3", "#FF7F00", "#FFFF33",
               "#A65628", "#F781BF", "#999999",
               "#A65628", "#F781BF", "#999999",
               "#A65628"
    )
  )

biom_tree.hc.yc <- mp_abad.r.hc.yc %>%
  mp_plot_diff_res(
    group.abun = TRUE,
    pwidth.abun=0.05,
    offset.effsize = 0.6
  ) +
  scale_fill_manual(values=c("#0F80FF", 
                            "#4C4C4C")) +
  scale_fill_manual(
    aesthetics = "fill_new", # The fill aes was renamed to "fill_new" for the abundance dotplot layer
    values = c("#4C4C4C", 
               "#0F80FF")
  ) +
  scale_fill_manual(
    aesthetics = "fill_new_new", # The fill aes for hight light layer of tree was renamed to 'fill_new_new'
    values = c("#E41A1C", "#377EB8", "#4DAF4A",
               "#984EA3", "#FF7F00", "#FFFF33",
               "#A65628", "#F781BF", "#999999",
               "#A65628", "#F781BF", "#999999",
               "#A65628"
    )
  )

ggarrange(biom_tree.ac.nc, biom_tree.ac.yc, biom_tree.hc.nc, biom_tree.hc.yc, nrow = 2, ncol = 2, align = "hv", common.legend = F)

biom_tree.r.nc <- mp_abad.r.nc %>%
  mp_plot_diff_res(
    group.abun = TRUE,
    pwidth.abun=0.05,
    offset.effsize = 0.6
  ) +
  scale_fill_manual(values = c("#62C7FA", 
                               "#0F80FF",
                               "#4C4C4C")) +
  scale_fill_manual(
    aesthetics = "fill_new", # The fill aes was renamed to "fill_new" for the abundance dotplot layer
    values = c("#4C4C4C", 
               "#0F80FF",
               "#62C7FA")
  ) +
  scale_fill_manual(
    aesthetics = "fill_new_new", # The fill aes for hight light layer of tree was renamed to 'fill_new_new'
    values = c("#E41A1C", "#377EB8", "#4DAF4A",
               "#984EA3", "#FF7F00", "#FFFF33",
               "#A65628", "#F781BF", "#999999",
               "#A65628", "#F781BF", "#999999",
               "#A65628"
    )
  )

biom_tree.r.yc <- mp_abad.r.yc %>%
  mp_plot_diff_res(
    group.abun = TRUE,
    pwidth.abun=0.01,
    offset.effsize = 0.6
  ) +
  scale_fill_manual(values = c("#62C7FA", 
                               "#0F80FF",
                               "#4C4C4C")) +
  scale_fill_manual(
    aesthetics = "fill_new", # The fill aes was renamed to "fill_new" for the abundance dotplot layer
    values = c("#4C4C4C", 
               "#0F80FF",
               "#62C7FA")
  ) +
  scale_fill_manual(
    aesthetics = "fill_new_new", # The fill aes for hight light layer of tree was renamed to 'fill_new_new'
    values = c("#E41A1C", "#377EB8", "#4DAF4A",
               "#984EA3", "#FF7F00", "#FFFF33",
               "#A65628", "#F781BF", "#999999",
               "#A65628", "#F781BF", "#999999",
               "#A65628"
    )
  )

ggarrange(biom_tree.r.nc, biom_tree.r.yc, nrow = 1, ncol = 2, align = "hv", common.legend = F)

# 5) Metabolic inference

# 5.1) Exporting files for PICRUSt2

library(biomformat)

ps_diff_ac_nc <- as.phyloseq(mp_a.r.ac.nc)
bm_diff_ac_nc <- make_biom(data = as.matrix(otu_table(ps_diff_ac_nc))
                           , sample_metadata = meta(ps_diff_ac_nc))
write_biom(bm_diff_ac_nc, "biom_ac_nc.biom")
seqs_diff_ac_nc <- rownames(otu_table(ps_diff_ac_nc))
write.fasta(as.list(seqs_diff_ac_nc), names = seqs_diff_ac_nc, file.out = "seqs_ac_nc.fna")


bm_r_nc <- make_biom(data = as.matrix(t(otu_table(ps_r_nocu))))
                           # , sample_metadata = meta(ps_r_nocu))
write_biom(bm_r_nc, "biom_r_nc.biom")
seqs_r_nc <- rownames(t(otu_table(ps_r_nocu)))
write.fasta(as.list(seqs_r_nc), names = seqs_r_nc, file.out = "seqs_r_nc.fna")

bm_r_yc <- make_biom(data = as.matrix(t(otu_table(ps_r_yescu))))
                     # , sample_metadata = meta(ps_r_yescu))
write_biom(bm_r_yc, "biom_r_yc.biom")
seqs_r_yc <- rownames(t(otu_table(ps_r_yescu)))
write.fasta(as.list(seqs_r_yc), names = seqs_r_yc, file.out = "seqs_r_yc.fna")


## PICRUSt2 runs in conda3. So the output has to be imported back

# 5.2) Processing PICRUSt2 output

library(ggpicrust2)
library("metagenomeSeq")
library("ALDEx2")
library("lefser")
library("Maaslin2")

picrustres <- ggpicrust2(file = "pred_metagenome_unstrat.tsv",
           metadata = meta(ps_diff_ac_nc),
           group = "Fertilization",
           pathway = "EC",
           daa_method = "LinDA",
           ko_to_kegg = TRUE,
           order = "pathway_class",
           p_values_bar = TRUE,
           x_lab = "pathway_name")

pc_res_r_nc <- ggpicrust2(file = paste(getwd(),"/r_nc_picrust/EC_metagenome_out/pred_metagenome_unstrat.tsv", sep = ""),
                         metadata = meta(ps_r_nocu),
                         group = "Fertilization",
                         pathway = "EC",
                         daa_method = "LinDA",
                         ko_to_kegg = F,
                         order = "pathway_class",
                         p_values_bar = TRUE,
                         x_lab = "pathway_name",
                         reference = "No"
                         )


daa_results_list <- ggpicrust2(file = "pred_metagenome_unstrat.tsv",
    metadata = meta(ps_diff_ac_nc),
    group = "Fertilization",
    pathway = "EC",
    daa_method = "LinDA",
    p_values_bar = TRUE,
    # order = "Fertilization",
    ko_to_kegg = FALSE,
    x_lab = "description",
    p.adjust = "BH",
    select = NULL,
    reference = NULL
  )