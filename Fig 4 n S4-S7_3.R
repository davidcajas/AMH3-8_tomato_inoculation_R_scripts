# Auxiliar objects and functions

library(ggprism)  

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
            'No_plant' = 'Bulk Soil')
tr_leg <- legend(x="bottom", horiz = TRUE, legend = c('Control','Hoagland','AMH3-8','Bulk Soil'), bty = "n", pch=22 , fill=c("#4C4C4C", "#0F80FF", "#62C7FA", "#A1A3A8"), border = "black", text.col = "black", pt.lwd = 0, cex=1.1, pt.cex=1.1)

# Copper
cop_col <- c('FALSE' = "#028202",
            'TRUE' = "#DE7920")
cop_lab <-c('FALSE' = 'No Copper',
            'TRUE' = 'Copper')
mp_extract_abundance(mp_aba)
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

# Figures

# F4

# A

F4A <- ggarrange(pcoa_nc,pcoa_yc, nrow=1, ncol=2, legend = "none", common.legend = TRUE)

# B

F4B <- ggarrange(ab_genphyl_h_parts) # o ab_genphyl_v o ab_genphyl_h

# C

F4C <- ggarrange(biomarker) # o biomarker_phy

# Plot

ggarrange(F4A
          ,F4B
          ,biomarker_phy
          , ncol=1
          , align = "hv"
          , heights = c(0.9,1.5,1.1)
          , labels = "AUTO", font.label = list(size = 28))
par(family = "Arial", font=1, mar = c(0,0,1,0), oma = c(0.1, 0, 0, 0), fig = c(0, 1, 0, 1), mfrow = c(1,1), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend(x="bottom", horiz = TRUE, legend = c('Control','Hoagland','AMH3-8','Bulk soil'), bty = "n", pch=22 , fill=c("#4C4C4C", "#0F80FF", "#62C7FA", "#A1A3A8"), border = "black", text.col = "black", pt.lwd = 0, cex=1.1, pt.cex=1.1)

# export in 16 x 20 in portrait
  
# FS4

# A
FS4A <- ggarrange(adiv_nocu  + pr + nx + guides(y = "prism_offset_minor")
                  , adiv_yescu  + pr + ny + nx, ncol=2
                  , legend = "bottom", common.legend = TRUE, hjust = -0.7)

# B

FS4B <- (upsetbars| ((upset_nocu/upset_yescu))) # alt (upsetbars/ (upset_nocu+upset_yescu))

# Plot
pdf("FS4.pdf", width = 14, height = 16)
png("FS4.png", width = 2800, height = 3200, res = 200)
ggarrange(FS4A + theme(plot.margin = margin(0,300,0,300))
          , FS4B
          , nrow=2
          , heights = c(1,2.2)
          , labels = "AUTO", font.label = list(size = 28))
dev.off()


# export in 14 x 14 inch

# FS5

# Plot

# ggarrange(coreplot_nocu
#           , coreplot_yescu 
#           , coreplot_genus_nocu
#           , coreplot_genus_yescu
#           , nrow=2
#           , ncol=2
#           , heights = c(1,2)
#           , labels = c("A","","B",""), font.label = list(size = 26), align = "hv")

# export in 15 x 15 inch
# OR

ggarrange(coreplot_phy
          , coreplot_gen 
          , nrow=2
          # , ncol=2
          , heights = c(0.5,2)
          , common.legend = T
          , legend = "right"
          , labels = "AUTO", font.label = list(size = 26), align = "hv")

# export in 10 x 18 inch portrait

# FS6

# Plot

ggarrange(clustering_map
          , pcoa_all
          , ncol=2
          , labels = "AUTO", font.label = list(size = 18))

# export in 4 x 10 inch landscape

# FS7

# Plot

ggarrange(as.ggplot(diffbox_ac_nc) + theme(plot.margin = margin(t = 20))
          , as.ggplot(diffbox_ac_yc)
          , as.ggplot(diffbox_hc_nc)
          , as.ggplot(diffbox_hc_yc)
          ,NULL,NULL
          , nrow=3
          , ncol=2
          , heights = c(1,1,0.1)
          , labels = c("No Copper","Copper","",""), font.label = list(size = 20), label.x = 0.5, align = "h")
par(family = "Arial", font=1, mar = c(0,0,1,0), oma = c(0.1, 0, 0, 0), fig = c(0, 1, 0, 1), mfrow = c(1,1), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(x="bottom", horiz = TRUE, legend = c('Control','Hoagland','AMH3-8'), bty = "n", pch=22 , fill=c("#4C4C4C", "#0F80FF", "#62C7FA"), border = "black", text.col = "black", pt.lwd = 0, cex=1.1, pt.cex=1.1)

# OR
# ggarrange(diffplot_ac_nocu + theme(plot.margin = margin(t = 20))
#           , diffplot_ac_yescu
#           , diffplot_hc_nocu
#           , diffplot_hc_yescu
#           ,NULL,NULL
#           , nrow=3
#           , ncol=2
#           , heights = c(1,1,0.1)
#           , labels = c("Sin cobre","Con cobre","",""), font.label = list(size = 20), label.x = 0.5, align = "h")
# par(family = "Arial", font=1, mar = c(0,0,1,0), oma = c(0.1, 0, 0, 0), fig = c(0, 1, 0, 1), mfrow = c(1,1), new = TRUE)
# plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
# legend(x="bottom", horiz = TRUE, legend = c('Control','Hoagland','AMH3-8'), bty = "n", pch=22 , fill=c("#4C4C4C", "#0F80FF", "#62C7FA"), border = "black", text.col = "black", pt.lwd = 0, cex=1.1, pt.cex=1.1)

# export in 15 x 15 inch