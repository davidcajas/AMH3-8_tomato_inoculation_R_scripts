## 1) DESCARGAR DEPENDENCIAS

# install.packages("fmsb")

## 2) CARGAR LIBRERÍAS

library(fmsb) # paquete para spider plot

## 3) CARGAR Y LIMPIAR DATOS

#df.no.cu <- read.delim2("~/Library/CloudStorage/OneDrive-usach.cl/Académicos/Postgrado/PUCV/Tesis/Resultados/Biofertilización/2/4-Cultivos/Ecoplates/R/Diversity/diversidad_no-cu.txt")

df.no.cu <- read.delim2("diversidad_no-cu.txt")

#df.yes.cu <- read.delim2("~/Library/CloudStorage/OneDrive-usach.cl/Académicos/Postgrado/PUCV/Tesis/Resultados/Biofertilización/2/4-Cultivos/Ecoplates/R/Diversity/diversidad_no-cu.txt")

df.yes.cu <- read.delim2("diversidad_yes-cu.txt")

# df.no.cu <- df.no.cu[,-1]

rownames(df.no.cu) <- c("Control",
                        "Hoagland",
                        "AMH3-8",
                        "Bulk soil")

colnames(df.no.cu) <- c("AWCD",
                        "H'",
                        "R",
                        "E")
# View(df.no.cu)

# df.yes.cu <- df.yes.cu[,-1]

#Fijar nombres de filas y columnas

rownames(df.yes.cu) <- c("Control",
                        "Hoagland",
                        "AMH3-8",
                        "Bulk soil")

colnames(df.yes.cu) <- c("AWCD",
                        "H'",
                        "R",
                        "E")
# View(df.yes.cu)

## 4) Spider plot

# Hay que agregar líneas para max y min de los valores de cada variable en el dataframe.
df.no.cu <- rbind(rep(1,5) , rep(0,5) , df.no.cu)

df.yes.cu <- rbind(rep(1,5) , rep(0,5) , df.yes.cu)

# Plot

par(family = "Arial", font=2, mar = c(0, 1, 0, 0.5), mfrow = c(1,2)) #Matriz de 1x2, 

colors <- c("#4C4C4C", 
            "#0F80FF", 
            "#62C7FA", 
            "#A1A3A8") #paleta de colores para Control, Hoagland, AMH3-8 y Suelo NR 
#gráfico sin cobre
radarchart(df.no.cu  , axistype=1 , 
            pcol= colors, #colores antes fijados
            pfcol=NULL , plwd=3.5 , plty=1, # líneas segmentadas de ancho 3
            cglcol="black", cglty=1, axislabcol="White", caxislabels=seq(0,1,0.25), cglwd=0.45, # sin leyenda de eje interno
            vlcex=1 #tamaño de letra ejes variables
) 
title(main="No Copper", line = -3) # título de subgráfico
# text(-1, 1.7, "A", cex=1.5) # Letra de subgráfico
#gráfico con cobre
radarchart(df.yes.cu  , axistype=1 , 
           pcol= colors, #colores antes fijados
           pfcol=NULL , plwd=3.5 , plty=1, # líneas segmentadas de ancho 3
           cglcol="black", cglty=1, axislabcol="White", caxislabels=seq(0,1,0.25), cglwd=0.45, # sin leyenda de eje interno
           vlcex=1 #tamaño de letra ejes variables
) 
title(main="Copper", line = -3) # título de subgráfico
# text(-1, 1.7, "B", cex=1.5) # Letra de subgráfico
#leyenda
par(family = "Arial", font=1, mar = c(0,0,0,0), oma = c(3, 0, 0, 0), fig = c(0, 1, 0, 1), mfrow = c(1,1), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(x="bottom", horiz = TRUE, legend = rownames(df.yes.cu[-c(1,2),]), bty = "n", pch=20 , col=colors, text.col = "black", cex=1.1, pt.cex=3)
# Esport as 8 x 6 inch. in landscape mode
