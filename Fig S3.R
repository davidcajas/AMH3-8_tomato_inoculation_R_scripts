## 1) DESCARGAR DEPENDENCIAS

# install.packages("textshape")
# install.packages("devtools")
# library(devtools)
# install_github("vqv/ggbiplot")
# install.packages("ggrepel")          # Install ggrepel package


## 2) CARGAR LIBRERÍAS
library(readxl) # necesario para leer archivos Excel
library(textshape) # herramienta para modificar el dataframe
library(stats) # contiene función prcomp
library(FactoMineR) # contiene función PCA
library(factoextra) # herramientas de visualización y exploración de datos en objeto PCA
library("ggrepel") # añadir texto sin sobrelape
library(ggbiplot)

## 3) CARGAR Y LIMPIAR DATOS

setwd("/Users/davidcajasmunoz/Library/CloudStorage/GoogleDrive-dadavid.cajas@gmail.com/Mi unidad/Academia/Postgrado/PUCV/Tesis/Resultados/Biofertilización/2/4-Cultivos/Ecoplates/R/PCA")
dataframe <- read.delim2("dataframe.txt", header=TRUE)
# View(dataframe)

# Si es necesario colocar la primera columna como nombre de filas. En este caso con números está bien.
# dataframe<-textshape::column_to_rownames(dataframe, loc = 1)

# Generar un dataframe a partir del archivo
df<-as.data.frame(dataframe)

# Quitar del df las variables no numéricas
df <- subset(df, select = -c(Group,Treatment,Copper) )
# View(df)
# Quitar del df las variables no globales de diversidad, dejando sólo las variables individuales
df <- subset(df, select = -c(Richness..R.,Diversity..H..,Evenness..E..,AWCD) )
# View(df)
# print(df)

## 4) PCA

# PCA con el paquete stats

respca <- prcomp(x = df, center = TRUE, scale. = TRUE)

# Crear objetos factor a partir del dataframe original para incluir las variables categóricas en el gráfico

# dataframe[dataframe["Copper"] == "Sin cobre", "Copper"] = "No Copper"
# dataframe[dataframe["Copper"] == "Con cobre", "Copper"] = "Copper"
# 
# dataframe[dataframe["Treatment"] == "Suelo NR", "Treatment"] = "Bulk soil"

Copper<-as.factor(dataframe$Copper)
Treatment<-as.factor(dataframe$Treatment)
Group<-as.factor(dataframe$Group)

# PCA con el paquete FactoMineR
respca2 <- PCA(X = df, scale.unit = TRUE, ncp = 6, graph = TRUE)

# Fijar colores para variables

# Treatments
group_col <- c('Control Cu-' = "#4C4C4C",
               'Control Cu+' =  "#4C4C4C",
             'Hoagland Cu-'= "#0F80FF",
             'Hoagland Cu+'= "#0F80FF",
             'AMH3-8 Cu-' = "#62C7FA", 
             'AMH3-8 Cu+' = "#62C7FA", 
             'Bulk soil Cu-' = "#A1A3A8", 
             'Bulk soil Cu+' = "#A1A3A8")
group_lab <- c('Control Cu-' = "Control",
               'Control Cu+' =  "",
               'Hoagland Cu-'= "Hoagland",
               'Hoagland Cu+'= "",
               'AMH3-8 Cu-' = "AMH3-8", 
               'AMH3-8 Cu+' = "", 
               'Bulk soil Cu-' = "Bulk soil", 
               'Bulk soil Cu+' = "")

tr_leg <- legend(x="bottom", horiz = TRUE, legend = c('Control','Hoagland','AMH3-8','Bulk soil'), bty = "n", pch=22 , fill=c("#4C4C4C", "#0F80FF", "#62C7FA", "#A1A3A8"), border = "black", text.col = "black", pt.lwd = 0, cex=1.1, pt.cex=1.1)

# Copper
cop_col <- c('No Copper' = "#028202",
             'Copper' = "#DE7920")

# Carbon source

CS <- rownames(t(df))
colors_cs <- c("ca","pol","pol","ch.p","ch.p","ch","ch","ch","ch","ch","ch","ch","ca.aa","ch","ch","ch.ca","ch.ca","ca","ca","ca","ca","ca","ca","aa","aa","aa","aa","aa","aa","a","a")
palette <- colorRampPalette(colors = c("#b55eb2", "#148F77"), space = "Lab")(length(unique(colors_cs)))
color_mapping <- setNames(palette, unique(colors_cs)[c(2,3,4,6,1,5,7,8)])
color_mapping[8] <- "darkred"

## 5) EXPLORACIÓN DE DATOS.

# Ver ejemplos de como visualizar aquí:
# http://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining
print(respca2) # Ver variables entregadas por el PCA

head(respca2$eig) # Proporción de la varianza explicada por cada componente
fviz_eig(respca2) # Visualizar eigenvalores
fviz_contrib(respca2,choice = "var") # Representa la contribución de variables a a los PC

get_pca_var(respca2) # Qué información podemos extraer sobre las variables?

get_pca_ind(respca2) # Qué información podemos extraer sobre las observaciones?

fviz_pca_ind(respca2, # Grafica las observaciones en los ejes PC1 y PC2.
             col.ind = "cos2", # Color según calidad de representación
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # pseudocolor de cálido a frío
             repel = FALSE     # No evitar sobrelape de texto
)

fviz_pca_var(respca2, # Grafica las variables del df en los ejes PC1 y PC2.
             col.var = "contrib", # Color según contribución a los PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # pseudocolor de cálido a frío
             repel = TRUE     # Evitar sobrelape de texto
)

fviz_pca_biplot(respca2, repel = TRUE, # Representa en un mismo gráfico las observaciones y las variables
                col.var = "#2E9FDF", # Un color para las variables
                geom.ind = "point", # Dibujar observaciones sin colocar etiquetas
                col.ind = "#696969"  # Un color para las observaciones
)

## 6) GRÁFICO FINAL.

# 6.1)
# Creé una función ggbiplot modificada, a la cual nombré ggbiplot2. Esta incluye la posibilidad de modificar los siguientes parámetros:
# arrow.color = color de las flechas de variables
# arrow.linetype = tipo de línea de las flechas de variables
# arrow.alpha = transparencia de las flechas de variables 
# arrow.size = escalador del tamaño de los vectores de flechas de variables
# varname.angle = ángulo de la etiqueta de las flechas de variables
# varname.color = color de la etiqueta de las flechas de variables
# varname.border = agrega borde a la etiqueta de las flechas de variables 
# varname.fill = agrega relleno a la etiqueta de las flechas de variables
# varname.ff = elige fontface de la etiqueta de las flechas de variables
# ellipse.linetype = tipo de línea de los elipses
# ellipse.linewidth = ancho de línea de los elipses
# ellipse.fill = aún no funciona
# intercept = si TRUE, coloca líneas en intercepto 0,0

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
  library(ggrepel)
  
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
      geom_label_repel(data = df.v, 
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

# 6.1) LE GRÁFICO

FS2 <- ggbiplot2(respca, 
         groups = dataframe$Group, # Procesamiento en función de factor Group en dataframe, que incluye condicoines y tratamientos
         ellipse = TRUE, # Crear elipses
         ellipse.prob = 0.70, # Intervalo de confianza para la elipse
         ellipse.linewidth = 0.5, # Ancho de línea para la elipse
         obs.scale = 1, # Ajuste escala de puntos
         var.scale = 1, # De alguna forma modifica la escala del eje X. Ni pico idea cómo funciona
         arrow.size = 1.85, # Escala de las flechas de variables para facilitar visualización
         arrow.color = color_mapping[colors_cs], # Color de flechas
         arrow.alpha = 0.7, # opacidad de flechas
         varname.color = color_mapping[colors_cs], # Color de etiquetas de variables
         varname.adjust = 1.2, # Aleja un poco las etiquetas de las flechas
         varname.size = 2.6, # Tamaño de las etiquetas de las variables
         varname.ff = "bold", # etiquetas en negrita
         intercept = TRUE
) + theme_bw() +
  scale_color_manual(name="Ellipses", values=group_col, #paleta de colores para AMH3-8, Control, Hoagland y Suelo NR
                      labels= group_lab) + 
  scale_shape_manual(values = c(
                                'AMH3-8' = 21, 
                                'Control' = 22, 
                                'Hoagland' = 23, 
                                'Bulk soil' = 24)
                     ) + #formas para los tratamientos
  scale_fill_manual(name = "Soil", values = c("#DE7920", 
                               "#028202")) + #paleta de colores para con cobre, sin cobre
  geom_point(aes(shape = Treatment, fill = Copper), size = 3) + # Fija objeto Tratamientos para las formas y Cobre para los colores de relleno
  guides(fill = guide_legend(override.aes = list(shape = 21) ), # Cambia formas de leyenda de "Fill" a formas con relleno para que muestre los colores
         shape = guide_legend(override.aes = list(fill = "black") ) ) + # Cambia el relleno de las formas de "Shape" para que se vean negras
  theme(panel.border = element_rect(color ="black", linewidth = 1), # Borde en todo el gráfico
        axis.text = element_text(face = "bold", size = 12), # Formato fuente de los números de los ejes
        axis.title = element_text(face = "bold", size = 12), # Formato fuente de los títulos de los ejes
        legend.key.size = unit(0.5, "cm")) # Espacios entre legendas
print(FS2)
# exportar en 9 x 8 in portrait
  
  