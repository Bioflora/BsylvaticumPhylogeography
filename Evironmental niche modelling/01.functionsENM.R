computeAUC <- function(p, a) {
  p <- stats::na.omit(p)
  a <- stats::na.omit(a)
  np <- length(p)
  na <- length(a)
  
    if (length(p) > 1000) {
      tr <- as.vector(quantile(p, 0:1000/1000))
    } else {
      tr <- p
    }
    if (length(a) > 1000) {
      tr <- c(tr, as.vector(quantile(a, 0:1000/1000)))
    } else {
      tr <- c(tr, a)
    }
    tr <- sort(unique( round(tr, 8)))
    tr <- c( tr - 0.0001, tr[length(tr)] + c(0, 0.0001))
  
  N <- na + np
  
  R <- sum(rank(c(p, a))[1:np]) - (np*(np+1)/2)
  AUC <- R / (as.numeric(na) * as.numeric(np))
  
  return(AUC)
  
}



#FUNCIÓN PARA EVALUAR MODELOS MAXENT MEDIANTE BOOTSTRAP
bootstrapMaxent <- function(data = NULL,
                            columna.presencia = NULL,
                            columnas.variables = NULL,
                            iteraciones = 3,
                            porcentaje.evaluacion = 25,
                            variables = NULL){

  require(raster)
  require(dismo)
  require(maxnet)

  #cambiando nombre del objeto data
  presencia.bootstrap<-data
  
  #separamos los datos de presencia de los de background
  datos.presencia <- data[data[,columna.presencia] == 1, ]
  datos.background  <- data[data[,columna.presencia] == 0, ]
  
  #computamos el número de filas a usar como test en cada iteracion
  datos.presencia.test.n <- round((porcentaje.evaluacion*nrow(datos.presencia))/100, 0)
  datos.background.test.n <- round((porcentaje.evaluacion*nrow(datos.background))/100, 0)
  
  #objetos para guardar auc y mapas
  auc.bootstrap<-vector()
  mapas.bootstrap<-brick()
  
  #iteramos por los grupos
  for (iteracion in 1:iteraciones){
    
    print(paste("Iteracion:", iteracion))
    
    #muestreo aleatorio de los datos de presencia (esto solo retorna índices de casos que se usan más abajo)
    sampling.presencia<-sample(x=nrow(datos.presencia), size=datos.presencia.test.n, replace=TRUE)
    sampling.background<-sample(x=nrow(datos.background), size=datos.background.test.n, replace=TRUE)
    
    #seleccionamos los datos de entrenamiento
    presencias.entrenamiento<-datos.presencia[-sampling.presencia, ]
    background.entrenamiento<-datos.background[-sampling.background, ]
    #los unimos
    datos.entrenamiento <- rbind(presencias.entrenamiento, background.entrenamiento)
    
    #seleccionamos los datos de evaluacion
    presencias.evaluacion<-datos.presencia[-sampling.presencia, ]
    background.evaluacion<-datos.background[-sampling.background, ]
    #los unimos
    datos.evaluacion <- rbind(presencias.evaluacion, background.evaluacion)

    #ajustamos el modelo (cualquier modelo puede ir aquí)
    modelo.entrenado <- maxnet::maxnet(p=datos.entrenamiento[,columna.presencia], data=datos.entrenamiento[, columnas.variables])
    
    #predecimos a mapa y lo guardamos en el brick
    mapas.bootstrap <- raster::addLayer(mapas.bootstrap, raster::predict(variables, modelo.entrenado, type="logistic"))
    
    #calculamos los valores predichos por el modelo para las presencias
    prediccion.presencias <- raster::predict(modelo.entrenado, datos.evaluacion[datos.evaluacion[,columna.presencia]==1, ], type="logistic")
    prediccion.background <- raster::predict(modelo.entrenado, datos.evaluacion[datos.evaluacion[,columna.presencia]==0, ], type="logistic")
    
    #calculamos AUC y lo guardamos
    auc.bootstrap[iteracion]<-computeAUC(p=prediccion.presencias, a=prediccion.background)
    
  } #final de las iteraciones
  
  #preparando resultados
  output.list <- list()
  output.list$auc.mean <- mean(auc.bootstrap)
  output.list$auc.sd <- sd(auc.bootstrap)
  output.list$models.mean <- calc(mapas.bootstrap, mean)
  output.list$models.sd <- calc(mapas.bootstrap, sd)
  
  return(output.list)
  
}#final de la función



  
  
  
  
nicheOverlapFromMaps <- function(map.1, map.2){

    pred1 <- na.omit(as.vector(map.1))
    pred2 <- na.omit(as.vector(map.2))

   d <- 1 - sum(abs(pred1/sum(pred1) - pred2/(sum(pred2))))/2
   i <- 1 - sum((sqrt(pred1/sum(pred1)) - sqrt(pred2/sum(pred2)))**2)/2
   r <- cor(pred1, pred2, method = "pearson")

  output <- list(D = d,
                 I = i,
                 env.cor = r)


return(output)

}




#Applica niche overlap analysis y niche equivalency test a dos especies A y B
nicheOverlap = function(A, B, iterations=100){

  #PCA con todos los datos (fíjate que usamos names(bioclimaticas) para seleccionar solo las variables predictivas)
  pca.A.B = dudi.pca(rbind(presencia.A ,presencia.B)[,names(bioclimaticas)], scannf=F, nf=2)

  #PCA scores area completa
  scores.A.B = pca.A.B$li

  #PCA scores para la presencia de A
  scores.A.presencia = suprow(pca.A.B, presencia.A[which(presencia.A$presencia==1), names(bioclimaticas)])$li

  # PCA scores para el area de A
  scores.A.area = suprow(pca.A.B, presencia.A[,names(bioclimaticas)])$li

  #PCA scores presencias NA
  scores.B.presencia = suprow(pca.A.B, presencia.B[which(presencia.B$presencia==1), names(bioclimaticas)])$li

  # PCA scores para el area NA
  scores.B.area = suprow(pca.A.B, presencia.B[,names(bioclimaticas)])$li

  #PCA scores de cada sitio a grid
  scores.A.grid = ecospat.grid.clim.dyn(glob=scores.A.B, glob1=scores.A.area, sp=scores.A.presencia, R=100, th.sp=0)
  scores.B.grid = ecospat.grid.clim.dyn(glob=scores.A.B, glob1=scores.B.area, sp=scores.B.presencia, R=100, th.sp=0)

  #plot
  ecospat.plot.niche.dyn(scores.A.grid, scores.B.grid, quant=0.25, interest=2, title= "Niche Overlap", name.axis1="PC1", name.axis2="PC2")

  #niche equivalency test
  niche.equivalency.test <- ecospat.niche.equivalency.test(scores.A.grid, scores.B.grid, rep=iterations, alternative = "greater")

  #print
  print(paste("Schoener's D = ", round(niche.equivalency.test$obs$D, 4), "; p-value = ", round(niche.equivalency.test$p.D, 4), sep=""))
  print(paste("Niche equivalency = ", round(niche.equivalency.test$obs$I, 4), "; p-value = ", round(niche.equivalency.test$p.I, 4), sep=""))

  return(niche.equivalency.test)

}



#función para calcular niche breadth según Levin (1968).
#devuelve una lista con los valores de B1 y B2
#model map es un mapa raster con un modelo de distribución
nicheBreadth = function(model.map){

  #estandarizar mapa
  model.map.standardized = model.map/cellStats(model.map, stat=sum)

  #convertir valores 0 a "casi 0" 10*e-40
  model.map.standardized[which(getValues(model.map.standardized) == 0)] = 1e-40

  #calcular B1
  #-----------
  x = model.map.standardized
  x = x[!is.na(x)]
  x = x/sum(x)

  #reemplazar valores que quedan por debajo de la precisión del ordenador
  x[x < .Machine$double.xmin] = .Machine$double.xmin

  #termina de calcular B1
  max.B1 = length(x) * (1/length(x)) * log(1/length(x))
  B1 = sum(x * log(x))/max.B1

#calcular B2
#----------
x = model.map.standardized
x = x[!is.na(x)]
x = x/sum(x)
B2 = (1/sum(x^2) - 1) / (length(x)-1)

print(paste("B1 = ", round(B1, 4), sep=""))
print(paste("B2 = ", round(B2, 4), sep=""))

results = list(B1 = B1, B2 = B2)
return(results)

}





#Genera background restringido a una distancia alrededor de los puntos conocidos de presencia
restrictedBackground = function(presencia, variables, buffer.distance.km, sampling.percentage){

  #librerías requeridas
  require(raster)
  require(rgeos)

  #convertir la tabla de presencia en objeto "sp" (clase SpatialPoints)
  presencia.sp=presencia[ , c("x", "y")]
  coordinates(presencia.sp)=c("x", "y")

  #hacer buffer
  buffer=circles(presencia.sp, d=buffer.distance.km*1000, lonlat=TRUE)
  buffer.dissolve=gUnaryUnion(buffer@polygons)

  #extraemos los identificadores de las celdas
  celdas=unlist(cellFromPolygon(object=variables, p=buffer.dissolve))

  #tomamos uno de los mapas de variables como plantilla para la máscara
  mascara.temp=raster(variables[[1]])

  #creamos vector de valores
  valores=rep(NaN, ncell(mascara.temp))
  valores[celdas]=1

  #mascara a partir de valores
  mascara.temp=setValues(mascara.temp, values=valores)
  mascara.temp=mask(mascara.temp, mask=variables[[1]])

  #generamos background
  background = data.frame(randomPoints(mask=mascara.temp, n=floor(length(celdas)/(100/sampling.percentage))))

  #ploteamos todo para ver que ha salido bien (background en rojo)
  plot(variables[[1]])
  points(background$x, background$y, cex=0.1, col="red")

  #extraemos los valores de las variables sobre los puntos
  background.variables=data.frame(extract(variables, background))

  #unimos las coordenadas con los valores de las variables
  background=cbind(background, background.variables)

  #le añadimos la columna de presencia
  background$presencia=0

  #unimos con las presencias
  presencia.background=rbind(presencia[,colnames(presencia) %in% colnames(background)], background)

  return(presencia.background)
}






#############################################################################
#ReduceSpatialClustering
#This function reduces the spatial clustering of a set of presence records. It is intended to reduce spatial autocorrelation, and reduce sampling bias, specially at larger geographical scales.

#It requires two different arguments:
#data.table: a table with two fields representing latitude (named 'y') and longitude (named 'x')

#a minimum.distance value, provided in the same units as the coordinates, that will define the search radius when looking for pairs of coordinates within search distance to get rid of. Hint: the minimum distance can be extracted from the resolution of a raster containint the environmental factors, like "minimum.distance=xres(v.brick.20km)"
ReduceSpatialClustering = function(data, minimum.distance){

#count rows
row=1


#repite la operación hasta que se cumple la condición de salida
repeat{

  #contenido de la fila (para no tirar de toda la tabla en todas las operaciones)
  f=data[row, ]

  #genera los límites de la cuadrícula de búsqueda
  ymax=f$latitude + minimum.distance
  ymin=f$latitude - minimum.distance
  xmax=f$longitude + minimum.distance
  xmin=f$longitude - minimum.distance

  #selecciona de la tabla los datos con coordenadas dentro del rectángulo que no tienen las mismas coordenadas que la fila con la que estamos trabajando, y las elimina de la tabla
  data=data[!((data$latitude <= ymax) & (data$latitude >= ymin) & (data$longitude <= xmax) & (data$longitude >= xmin) & (data$latitude != f$latitude | data$longitude != f$longitude)), ]

  #estima de filas por procesar
  print(paste("Processed rows: ", row, " out of ", nrow(data), sep=""))

  #suma 1 al contador de la fila
  row=row+1

  #condición de salida cuando llega a la última fila
  if(row>=nrow(data))break
}

return(data)

}






#############################################################################
#WEIGHT PRESENCE/BACKGROUND DATA
WeightPresenceBackground=function(presence.column){

  #computing weight for presences
  n.presences=sum(presence.column)
  print(paste("Presence points = ", n.presences, sep=""))
  weight.presences=1/n.presences
  print(paste("Weight for presences = ", weight.presences, sep=""))

  n.background=length(presence.column)-n.presences
  print(paste("Background points = ", n.background, sep=""))
  weight.background=1/n.background
  print(paste("Weight for background = ", weight.background, sep=""))

  #generamos un vector con los los pesos
  weights=c(rep(weight.presences, n.presences), rep(weight.background, n.background))
  #return(weights)
}


#############################################################################
#http://modtools.wordpress.com/2013/08/14/dsquared/
# Linear models come with an R-squared value that measures the proportion of variation that the model accounts for. The R-squared is provided with summary(model) in R. For generalized linear models (GLMs), the equivalent is the amount of deviance accounted for (D-squared; Guisan & Zimmermann 2000), but this value is not normally provided with the model summary. The Dsquared function, now included in the modEvA package (Barbosa et al. 2014), calculates it. There is also an option to calculate the adjusted D-squared, which takes into account the number of observations and the number of predictors, thus allowing direct comparison among different models (Weisberg 1980, Guisan & Zimmermann 2000).
Dsquared = function(model, adjust = TRUE) {
  # version 1.1 (13 Aug 2013)
  # calculates the explained deviance of a GLM
  # model: a model object of class "glm"
  # adjust: logical, whether or not to use the adjusted deviance taking into acount the nr of observations and parameters (Weisberg 1980; Guisan & Zimmermann 2000)
  d2 = (model$null.deviance - model$deviance) / model$null.deviance
  if (adjust) {
    n = length(model$fitted.values)
    p = length(model$coefficients)
    d2 = 1 - ((n - 1) / (n - p)) * (1 - d2)
  }
  return(d2)
}  # end Dsquared function


#############################################################################
#WHITENING
WhiteningEnsemble = function(average, deviation, plot.points, points){

  #This code is derived from the one written by Tomislav Hengl (available here: http://spatial-analyst.net/wiki/index.php?title=Uncertainty_visualization). The main difference is that my version doesn't rely on a spatial dataframe, but on raster maps (library raster).

  #average: raster map representing the average of a few spatial models
  #deviation: raster map representing the standard deviation of a few spatial models
  #points = two columns in x - y order to plot point coordinates
  #name = name of the analysis
  #path, without the final slash

  #required libraries
  require(colorspace)
  require(plotrix)
  #   require(SGDF2PCT)
  require(rgdal)

  #stacking the layers together
  ensemble=stack(average, deviation)
  names(ensemble)=c("average", "deviation")

  #STRECH THE AVERAGE VALUES ONLY IF THE MAXIMUM VALUE OF THE AVERAGE IS HIGHER THAN 1
  #   if (max(as.vector(ensemble[["average"]]), na.rm=TRUE) > 1){
  ensemble[["average"]]=setValues(ensemble[["average"]], plotrix::rescale(as.vector(ensemble[["average"]]), c(0,1)))
  #   }

  #STRECH THE VALUES OF THE NORMALIZED DEVIATION
  ensemble[["deviation"]]=setValues(ensemble[["deviation"]], plotrix::rescale(as.vector(ensemble[["deviation"]]), c(0,1)))

  #DERIVE HUE
  H=-90-as.vector(ensemble[["average"]])*300
  H=ifelse(as.vector(H)<=-360, as.vector(H)+360, as.vector(H))
  H=ifelse(as.vector(H)>=0, as.vector(H), (as.vector(H)+360))

  #DERIVE SATURATION
  S=1-as.vector(ensemble[["deviation"]])
  V=0.5*(1+as.vector(ensemble[["deviation"]]))

  #CONVERT TO RGB
  RGB = as(HSV(H, S, V), "RGB")

  #CREATES THE RGB LAYERS
  R=setValues(ensemble[["deviation"]], as.integer(ifelse(is.na(as.vector(ensemble[["average"]])), 255, RGB@coords[,1]*255)))
  G=setValues(ensemble[["deviation"]], as.integer(ifelse(is.na(as.vector(ensemble[["average"]])), 255, RGB@coords[,2]*255)))
  B=setValues(ensemble[["deviation"]], as.integer(ifelse(is.na(as.vector(ensemble[["average"]])), 255, RGB@coords[,3]*255)))
  #stack
  RGB=stack(R,G,B)
  names(RGB)=c("R", "G", "B")

  #PLOTTING THE MAP
  layout(matrix(c(1,1,1,1,1, 2,2,2, 1,1,1,1,1, 2,3,2, 1,1,1,1,1, 2,2,2), nrow = 3, ncol = 8, byrow = TRUE))

  plotRGB(RGB, 1, 2, 3)

  plot(0,type='n',axes=FALSE,ann=FALSE)

  #LEGEND (taken from Tomislav Hengl's code as it)
  #########
  legend.2D = expand.grid(x=seq(.01,1,.01),y=seq(.01,1,.01))
  # Hues
  legend.2D$tmpf1 = -90-legend.2D$y*300
  legend.2D$tmpf2 = ifelse(legend.2D$tmpf1<=-360, legend.2D$tmpf1+360, legend.2D$tmpf1)
  legend.2D$H = ifelse(legend.2D$tmpf2>=0, legend.2D$tmpf2, (legend.2D$tmpf2+360))
  # Saturation
  legend.2D$S = 1-legend.2D$x
  # Intensity
  legend.2D$V = 0.5+legend.2D$x/2

  gridded(legend.2D) = ~x+y
  legend.2D = as(legend.2D, "SpatialGridDataFrame")
  spplot(legend.2D["H"], col.regions=rev(gray(0:20/20)))
  spplot(legend.2D["S"], col.regions=rev(gray(0:20/20)))
  spplot(legend.2D["V"], col.regions=rev(gray(0:20/20)))

  legendimg = as(HSV(legend.2D$H, legend.2D$S, legend.2D$V), "RGB")
  #   plot(legendimg)
  legend.2D$red = as.integer(legendimg@coords[,1]*255)
  legend.2D$green = as.integer(legendimg@coords[,2]*255)
  legend.2D$blue = as.integer(legendimg@coords[,3]*255)

  #Display as a RGB image:
  legend.2Dimg = SGDF2PCT(legend.2D[c("red", "green", "blue")], ncolors=256, adjust.bands=FALSE)
  legend.2D$idx = legend.2Dimg$idx

  #
  image(legend.2D, "idx", col=legend.2Dimg$ct, main="Legend")
  axis(side=2, at=c(0, 0.25, 0.50, 0.75, 1), line=0, lwd=2)
  axis(side=1, at=c(0, 0.25, 0.50, 0.75, 1), line=-4, lwd=2)
  mtext("Habitat suitability", side=2, line=3, cex=1.5)
  mtext("Standard deviation", side=1, line=-1, cex=1.5)

}


