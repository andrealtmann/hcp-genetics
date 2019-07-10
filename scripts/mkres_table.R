
mk_vtk_table <- function(results, maxint=5, hem="Left"){
  fs.map <- read.csv("../data/FS_mapping.csv")
  rownames(fs.map) <- tolower(fs.map[,"FS_abbrev"])
  fs.map.hem <- subset(fs.map, Hemisphere==hem)

  usco <- results[rownames(fs.map.hem)]

  vtkcode <- fs.map.hem$code - 1000
  if (hem == "Right")
    vtkcode <- fs.map.hem$code - 2000

  fs.map.hem <- data.frame(fs.map.hem, vtkcode)

  fs.map.out <- data.frame(fs.map.hem, SCORE=usco[rownames(fs.map.hem)])

  omap <- matrix(NA,nrow=40, ncol=4)
  omap[,1] <- 1:40
  rownames(omap) <- 1:40
  
  cp <- colorRampPalette(c("red","yellow"))(100 * maxint)
  dummy <- floor(fs.map.out$SCORE*1000)
  dummy[dummy <= 0] <- 1
  dummy[dummy >100 * maxint] <- 100 * maxint
  hex.col <- cp[dummy]
  rgb.col <- data.frame(t(col2rgb(hex.col))) 
  rownames(rgb.col) <- fs.map.out$vtkcode

  omap[rownames(rgb.col),2] <- rgb.col[,1]
  omap[rownames(rgb.col),3] <- rgb.col[,2]
  omap[rownames(rgb.col),4] <- rgb.col[,3]
  colnames(omap) <- c("label","red","green","blue")
  omap["40",2:4] <- rep(192,3)
  omap[is.na(omap)] <- 0

  #gray out non sign results
  idx <- paste(fs.map.out[fs.map.out$SCORE<0.05,"vtkcode"])
  omap[idx,2] <- 192
  omap[idx,3] <- 192
  omap[idx,4] <- 192

  return(omap)
}