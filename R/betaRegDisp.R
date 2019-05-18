betaRegDisp <- function(y, x, xy.coords = NULL, ws = 3, 
                        method.1 = "jaccard",
                        method.2 = "ruzicka", 
                        method.3 = "ruzicka",
                        independent.data = FALSE,
                        illust.plot = FALSE){

  y <- y[order(x, decreasing = FALSE), ]
  if(!is.null(xy.coords))xy.coords <- xy.coords[order(x, decreasing = FALSE), ]
  x <- x[order(x, decreasing = FALSE)]
  N <- (ws)/2
  
  is.even <- function(ee){ ee %% 2 == 0 } 
    if(is.even(ws)==FALSE){
    N <- N - 0.5
  }
  
  size <- length(x)
  SEQ <- 1:(size-ws+1)
  if(is.even(ws)==FALSE){
    SEQ <- (N+1):(size-N)
    }
  
  if (independent.data) SEQ <- seq (1,size-(ws-1),ws)
  if (independent.data && is.even(ws)==FALSE) SEQ <- seq (N+1,size-N,ws)
  
  n <- length(SEQ)
  #print(n)
  result <- matrix(0, n, 10)
  colnames(result) <- c("grad", "mean.grad",
                        "mean.diss.pairs", "mean.diss.focal", "mean.dist.cent",
                        "SS.group", "SS.focal", 
                        "beta.TOT","beta.NES","beta.TUR")
  result[, "grad"] <- x[SEQ]
  rownames(result) <- names(x)[SEQ]
  disT <- vegdist(y, method = method.1)
  if(!is.null(xy.coords)){
     geo.dist<-matrix(0,n,2,
                      dimnames = list(NULL,c('mean.geodist','focal.geodist')))
  }
  
  
  count = 1
  for(i in SEQ){
    group <- rep("B", times = size)
    
    if(is.even(ws)==FALSE){
    sites <- (i-N):(i+N)
    }
    
    if(is.even(ws)==TRUE){
      sites <- i:(i+ws-1)
    }
    
    sites <- c(i, sites[sites!=i])
    
    ### Mean environmental gradient
    result[count, "mean.grad"] <- mean(x[sites]) 
    
    ### Mean euclidean (geo) distance
    if(!is.null(xy.coords)){
      d <- dist(xy.coords[sites,], method = "euclidean")####[1:ws]
      geo.dist[count,] <- c(mean(d), mean(d[1:(ws-1)]))
    }
    
    ### Illustrative plot
    if(illust.plot == TRUE){
      plot(1:size, x)
      points(sites, x[sites], cex = 2, col = "blue")
      points(i, x[i], cex = 3, col = "red", pch = 0)
    }
    
    
    group[sites] <- "A"
    mat <- y[sites, ]
    
    ### Mean pair-wise and focal dissimilarity
    dis <- vegdist(mat, method = method.1)
    result[count, "mean.diss.pairs"] <- mean(dis) #Mean for all sites
    result[count, "mean.diss.focal"] <- mean(dis[1:(ws-1)]) #Mean focal site
    
    
    ### Distance from group centroid
    mod <- betadisper(disT, group = group)
    d <- mod$distances
    result[count, "mean.dist.cent"] <- mean(d[group=="A"]) 
     
    
    ### Multisample beta diversity of Baselga (2010, 2017)
    beta.b <- beta.multi.abund(x = mat, index.family = method.3)
    nomes.m3 <- c("bray", "ruzicka")
    method.3 <- nomes.m3[pmatch(method.3, nomes.m3)]
    if(method.3 == "bray"){
      result[count, "beta.TOT"] <- beta.b$beta.BRAY
      result[count, "beta.TUR"] <- beta.b$beta.BRAY.BAL
      result[count, "beta.NES"] <- beta.b$beta.BRAY.GRA
  }
   if(method.3 == "ruzicka"){
      result[count, "beta.TOT"] <- beta.b$beta.RUZ
      result[count, "beta.TUR"] <- beta.b$beta.RUZ.BAL
      result[count, "beta.NES"] <- beta.b$beta.RUZ.GRA
   }
  
  
    ### Sum of Squares (SS) of Legendre and De Caceres (2013)
      res.SS <- beta.div(Y = mat, method = method.2, nperm = 0)
      result[count, "SS.group"] <- res.SS$beta["SStotal"]
      result[count, "SS.focal"] <- res.SS$LCBD[rownames(mat)[1]]*res.SS$beta["SStotal"]
      
      
    #cat (paste(n-count, 'points left to compute\n',sep=' '))
    count = count + 1
  }
  if(!is.null(xy.coords))result<-cbind(result,geo.dist)
  return(result)
}