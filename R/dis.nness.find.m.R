dis.nness.find.m <- function(comm, ness=FALSE){
   
min1 <- min(rowSums(comm))
if(ness==TRUE){
   if(min1 < 61){
      ms <- 1:(min1/2)
   } else{
      ms <- round(seq(1, (min1/2), length.out=30), digits=0)
   }   
}

if(ness==FALSE){
   if(min1 < 31){
      ms <- 1:min1
   } else{
      ms <- round(seq(1, min1, length.out=30), digits=0)
   }   
}


comp    <- (nrow(comm)*(nrow(comm)-1))/2
ms.n    <- length(ms)
dists.m <- matrix(, comp, ms.n)

for(i in 1:ms.n){
   dists.m[, i] <- as.vector(dis.nness(comm, m=ms[i], ness=ness))
	}
#print(dists.m)
kendall.resu <- cor(dists.m, method="kendall")
rownames(kendall.resu) <- ms
colnames(kendall.resu) <- ms
#print(kendall.resu)

dif1 <- kendall.resu[, 1] - kendall.resu[, ms.n] #correlations to m=1 and to the highest m.
dif2 <- abs(dif1) 
posi <- which.min(dif2) # the position of the m-value which makes cor to m=1 and to m=max most similar (minimal difference).
m <- ms[posi]
return(m)
}
