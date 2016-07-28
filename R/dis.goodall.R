dis.goodall <- function(comm, p.simi="steinhaus", approach="proportion"){ 
   
if(any(colSums(comm)==0)){warning("There is at least one species with no record in samples. It/they was/were removed.")}
   
comm <- comm[rowSums(comm)>0, colSums(comm)>0]

similarities <- c("gower","steinhaus")
simi         <- similarities[ pmatch(p.simi, similarities) ]

app   <- c("proportion", "chisquare")
approach <- app[ pmatch(approach, app)]
  
samples <- nrow(comm)
spp     <- ncol(comm)
n.pairs <- (samples*(samples-1))/2 # number of columns in table (a) of L&L.
index   <- combn(samples, 2) # 2-rows matrix with addresses of sample pairs. 


#                  >>>>> Legendre & Legendre Step (a) <<<<<
partial <- matrix(NA, spp, n.pairs) # matrix simi partial. 1 row for each sp. 
                                    # 1 column for each sample pair.
ranges <- apply(comm, 2, range)
R      <- abs( apply(ranges, 2, diff) )

for(j in 1:spp){
   for(b in 1:n.pairs){
      pair <- index[, b]
      
      if(all(comm[pair, j] == 0)) {partial[j, b] <- 0}
   	  else{
         
         if(simi=="gower"){ # This is s12j described in p259 of L&L.
            num <- abs(diff(comm[pair, j]))
            partial[j, b] <- 1-(num/R[j])
         }
         
         if(simi=="steinhaus"){ # This is BC similarity for 1 sp. S17 of L&L.
            partial[j, b] <- (2*min(comm[pair, j])) / sum(comm[pair, j])
         }
         
      } # closes else	         
   }
} 


#                  >>>>> Legendre & Legendre Step (b) <<<<<
prop.partial <- matrix(NA, spp, n.pairs)  # prop. values >= focal value for sp j.
for(j in 1:spp){
	 for(b in 1:n.pairs){ 
      prop.partial[j, b] <- sum(partial[j, ] >= partial[j, b]) / n.pairs
   } 
}  


#                  >>>>> Legendre & Legendre Step (c) <<<<<
# Notice L&L presents it as a dist (or matrix) object.
# However, it is best for computation to maintain the vector object.
product <- apply(prop.partial, 2, prod)


#                  >>>>> Legendre & Legendre Step (d and e) <<<<<
if(approach=="proportion"){ # Proportions of values that are >= focal values
  resu <- numeric(n.pairs)
  for(i in 1:n.pairs){
     resu[i] <- sum(product >= product[i])/n.pairs
  } 
}  
if(approach=="chisquare"){
   chi <- -2 * log(product)
   resu <- pchisq(chi, 2*spp) # Notice this already is the complement of p. 
}


#   >>>>> Transforming into dist object and dissimilarities <<<<<             
good <- matrix(NA, samples, samples)
rownames(good) <- rownames(comm)
for(i in 1:n.pairs){ good[index[2, i], index[1, i]] <- resu[i]  }

return(1-as.dist(good))  
}


