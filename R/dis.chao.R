dis.chao <- function(comm, index="jaccard", version="rare", freq=NULL){

indices <- c("jaccard", "sorensen")
index   <- indices[pmatch(index, indices)]    

versions <- c("probability", "rare")
version  <- versions[pmatch(version, versions)]   

sites0 <- nrow(comm)
S0 <- ncol(comm)  
comm <- comm[rowSums(comm)>0, colSums(comm)>0]
sites <- nrow(comm)
S <- ncol(comm)
if(sites0 < sites){warning("There is at least one sample with no species.")}
if(S0 < S){warning("There is at least one species with no record in samples. It/they was/were removed")}


if(version=="probability" & is.null(freq)==FALSE){
   warning("The probability version of the Chao index is not available for incidence-based data. Argument 'freq' was ignored and thus the abundance-based formulae was used.")
}


resu <- matrix(NA, sites, sites)

for(i in 1:(sites-1)){
   samp1 <- comm[i, ]
   for(j in (i+1):sites){
      samp2 <- comm[j, ]
      pair <- rbind(samp1, samp2)
      pair.pa <- ifelse(pair > 0, 1, 0)
      D12.TF  <- colSums(pair.pa) == 2 
      
      if(all(D12.TF == F)) {resu[j, i] <- 0} # in case there is no shared species, simi=0 
         else{
             D12.spp <- pair[ , D12.TF]
             D12.spp <- as.matrix(D12.spp) #because if a single sp is shared, object is a vector. 
             n <- sum(pair[1,])
             m <- sum(pair[2,])
                
             if(version == "probability"){
                U <- sum(D12.spp[1, ]/n) 
                V <- sum(D12.spp[2, ]/m) 
             }
             
             if(version == "rare"){
                f1. <- sum(D12.spp[1, ]==1)
                f.1 <- sum(D12.spp[2, ]==1)
              
                f2. <- sum(D12.spp[1, ]==2)
                if(f2.==0) {f2. <- 1} # because it is used in the denominator.
                f.2 <- sum(D12.spp[2, ]==2)
                if(f.2==0) {f.2 <- 1}
              
                U1 <- sum(D12.spp[1,]/n) 
                V1 <- sum(D12.spp[2,]/m) 

                U3 <- f.1/(2*f.2)
                V3 <- f1./(2*f2.)

                U4 <- sum( (D12.spp[1,]/n) * (D12.spp[2,] == 1)  )
                V4 <- sum( (D12.spp[2,]/m) * (D12.spp[1,] == 1)  )


                if(is.null(freq)){          # for abundance data
                   U2 <- ((m-1)/m)  
                   V2 <- ((n-1)/n)                  
                } 

                if(is.null(freq) == FALSE){  # for incidence data
                   if(is.vector(freq)==FALSE){stop("Argument 'freq' must be a numeric vector")}
                   if(length(freq)!=nrow(comm)){stop("Number of samples must be the same as length of 'freq'")}
                  
                   w <- freq[i]
                   z <- freq[j]
                  
                   U2 <- (z-1)/z                  
                   V2 <- (w-1)/w
                } 
                
                
                U  <- U1 + U2*U3*U4
                V  <- V1 + V2*V3*V4

                U <- ifelse(U>1, 1, U)
                V <- ifelse(V>1, 1, V)
             } # closes if(version=="shared")
        
             
             if (index == "jaccard") {resu[j, i] <- (U*V)/(U+V-U*V)}
             if (index == "sorensen"){resu[j, i] <- (2*(U*V))/(U+V)}
                  
           } # closes else (when there is/are shared species)      
    } # closes for j
} # closes for i
rownames(resu) <- rownames(comm)

return(1-(as.dist(resu)))  
}



