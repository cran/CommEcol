dis.nness <- function(comm, m=NULL, ness=FALSE){

comm   <- as.matrix(comm)
sites  <- nrow(comm)
S      <- ncol(comm)

comm.redu   <- comm[rowSums(comm)>0, colSums(comm)>0]
sites.redu <- nrow(comm.redu)
S.redu     <- ncol(comm.redu)  

if(sites > sites.redu){
   warning("There is at least one sample with no species. Dissimilarities were not computed for pairs including it/them.")}
if(S > S.redu){
   warning("There is at least one species with no record in samples. It/they were removed.")}

if(is.null(m)){ m <- dis.nness.find.m(comm, ness=ness) }

if(abs(m - round(m)) > .Machine$double.eps^0.5){
   stop("m must be a non-null integer.")}

comm <- comm[, colSums(comm)>0] # Notice empty sites are maintained.
S <- S.redu

if(any(abs(comm - round(comm, digits=0)) > .Machine$double.eps^0.5)){
   warning("There is/are non-integer(s) in the dataset. NNESS and NESS are most suitable to abundance data.")
   if(any(rowSums(comm)>1000)){
      comm <- ceiling(comm)
      warning("At least one sample contains more than 1000 individuals. Non-integer abundance values were rounded up to the next upper integer.")   
   }
}

if(ness==FALSE & any(m    > rowSums(comm))){
   warning("There is at least one sample in which total abundance is lower than m. NNESS dissimilarities were not computed for pairs including it/they.")}   

if(ness==TRUE  & any((m*2)> rowSums(comm))){
   warning("There is at least one sample in which total abundance is lower than m*2. NESS dissimilarities were not computed for pairs including it/they.")}   

if(any(rowSums(comm)>1000)){
   choo <- function(n, k)    {gmp::chooseZ(n, k)}
   nu <- function(x)         {gmp::asNumeric(x) }
}
   else{
      choo <- function(n, k) {choose(n, k) }
      nu <- function(x)      {as.numeric(x)}
   }



resu <- matrix(NA, sites, sites)
for(i in 1:(sites-1)){
   samp1 <- comm[i, ]
   N1 <- sum(samp1)
   N1.m <- choo(N1, m)
      
   for(j in (i+1):sites){
      samp2 <- comm[j, ]
      N2 <- sum(samp2)
      N2.m <- choo(N2, m)
      
      if(N1.m!=0 & N2.m!=0){
      
         ESS12.S <- numeric(S)
         ESS1.S  <- numeric(S)
         ESS2.S  <- numeric(S)
         for(k in 1:S){
            nume1 <- choo(nu(N1 - samp1[k]), m)
            nume2 <- choo(nu(N2 - samp2[k]), m)
            ESS12.S[k] <- nu((1-(nume1 / N1.m)) * (1-(nume2 / N2.m))) ## eq.4 of Grassle and Smith
         
            if(ness==TRUE){
               nume1.2m <- choo((N1 - samp1[k]), 2*m)
               nume2.2m <- choo((N2 - samp2[k]), 2*m)
               N1.2m <- choo(N1, 2*m)
               N2.2m <- choo(N2, 2*m)            
               ESS1.S[k] <- nu( 1-2*(nume1 / N1.m) + (nume1.2m / N1.2m)) # eq.5 of Grassle and Smith
               ESS2.S[k] <- nu( 1-2*(nume2 / N2.m) + (nume2.2m / N2.2m))
               # the 1-2* part above must not include parenthesis. 
               #        In fact, 1-(2*...)+() gives wrong results.
      
            
            } else{
               ESS1.S[k]  <- nu((1-(nume1 / N1.m))^2)
               ESS2.S[k]  <- nu((1-(nume2 / N2.m))^2)     
            }  
         } #closes for k
      
         ESS12 <- sum(ESS12.S)
         ESS1  <- sum(ESS1.S)
         ESS2  <- sum(ESS2.S)
         resu[j, i] <- (2*ESS12) / (ESS1 + ESS2)
      }# closes if(N1.m!=0 & N2.m!=0){
   } # closes for j
} # closes for i
      
resu <- ifelse(resu>1, 1, resu)
colnames(resu) <- rownames(comm)
rownames(resu) <- rownames(comm)

resu.dist <- 1-(as.dist(resu))
return(resu.dist)
}
