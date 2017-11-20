standExtent <- function(d.large, d.small, ini.large=NULL, ini.small=NULL){
d.large <- as.matrix(d.large)
d.small <- as.matrix(d.small)
n.large <- nrow(d.large)
n.small <- nrow(d.small)

if(is.null(rownames(d.large))){
   rownames(d.large) <- paste("l", 1:n.large, sep="")
}
if(is.null(rownames(d.small))){
   rownames(d.small) <- paste("s", 1:n.small, sep="")
}
   
names.large <- rownames(d.large)
names.small <- rownames(d.small)
   
if(is.null(ini.large)){
   ini.large <- sample(names.large, 1)
}
if(is.null(ini.small)){
   ini.small <- sample(names.small, 1)
}
   
resu <- matrix(NA, nrow=n.large, ncol=5)
resu <- as.data.frame(resu)
colnames(resu) <- c("site.large", "dist.mean.large", 
                    "site.small", "dist.mean.small", "percent.dif")
   
resu[1, "site.large"] <- ini.large
resu[1, "site.small"] <- ini.small
resu[1, "dist.mean.large"] <- 0
resu[1, "dist.mean.small"] <- 0
   
ini2.large <- sort(d.large[ini.large, ])[2] # the 1st one is the ini-ini = 0.
resu[2, "site.large"] <- names(ini2.large)
resu[2, "dist.mean.large"] <- ini2.large

i1s.pos <- which(rownames(d.small)==ini.small)
ini2.small <- sort(abs(d.small[i1s.pos, -i1s.pos] - ini2.large))[1] #the ini.small col is removed
resu[2, "site.small"] <- names(ini2.small)
resu[2, "dist.mean.small"] <- d.small[ini.small, names(ini2.small)]
   

#                      ***** LARGE *****
for(i in 3:n.large){
   set.large   <- resu[1:(i-1), "site.large"] # the set of sites selected so far.
   names.large <- setdiff(names.large, set.large) # remove names present in the set.
   
      #  Obtain average distances among all members of a new LARGE set
      #  formed by each non-set site (site la) and the previous set.
   temp.large <- matrix(, nrow=(n.large-(i-1)), ncol=2) # to store temp dists
   colnames(temp.large) <- c("site.large", "dist.mean.large")
   temp.large <- as.data.frame(temp.large)
   count <- 1
   for(la in names.large){
      temp.set.large <- c(set.large, la) # the selected set so far and the "test" site.
      temp.d.large  <- d.large[temp.set.large, temp.set.large]
      temp.large[count, "site.large"] <- la
      temp.large[count, "dist.mean.large"] <- sum(temp.d.large) / (i*(i-1))
      count <- count+1
   } # closes for la. 
      
   sel.large <- which.min(temp.large[, "dist.mean.large"])
   resu[i, "site.large"]      <- temp.large[sel.large, "site.large"]
   resu[i, "dist.mean.large"] <- temp.large[sel.large, "dist.mean.large"]
} # closes for i. This will build the final result for the LARGE set.
   
   

#                      ***** SMALL *****
if(n.large >= n.small){
   n.aglo <- n.small
}
   else{
      n.aglo <- n.large
   }

for(j in 3:n.aglo){
   set.small   <- resu[1:(j-1), "site.small"]
   names.small <- setdiff(names.small, set.small)
   
      #  Obtain average distances among all members of a new SMALL set
      #  formed by each non-set site (site sm) and the previous set.
   temp.small <- matrix(, nrow=(n.small-(j-1)), ncol=2)
   colnames(temp.small) <- c("site.small", "dist.mean.small")
   temp.small <- as.data.frame(temp.small)
   count <- 1
   for(sm in names.small){
      temp.set.small <- c(set.small, sm) # the selected set so far and the "test" site.
      temp.d.small   <- d.small[temp.set.small, temp.set.small]
      temp.small[count, "site.small"] <- sm
      temp.small[count, "dist.mean.small"] <- sum(temp.d.small) / (j*(j-1))
      count <- count+1
   } # closes for sm. 
   
   sel.small <- which.min(abs(temp.small[ , "dist.mean.small"] - 
                                    resu[j, "dist.mean.large"]  ))
   resu[j, "site.small"]      <- temp.small[sel.small, "site.small"]
   resu[j, "dist.mean.small"] <- temp.small[sel.small, "dist.mean.small"]
} # closes for j. This will build the final result for the SMALL set.


resu[, "percent.dif"] <- round(((resu[, "dist.mean.large"] - 
                          resu[, "dist.mean.small"] ) /
                          resu[, "dist.mean.large"])*100, 5)
resu[1, "percent.dif"] <- 0
return(resu)
}