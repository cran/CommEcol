## tree.nodf.one ==>> internal function ##
treeNodfOne<-function(comm, tree, order){
sites <- nrow(comm)
spp   <- ncol(comm)
tree  <- as.phylo(tree)
pds <- pd(comm, tree)$PD

if(order==TRUE){ 
   ord  <- order(pds, decreasing=TRUE)
   comm <- comm[ord, ]
   pds  <- pd(comm, tree)$PD
}
  
tree.nodf.dist <- matrix(NA, sites, sites)        ## To store dist. results
rownames(tree.nodf.dist) <- rownames(comm)
colnames(tree.nodf.dist) <- rownames(comm)

s.fraction.dist<-tree.nodf.dist

for(i in 1:(sites - 1)) {
   pd.set <- pds[i]

   for(j in (i + 1):sites) {
      pd.subset <- pds[j]

      if(pd.subset >= pd.set) {
         one.tree.nodf.par  <- 0
         one.s.fraction.par <- 0
      } else{
           comb <- colSums(comm[c(i,j), ])
           s.shared <- sum(comb > 1)
           one.s.fraction.par <- s.shared / sum(comm[j, ])
           
           comb <- ifelse(comb > 0, 1, 0)
           comb <- matrix(comb, 1)
           colnames(comb) <- colnames(comm)
           pd.comb <- pd(comb, tree)$PD
           pd.shared <- pd.set + pd.subset - pd.comb
           one.tree.nodf.par <- pd.shared/pd.subset             
         }#closes else
   
   tree.nodf.dist[j, i]  <- one.tree.nodf.par
   s.fraction.dist[j, i] <- one.s.fraction.par
   }# closes for j
}# closes for i

tree.nodf.dist  <- as.dist(tree.nodf.dist * 100)
s.fraction.dist <- as.dist(s.fraction.dist * 100)
return(list(tree.nodf.dist=tree.nodf.dist, s.fraction.dist=s.fraction.dist))
}

