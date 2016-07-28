part.p.tree<-function (comm, tree, index.family="sorensen") {
if (is.null(tree$edge.length)) {stop("Tree has no branch lengths, cannot compute pd")}
if (!is.rooted(tree))   {stop("Rooted tree required for phylosor calculation")}

families<-c("sorensen","jaccard")
index<-families[pmatch(index.family, families)]

comm <- as.matrix(comm)
s <- nrow(comm)
   
phylodist<-matrix(NA, s, s)            ## To store dist. results
rownames(phylodist) <- rownames(comm)
colnames(phylodist) <- rownames(comm)

phylodist.turn<-phylodist
phylodist.nest<-phylodist
       
samp_comb <- matrix(NA, s * (s - 1)/2, ncol(comm))
colnames(samp_comb) <- colnames(comm)
   
i <- 1
for (l in 1:(s - 1)) {
   for (k in (l + 1):s) {
      samp_comb[i, ] <- comm[l, ] + comm[k, ] #1 sample including all spp. of samples l and k.
      i <- i + 1
   }
}

pdsamp <- pd(comm, tree)
pdsamp_comb <- pd(samp_comb, tree)
i <- 1

for (l in 1:(s - 1)) {
   pdl <- pdsamp[l, "PD"]
   for (k in (l + 1):s) {
      pdk <- pdsamp[k, "PD"]
      pdcomb <- pdsamp_comb[i, "PD"]
      pdsharedlk <- pdl + pdk - pdcomb
            
      a<-pdsharedlk
      b<-pdl-pdsharedlk
      c<-pdk-pdsharedlk

      if(index=="sorensen"){
         phylodist[k, l] <-     (b+c) / (2*a+b+c) # changed by ASM.
         phylodist.turn[k, l]<- (min(b,c)) / (a + min(b,c) )
         phylodist.nest[k, l]<- ((max(b,c)-min(b,c))/(2*a+b+c)) * (a/(a+min(b,c)))
      }

      if(index=="jaccard"){
         phylodist[k, l] <-     (b+c) / (a+b+c) 
         phylodist.turn[k, l]<- (2*min(b,c) ) / (a+2*min(b+c) )     
         phylodist.nest[k, l]<- ((max(b,c)-min(b,c))/(a+b+c)) * (a/(a+2*min(b,c)))
      }

      i <- i + 1     
   }
}
#phylodist.nest<-phylodist-phylodist.simp

if(index=="sorensen"){
resu<-list(sor=as.dist(phylodist), sim=as.dist(phylodist.turn), sne=as.dist(phylodist.nest))}

if(index=="jaccard"){
resu<-list(jac=as.dist(phylodist), jtu=as.dist(phylodist.turn), jne=as.dist(phylodist.nest))}

return(resu)
}
