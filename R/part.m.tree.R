part.m.tree <- function(comm, tree, index.family="sorensen"){
if (is.null(tree$edge.length)) {stop("Tree has no branch lengths, cannot compute pd")}
if (!is.rooted(tree))   {stop("Rooted geny required for phylosor calculation")}

families<-c("sorensen","jaccard")
index<-families[pmatch(index.family, families)]
s<-nrow(comm)

#********************************************************* SHARED matrix
samp.comb <- matrix(NA, s * (s - 1)/2, ncol(comm)) #
colnames(samp.comb) <- colnames(comm)
i <- 1
   for (l in 1:(s - 1)) {
      for (k in (l + 1):s) {
         samp.comb[i, ] <- colSums(comm[c(l,k), ]) #1 sample will all spp in samples l and k.
         i <- i + 1
        }
    }
    pdsamp <- pd(comm, tree)
    pdsamp.comb <- pd(samp.comb, tree)
    
SHARED<-matrix(0, s, s)
diag(SHARED)<-pdsamp[,"PD"]
i <- 1
for (l in 1:(s - 1)) {
   pdl <- pdsamp[l, "PD"]   
   for (k in (l + 1):s) {
       pdk <- pdsamp[k, "PD"]
       pdcomb <- pdsamp.comb[i, "PD"]
       SHARED[k,l]<- pdl + pdk - pdcomb
       SHARED[l,k]<- pdl + pdk - pdcomb
       i <- i + 1
    }
}
#********************************************************* End of SHARED matrix.

#********************************************************* NOT SHARED.
not.shared <-  abs(sweep(SHARED, 2, diag(SHARED)))
sum.not.shared <- not.shared + t(not.shared)
max.not.shared <- pmax(not.shared, t(not.shared))
min.not.shared <- pmin(not.shared, t(not.shared))
    
                
sumSi <- sum(diag(SHARED)) # species by site richness     --> SUM OF PDs
St<-pd(rbind(colSums(comm),colSums(comm)), tree)[1,1] # BECAUSE pd ONLY WORKS ON 2 DIMENSIONS
                                                      # St = TOTAL PD IN THE POOLED SAMPLE
#St <- sum(colSums(comm) > 0)  # regional species richness 

a <- sumSi - St            # multi site shared species term
maxbibj <- sum(max.not.shared[lower.tri(max.not.shared)])
minbibj <- sum(min.not.shared[lower.tri(min.not.shared)])
#********************************************************* End of NOT SHARED.


#********************************************************* Indices and resu.

if(index=="sorensen"){
  beta.turn <- minbibj / (minbibj + a)
  beta.nest <- (a / (minbibj + a)) * ((maxbibj - minbibj) / ((2 * a) + maxbibj + minbibj))
  beta <- (minbibj + maxbibj) / (minbibj + maxbibj + (2 * a))
}

if(index=="jaccard"){
  beta.turn <- (2*minbibj) / (2*minbibj + a)
  beta.nest <- (a / (2*minbibj + a)) * ((maxbibj - minbibj) / (a + maxbibj + minbibj))
  beta <- (minbibj + maxbibj) / (minbibj + maxbibj +  a)
}


if(index=="sorensen"){ resu<-list(n.sites=s, SOR=beta, SIM=beta.turn, SNE=beta.nest)}
if(index=="jaccard") { resu<-list(n.sites=s, JAC=beta, JTU=beta.turn, JNE=beta.nest)}

return(resu)}
