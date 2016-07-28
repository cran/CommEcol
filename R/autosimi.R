autosimi <- function(comm, method="bray", binary=FALSE, log.transf=FALSE, 
                     simi=TRUE, permutations=50){

n.samp <- nrow(comm)
spp    <- ncol(comm)
size   <- n.samp%/%2
sizes  <- 1:size

resu.m <- matrix(NA, size, permutations)

for(i in 1:permutations){
   for(j in sizes){
      if(j==1){ pair<-comm[sample(1:n.samp,2), ]}
      if(j>1) {
          temp<-comm[sample(1:n.samp, j*2),]
          temp1<-temp[1:j,]
          temp2<-temp[(j+1):(j*2),]
          temp1<-colSums(temp1)
          temp2<-colSums(temp2)
          pair<-rbind(temp1,temp2)
          }

      if(binary==TRUE){pair<-ifelse(pair>0,1,0)}# close if binary
      if(binary==FALSE & log.transf==TRUE){pair<-log(pair+1)}
      if(binary==TRUE & log.transf==TRUE){stop("You can not log-transform presence/absence data")}
      
      if(simi==TRUE){resu.m[j,i]<-as.numeric(1-vegdist(pair, method, binary))}
      if(simi==FALSE){resu.m[j,i]<-as.numeric(vegdist(pair, method, binary))}
   }#close for j
}# close for i

mean.perm <- rowMeans(resu.m)
resu <- data.frame(sizes, mean.perm)
indexes <- c("manhattan", "euclidean", "canberra",      "bray", "kulczynski", "jaccard", "gower",
              "altGower",  "morisita",     "horn", "mountford", "raup" , "binomial", "chao")
colnames(resu) <- c("sample.size", indexes[pmatch(method,indexes)])
return(resu)
}


