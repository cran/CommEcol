treeNodfTest<-function(comm, col.tree, order.rows=FALSE, 
                             row.tree, order.cols=FALSE, 
                             null.model="perm.rows", permutations=999){
   
models<-c(    "perm.rows",     "perm.cols",     "perm.rc", 
          "perm.tip.rows", "perm.tip.cols", "perm.tip.rc", "ff")

if (is.na(pmatch(null.model, models))) {
   stop("You should specify a null model.")
}

sites <- nrow(comm)
spp   <- ncol(comm)

obs <- treeNodf(comm, col.tree, order.rows, row.tree, order.cols)

rows.aleats           <- matrix(NA, nrow=1, ncol=3)
colnames(rows.aleats) <- c('tree.nodf', 's.fraction', 'topo.nodf')
cols.aleats           <- rows.aleats
mat.aleats            <- rows.aleats


if(missing(col.tree)==FALSE){
  rows.aleats           <- matrix(NA, nrow=permutations, ncol=3)
  colnames(rows.aleats) <- c('tree.nodf', 's.fraction', 'topo.nodf')
  col.tree.aleat <- col.tree
}

if(missing(row.tree)==FALSE){
   cols.aleats           <- matrix(NA, nrow=permutations, ncol=3)
   colnames(cols.aleats) <- c('tree.nodf', 's.fraction', 'topo.nodf')
   row.tree.aleat <- row.tree
}

if(missing(col.tree)==FALSE & missing(row.tree)==FALSE) {
   mat.aleats           <- matrix(NA, nrow=permutations, ncol=3)
   colnames(mat.aleats) <- c('tree.nodf', 's.fraction', 'topo.nodf')
}


comm.aleat <- comm

for(i in 1:permutations){
   
## Null models ##
   if(null.model=="perm.rows") {
      comm.aleat <- comm[sample(1:sites), ] 
   }
   
   if(null.model=="perm.cols") {
      comm.aleat <- comm[, sample(1:spp)]
   }
   
   if(null.model=="perm.rc")   {
      comm.aleat <- comm[sample(1:sites), ]
      comm.aleat <- comm.aleat[, sample(1:spp)]
   }
   
   if(null.model=="ff")        {
      comm.aleat <- permatswap(comm, method="quasiswap",
                               mtype="prab", times=1)$perm[[1]]
   }

   if(null.model=="perm.tip.cols"){
     col.tree.aleat$tip.label <- col.tree$tip.label[sample(length(col.tree$tip.label))]
   } 
  
   if(null.model=="perm.tip.rows"){
     row.tree.aleat$tip.label <- row.tree$tip.label[sample(length(row.tree$tip.label))]
   } 

   if(null.model=="perm.tip.rc"){
     col.tree.aleat$tip.label <- col.tree$tip.label[sample(length(col.tree$tip.label))]
     row.tree.aleat$tip.label <- row.tree$tip.label[sample(length(row.tree$tip.label))]} 


## Calculations of statistics ##
   if(missing(col.tree)==FALSE & missing(row.tree)==FALSE) {
     resu.m<- treeNodf(comm=comm.aleat,
                       col.tree=col.tree.aleat, order.rows=FALSE,
                       row.tree=row.tree.aleat, order.cols=FALSE)
     rows.aleats[i, ] <- resu.m$rows
     cols.aleats[i, ] <- resu.m$cols
     mat.aleats[i, ]  <- resu.m$mat     
   } else {
     if(missing(col.tree)==FALSE){
       rows.aleats[i, ] <- treeNodf(comm=comm.aleat, col.tree=col.tree.aleat, 
                                    order.rows=FALSE)$rows
     }
     
     if(missing(row.tree)==FALSE){
       cols.aleats[i, ] <- treeNodf(comm=comm.aleat, col.tree, order.rows=FALSE,
                                    row.tree=row.tree.aleat, order.cols=FALSE)$cols  
                                 # although info on col.tree is not used, 
                                 # it should be present here in the list of arguments.
     } 
   } # closes else
   
} # close for permutations.


resu<-c(obs, list(rows.aleats=rows.aleats, cols.aleats=cols.aleats, 
        mat.aleats=mat.aleats, permutations=permutations))

class(resu)<-"treeNodfTest"
resu
}

