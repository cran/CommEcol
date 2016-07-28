treeNodf<-function (comm, col.tree, order.rows=FALSE, row.tree, order.cols=FALSE) {
comm <- ifelse(comm > 0, 1, 0)
comm <- as.matrix(comm)
rows <- rep(NA, 3)
names(rows)<-c('tree.nodf', 's.fraction', 'topo.nodf')
cols <- rows
mat  <- rows

resu <- list(
         rows = rows,
         cols = cols,
         mat  = mat,
         
         tree.nodf.rows.dist  = NA,
         tree.nodf.cols.dist  = NA,
         s.fraction.rows.dist = NA,
         s.fraction.cols.dist = NA,
         topo.nodf.rows.dist  = NA,
         topo.nodf.cols.dist  = NA
         )

# treeNodf to rows (sites).
if(missing(col.tree)==FALSE){
  if (is.null(col.tree$edge.length)) {
     stop("Tree has no branch lengths, cannot compute BL-diversity")
  }
  
  if (!is.rooted(col.tree)) {
     stop("Rooted tree required for treeNodf calculation")
  }
  
  one.rows <- treeNodfOne(comm, tree=col.tree, order=order.rows)
  resu$rows['tree.nodf']  <- mean(one.rows$tree.nodf.dist)
  resu$rows['s.fraction'] <- mean(one.rows$s.fraction.dist)
  resu$rows['topo.nodf']  <- mean(one.rows$tree.nodf.dist - one.rows$s.fraction.dist) 
  resu$tree.nodf.rows.dist  <- one.rows$tree.nodf.dist
  resu$s.fraction.rows.dist <- one.rows$s.fraction.dist
  resu$topo.nodf.rows.dist  <- one.rows$tree.nodf.dist - one.rows$s.fraction.dist
}

# treeNodf to colums (species)
if(missing(row.tree)==FALSE){
  if (is.null(row.tree$edge.length)) {
     stop("Tree has no branch lengths, cannot compute BL-diversity")
  }
  
  if (!is.rooted(row.tree)) {
     stop("Rooted tree required for treeNodf calculation")
  }
  
  one.cols<-treeNodfOne(t(comm), tree=row.tree, order=order.cols) 
  resu$cols['tree.nodf']  <- mean(one.cols$tree.nodf.dist)
  resu$cols['s.fraction'] <- mean(one.cols$s.fraction.dist)
  resu$cols['topo.nodf']  <- mean(one.cols$tree.nodf.dist - one.cols$s.fraction.dist) 
  resu$tree.nodf.cols.dist  <- one.cols$tree.nodf.dist
  resu$s.fraction.cols.dist <- one.cols$s.fraction.dist
  resu$topo.nodf.cols.dist  <- one.cols$tree.nodf.dist - one.cols$s.fraction.dist
}

# treeNODF to both margins
if(missing(col.tree)==FALSE & missing(row.tree)==FALSE) {
  
  tot.tree       <- sum(sum(one.rows$tree.nodf.dist),  sum(one.cols$tree.nodf.dist))
  tot.s.fraction <- sum(sum(one.rows$s.fraction.dist), sum(one.cols$s.fraction.dist))
  ndist <- sum(length(one.rows$tree.nodf.dist), length(one.cols$tree.nodf.dist))
  resu$mat['tree.nodf']  <- tot.tree       / ndist 
  resu$mat['s.fraction'] <- tot.s.fraction / ndist
  resu$mat['topo.nodf']  <- (tot.tree      / ndist ) - (tot.s.fraction / ndist) 
  }

class(resu) <- "treeNodf"
resu
}





