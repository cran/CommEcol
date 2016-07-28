print.treeNodf<- function (x, digits=3, ...){
resu<-matrix(NA, 3, 3)
rownames(resu) <- c("rows", "columns", "matrix")
colnames(resu) <- c("treeNODF", "S.fraction", "topoNODF")
resu[1,] <- x$rows
resu[2,] <- x$cols
resu[3,] <- x$mat

print(resu)
invisible(x)
}

