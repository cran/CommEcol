print.treeNodfTest<- function (x, digits=4, ...){

permutations  <- x$permutations

if(is.na(x$rows[1])==FALSE){
   resu.rows <- as.data.frame(matrix(NA, 3, 5))
   rownames(resu.rows) <- c("treeNODF.rows","S.Fraction.rows","topoNODF.rows")
   colnames(resu.rows) <- c("Obs","M.aleat","SD.aleat","Z","Prob")
   
   resu.rows[, 1] <- x$rows
   resu.rows[, 2] <- colMeans(x$rows.aleats)
   resu.rows[, 3] <- apply(x$rows.aleats, 2, sd)
   resu.rows[, 4] <- (resu.rows[, 1] - resu.rows[, 2]) /  resu.rows[, 3]

   ab<-function(x, a) {x >= a}
   cases.rows      <- colSums(t(apply(x$rows.aleats, 1, ab, a=x$rows)))
   resu.rows[, 5]  <- (cases.rows+1) / (permutations+1)
}


if(is.na(x$cols[1])==FALSE){
   resu.cols <- as.data.frame(matrix(NA, 3, 5))
   rownames(resu.cols) <- c("treeNODF.cols","S.Fraction.cols","topoNODF.cols")
   colnames(resu.cols) <- c("Obs","M.aleat","SD.aleat","Z","Prob")
   
   resu.cols[, 1] <- x$cols
   resu.cols[, 2] <- colMeans(x$cols.aleats)
   resu.cols[, 3] <- apply(x$cols.aleats, 2, sd)
   resu.cols[, 4] <- (resu.cols[, 1] - resu.cols[, 2]) /  resu.cols[, 3]

   ab<-function(x, a) {x >= a}
   cases.cols      <- colSums(t(apply(x$cols.aleats, 1, ab, a=x$cols)))
   resu.cols[, 5]  <- (cases.cols+1) / (permutations+1)
}


if(is.na(x$rows[1])==FALSE & is.na(x$cols[1])==FALSE) {
   resu.mat <- as.data.frame(matrix(NA, 3, 5))
   rownames(resu.mat) <- c("treeNODF.mat","S.Fraction.mat","topoNODF.mat")
   colnames(resu.mat) <- c("Obs","M.aleat","SD.aleat","Z","Prob")
   
   resu.mat[, 1] <- x$mat
   resu.mat[, 2] <- colMeans(x$mat.aleats)
   resu.mat[, 3] <- apply(x$mat.aleats, 2, sd)
   resu.mat[, 4] <- (resu.mat[, 1] - resu.mat[, 2]) /  resu.mat[, 3]

   ab<-function(x, a) {x >= a}
   cases.mat      <- colSums(t(apply(x$mat.aleats, 1, ab, a=x$mat)))
   resu.mat[, 5]  <- (cases.mat+1) / (permutations+1)
}


if(is.na(x$rows[1])==FALSE & is.na(x$cols[1])==TRUE) {print(resu.rows)}
if(is.na(x$rows[1])==TRUE  & is.na(x$cols[1])==FALSE){print(resu.cols)}
if(is.na(x$rows[1])==FALSE & is.na(x$cols[1])==FALSE){print(list(resu.rows=resu.rows,
                                                resu.cols=resu.cols, resu.mat=resu.mat))}


invisible(x)
}

