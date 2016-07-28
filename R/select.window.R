select.window<-function(xf, yf, radius=1, xydata){

x<-xydata[,1]
y<-xydata[,2]

sele<-which(x>=(xf-radius) & x<=(xf+radius) &  y>=(yf-radius) & y<=(yf+radius) )
n<-length(sele)

if(n==1){selected<-matrix(xydata[sele,], nrow=1)
         columns<-c(1,2,which(selected>0))
         columns<-unique(columns)
         resu<-matrix(selected[1,columns], nrow=1)}
        

if(n>1){
   f<-which(x==xf & y==yf) ##focal cell
   sele.s<-sele[sele!=f]   ##remove focal cell
   resu0<-rbind(xydata[f,], xydata[sele.s,] ) #focal cell in the first row.
   columns<-c(1,2,which(colSums(resu0)>0))
   columns<-unique(columns)
   resu<-resu0[,columns]
   }

return(resu)
}

