\name{select.window}
\alias{select.window}

\title{Select cells inside a defined window of a gridded dataset}

\description{Selection of cells (or sites) around a focal cell. Cells must have have xy coordinates. Users should define a 'radius' from the focal cell that will define the window}

\usage{select.window(xf, yf, radius = 1, xydata)}
\arguments{
  \item{xf}{The x-coordinate of the focal cell.}
  \item{yf}{The y-coordinate of the focal cell.}
  \item{radius}{The radius that define the window to select neighbor cells.}
  \item{xydata}{A matrix-like object. The first and second columns of data must have x and y coordinates, respectively. The remaining columns must contain species or other atributes to be selected.}
}

\details{The function is intended to be used in gridded data, but should work on sites irregularly scattered in the xy-space.}

\value{A matrix-like object in which the first and second columns are the x and y coordinates. The remaining columns includes species (or attributes) observed in the selected cells.}

\author{Adriano Sanches Melo}


\examples{
x <- rep(1:5,each=5)
y <- rep(1:5,5)
spp <- matrix(1:100,25,4)
colnames(spp) <- c("sp1","sp2","sp3","sp4")
xyspp <- cbind(x,y,spp)
resu <- select.window(xf=3, yf=3, radius=1.1, xydata=xyspp)
resu

plot(x,y) 
# maintain the plot window open.
points(resu[,1:2],col=2,cex=2 ) 
# cells of the selected window in red.

# A reduced number of cells will be selected for focal cells located in margins.
resu <- select.window(xf=5, yf=5, radius=1.1, xydata=xyspp)
plot(x,y) 
# maintain the plot window open.
points(resu[,1:2],col=2,cex=2 ) 
# cells of the selected window in red.

# Unrecorded species in the selected window are removed from resulting 
#  dataframe (or matrix):
spp<-matrix(rep(0:1,each=50),25,4)
colnames(spp)<-c("sp1","sp2","sp3","sp4")
xyspp<-cbind(x,y,spp)
select.window(xf=3, yf=3, radius=1.1, xydata=xyspp)
}
