compas <- function(S, dims, am, clump=1, beta.R, coords, n.quanti=5, n.quali=0.1, add1){

if(clump < 1){ stop("Argument clump must be equal or higher than 1.")}
if(clump != round(clump)){ warning("Argument clump was not an integer and was rounded to an integer.")}
clump <- round(clump)

if(n.quanti <= 0){ 
   warning("Argument n.quanti must be higher than 0; it was substituted by 0.001 ")
   n.quanti <- 0.001
}


if(is.vector(coords)){
   n <- length(coords)
   coords <- cbind(coords, rep(0, n))}

n <- nrow(coords)
alfa <- matrix(runif(dims*S, 0.1, 5), S, dims)		
gama <- matrix(runif(dims*S, 0.1, 5), S, dims)
A <-    rlnorm(n=S, meanlog=am) 	#Abundance of each species in the modal point

#======  Position of the modal point in the gradient
m.temp <- runif((dims*S*clump), min=-50, max=150)
m <- matrix(m.temp, S*clump, dims)
if(clump > 1){
   m     <- apply(m, 2, sort)
   probs <- (1:(S*clump)) / (S*clump)
   probs <- probs^clump # exponential curve of probabilities.
   m     <- apply(m, 2, sample, size=S, prob=probs)
}
#plot(m, xlim=c(-50, 150), ylim=c(-50,150)); abline(v=c(0,100), h=c(0,100))


#==========
r <- matrix(NA, S, dims)
for(i in 1:dims){  #Range size of S spp.(rows) in dims gradients (cols)
  #r[,i]<-rnorm(S, mean=100*beta[i], sd=0.3*100*beta[i]) # version <1.7.1
   r[,i]<-rnorm(S, mean=100/beta.R[i], sd=0.3*100/beta.R[i])
} 

b <- alfa/(alfa+gama)
d.temp <- matrix(NA, S, dims)
for(g in 1:dims){d.temp[,g] <- (b[,g]^alfa[,g])*((1-b[,g])^gama[,g])}

d<- matrix(NA,S)
for(i in 1:S){d[i]<-prod(d.temp[i,])}


# Create array3D with (rows=samples, col=spp, layer=gradient) values from right-hand portion 
# (after Pi symbol) of  equation 5 in Minchin (1987).
resu.temp <- array(NA, c(n, S, dims))
for(g in 1:dims){
 for(i in 1:S){
  for(k in 1:n){
    resu.temp[k,i,g]<-(((((coords[k,g]-m[i,g])/r[i,g])+b[i,g])^alfa[i,g]) * ((1-(((coords[k,g]-m[i,g])/r[i,g])+b[i,g]))^gama[i,g]))
		}}}

resu.temp <- ifelse(is.na(resu.temp),0,resu.temp)## Substitute NAs by 0.
resu<-matrix(NA, n, S)
for(i in 1:S){
   for(k in 1:n){
     resu[k, i] <- prod(resu.temp[k, i,]) #The products of Eq.5 in Minchin (1987).
}  }        

for(k in 1:n){
   resu[k, ] <- (A/d)*resu[k,]}# Final values obtained by multiplying left- to hand-side Eq.5.
	
# Quantitative Noise
# Each x value in the dataset is replaced by a random value from a 
#   Negative Binomial distribution with mean = x and var = x*(1+n.quanti).
for(k in 1:n){
   for(i in 1:S){
      #resu[k, i] <- rpois(1, resu[k,i])}}
      value <- resu[k,i]
      if(value > 0){ # because zero produces value.size = NaN.
         vari <- value*(1+n.quanti) 
         value.size <- value^2 / (vari-value)
         resu[k, i] <- rnbinom(1, mu=value, size = value.size)
      }  
   }
}   

# Qualitative Noise
# sp has probability '1-n.quali' of occurring in a site within its range. 
# The command replaces n.quali*100% of the values by "0".
for(i in 1:S){
   resu[, i] <- rbinom(n,1,(1-n.quali))*resu[,i]
}

resu <- ifelse(resu<0, 0, resu)             #Exclude negative values.
S.pres <- sum(ifelse(colSums(resu)>0, 1, 0))#Number of species present in the simulated sample.
cols.pres <- ifelse(colSums(resu)>0, 1:S, 0)#'Names'(numbers) of columns containing observed species.
resu.f <- matrix(NA, n, S.pres)
resu.f <- resu[, cols.pres]                 #Exclude columns without observed species.

# Add (add1*100)% of 'marginal spp.' occuring randomly with 1 individual in the dataset.
single <- matrix(0, n, round(add1*ncol(resu.f)))
for(z in 1:round(add1*ncol(resu.f))){
  single[sample(1:n,1), z]<-1}

resu.f <- cbind(resu.f, single)
colnames(resu.f) <- paste("sp",   1:ncol(resu.f), sep=".")
rownames(resu.f) <- paste("site", 1:nrow(resu.f), sep=".")


# Response Curves without noises for 1 gradient simulations.
# The function does not create plot for dims>1.
if(dims==1){
  plot(x=0, y=0,xlab="Gradient",ylab="Abundance",xlim=c(0,100),ylim=c(0,max(A)),
        main="Response curves WITHOUT quanti- and qualitative noises",type="n") 
  for(i in 1:S){
    species<-function(x) A[i]/d[i] * (((x-m[i])/r[i])+b[i])^alfa[i] * (1-(((x-m[i])/r[i])+b[i]))^gama[i]
    curve(species, from=0, to=100, add=TRUE, col=(sample(1:10,1)))}}

return(resu.f)
}


