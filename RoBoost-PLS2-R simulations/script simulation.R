
library(R.matlab)
library(MASS)
library(plotly)

library(ggplot2)

#### basic functions
bruit.struc <- function(ncol = NULL,ngauss = NULL,nstruc = NULL,var=NULL){
  
  if(is.null(var)){var = ncol}
  
  seqlast <- function (from, to, by) 
  {
    vec <- do.call(what = seq, args = list(from, to, by))
    if ( tail(vec, 1) != to ) {
      return(c(vec, to))
    } else {
      return(vec)
    }
  }
  
  
  
  mat.bruit <- matrix(ncol = ncol, nrow = nstruc) 
  
  for (a in seq(1:nstruc)){
    
    div = round((c(length(var)/ngauss)),0)
    x <- c(1:ncol)[var[sample(1:length(var),ngauss,replace= TRUE)]]
    
    pos <- sample(c(1),size = length(x),replace =TRUE)
    sd<- sample(10:200,size = length(x),replace =TRUE)
    vec <- rep(0, ncol)
    z = 1:ncol
    for (i in seq(1:length(x))){
      
      gauss <- dnorm(x = z,mean = x[i],sd = sd[i])
      vec <- (vec + gauss )*pos[i]
    }
    
    mat.bruit[a,] <- vec
  }
  return(mat.bruit)
}


sim.bruit <- function(n,Xr=Xr,sd.b = sb.d) {
  
  l <- ncol(Xr)
  d <- lapply(1:l,FUN = function(x,n,Xr){
    mean <- 0
    sd <- (max(Xr[1,])-min(Xr[1,]))/100*sd.b
    rnorm(n = n,mean = mean,sd = sd)
  },Xr=Xr,n=n)
  
  d <- as.data.frame(d)
  colnames(d) <- NULL
  return(d)
}

sim.bruitY = function(n,X,sd.b) {
  mean = 0
  d = matrix(nrow = nrow(X),ncol = ncol(X))
  for (i in (1:ncol(d))){
    sd = sd(X[,i])*sd.b
    d[,i] = rnorm(n = n,mean = mean,sd = sd)
  }
  
  d <- as.data.frame(d)
  colnames(d) <- NULL
  return(d)
}




#### Simulation 1 ####

## pur spectra extraction

data.al <- readMat("data_alcool.mat")
Y <- matrix(ncol = 3 , nrow = 22)
Y[,1] <- data.al$Y.EtOH
Y[,2] <- data.al$Y.Glu
Y[,3] <- data.al$Y.Eau
X <- data.al$x
colnames(X) <- data.al$lbd
Y <- Y[,1:3]/c(1000,1000,1000)
Y <- as.data.frame(Y)
Y$inter1.3 <- (Y[,1] * Y[,3])
K = t(X)%*%t(ginv(as.matrix(Y)))
K <- t(K)
colnames(K) <- data.al$lbd


# T(scores) simulation for the pur spectra

Ts <- as.data.frame((as.matrix(rnorm(1000,20,5)/10)))
Ts$"2" <-  (as.matrix(rnorm(1000,20,5)/20))
Ts$"3" <- (as.matrix(rnorm(1000,20,5)))
T <- Ts


# X simulation with the pur spectra

X = as.matrix(T[,1:3])%*%K[c(1,2,3),]
X = X/100

# artifical spectra simulation

nstruc = 10
struc <- bruit.struc(ncol = 1001,ngauss = 7,nstruc = nstruc,var = c(400:1001))
colnames(struc) <- data.al$lbd

# addition of the simulated spectra to the simulation of X
X = X + matrix(abs(rnorm(1000*nstruc,0,1)),nrow = 1000,ncol = nstruc)%*%struc
X = X  + as.matrix(sim.bruit(1000,Xr = X,sd.b = 0.1))


# noise simulation of Y

Y = T[,1:3]*10
b = sim.bruitY(1000,X = Y,sd.b = 0.15)

# outliers creation

Y[1:200,3] = T[1:200,3]*-5
Y = Y + b


save(X,file = "X.sim1.rda")
save(Y,file = "Y.sim1.rda")
save(b,file =  "b.sim1.rda")


##### Simulation 2 ####

## pur spectra extraction

data.al <- readMat("data_alcool.mat")
Y <- matrix(ncol = 3 , nrow = 22)
Y[,1] <- data.al$Y.EtOH
Y[,2] <- data.al$Y.Glu
Y[,3] <- data.al$Y.Eau
X <- data.al$x
colnames(X) <- data.al$lbd
Y <- Y[,1:3]/c(1000,1000,1000)
Y <- as.data.frame(Y)
Y$inter1.3 <- (Y[,1] * Y[,3])
K = t(X)%*%t(ginv(as.matrix(Y)))
K <- t(K)
colnames(K) <- data.al$lbd



# T(scores) simulation for the pur spectra

Ts <- as.data.frame((as.matrix(rnorm(1000,20,5)/10)))
Ts$"2" <-  (as.matrix(rnorm(1000,20,5)/20))
Ts$"3" <- (as.matrix(rnorm(1000,20,5)))
T <- Ts

# X simulation with the pur spectra

X = as.matrix(T[,1:3])%*%K[c(1,2,3),]
X = X/100


# artifical spectra simulation

nstruc = 10
X2 <- X
struc <- bruit.struc(ncol = 1001,ngauss = 7,nstruc = nstruc,var = c(400:1001))
colnames(struc) <- data.al$lbd

# addition of the simulated spectra to the simulation of X

X = X + matrix(abs(rnorm(1000*nstruc,0,1)),nrow = 1000,ncol = nstruc)%*%struc

X = X  + as.matrix(sim.bruit(1000,Xr = X,sd.b = 0.1))


# Y simulation and outliers generation

Y = T[,1:3]*10

b = sim.bruitY(1000,X = Y,sd.b = 0.15)
Y[1:200,2] = T[1:200,2]*-10
Y = Y + b

save(X,file = "X.sim2.rda")
save(Y,file = "Y.sim2.rda")
save(b,file =  "b.sim2.rda")

##### Simulation 3 ####

## pur spectra extraction

data.al <- readMat("data_alcool.mat")
Y <- matrix(ncol = 3 , nrow = 22)
Y[,1] <- data.al$Y.EtOH
Y[,2] <- data.al$Y.Glu
Y[,3] <- data.al$Y.Eau
X <- data.al$x
colnames(X) <- data.al$lbd
Y <- Y[,1:3]/c(1000,1000,1000)
Y <- as.data.frame(Y)
Y$inter1.3 <- (Y[,1] * Y[,3])
K = t(X)%*%t(ginv(as.matrix(Y)))
K <- t(K)
colnames(K) <- data.al$lbd



# T(scores) simulation for the pur spectra


Ts <- as.data.frame((as.matrix(rnorm(1000,20,5)/10)))
hist(Ts$V1)



Ts$"2" <-  (as.matrix(rnorm(1000,20,5)/20))
Ts$"3" <- (as.matrix(rnorm(1000,20,5)))
T <- Ts

# X simulation with the pur spectra

X = as.matrix(T[,1:3])%*%K[c(1,2,3),]
X = X/100

# artificial spectral signatures simulation and X simulation

nstruc = 10
X2 <- X
struc <- bruit.struc(ncol = 1001,ngauss = 7,nstruc = nstruc,var = c(400:1001))
colnames(struc) <- data.al$lbd


X = X + matrix(abs(rnorm(1000*nstruc,0,1)),nrow = 1000,ncol = nstruc)%*%struc

X = X  + as.matrix(sim.bruit(1000,Xr = X,sd.b = 0.1))


# outliers simulation

nstruc = 10
struc <- bruit.struc(ncol = 1001,ngauss = 3,nstruc = nstruc,var = c(400:1001))
colnames(struc) <- data.al$lbd
X[1:200,] = X[1:200,] + matrix(abs(rnorm(1000*nstruc,0,3)),nrow = 200,ncol = nstruc)%*%struc



# Y simulation 

Y = T[,1:3]*10
b = sim.bruitY(1000,X = Y,sd.b = 0.15)
Y = Y + b

save(X,file = "X.sim3.rda")
save(Y,file = "Y.sim3.rda")
save(b,file =  "b.sim3.rda")
