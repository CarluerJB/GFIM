##2 functions :

# [generate.quad.design]  generates a "quadratic" design matrix"
## with simple, interaction and pure quadra effects
## by default the simple effect matrix is binomial

#[sparse.param.quad.model.many] construit un param?tre sparse en
#ventilant les coord non nulles entre effets simple, d'interaction, quadra et mixte

#[sparse.param.quad.model] constructs a sparse parameter model by 
# breaking down the non-zero coordinates into simple, interaction and quadratic effects

sparse.param.quad.model.many<-function(s,m, theta_range=theta_range){
  K<-m
  pu.fi.or.rate = 0.4
  pu.int.rate = 0.20
  mix.fi.or.int.rate = 0.20
  pu.int.mul.rate = 0.20
  
  

  pu.fi.or<-floor(s*pu.fi.or.rate) #pure first order
  pu.int = floor(s*pu.int.rate) #pure interaction
  mix.fi.or.int = floor(s*mix.fi.or.int.rate)  #mix between first order and interaction
  
  pu.int.mul<- s-(pu.fi.or+pu.int+mix.fi.or.int)
  
  batch_available_2D = c((K+1):(K*(K+1)/2))
  batch_available_1D = c(1:K)
  
  sup4<-c(4039, 4056) # FIXED FOR BETTER VISUALISATION
  sup3 = c(1605, 58) # FIXED FOR BETTER VISUALISATION
  
  
  
  sup1<-sample(batch_available_1D,pu.fi.or)
  sup2<-sample(batch_available_2D,pu.int)
  for(elem in sup2){
    batch_available_2D = batch_available_2D[batch_available_2D!=elem]
  }
  theta <- matrix(0,K+K*(K-1)/2, 1)
  sup<-sort(c(sup1,sup2,sup3,sup4))
  theta[sup]<-sample(c(-1,1),s,replace=TRUE)

  return(theta)
}

sparse.param.quad.model<-function(s,m, interOnly=FALSE, theta_range=theta_range){
  K<-m
  if(interOnly){
    # only diag pure interaction
    pu.int<-s #pure interaction
    pu.fi.or<-0 #pure first order
  }else{
    #mixed model
    pu.int<-floor(s/2) #pure interaction
    pu.fi.or<-s-floor(s/2) #pure first order
  }
  
  sup1<-sample(1:K,pu.fi.or)
  sup2<-sample(((K+1):(K*(K+1)/2)),pu.int)
  theta <- matrix(0,K+K*(K-1)/2, 1)
  sup<-sort(c(sup1,sup2))
  theta[sup]<-sample(c(-1,1),s,replace=TRUE)
  return(theta)
}

generate.quad.design<-function(m,n){
  
  K<-m # K is nb of SNP
  idx=c(1:K)
  
  X<-matrix(rbinom(K*n,1,p),nrow=n) # X binomial matrix with p the success rate (centered using p)
  
  X.int<-matrix(data = NA, nrow = dim(X)[1], byrow = FALSE) # Create empty interaction matrix
  idx.int<-matrix(data = NA, nrow = 1, byrow = FALSE)
  X = scale(X, scale = F, center = T)
  print(X)
  for (i in 1:(K-1)){ # for each SNPs id
    temp<-matrix(X[,i],ncol=(K-i),nrow=dim(X)[1],byrow=F) # copy (K+1)-i time the SNP column with diag
    print(temp)
    X.int<-cbind(X.int,temp) # save it to tab
    temp_idx<-matrix(idx[i], ncol=K-i, nrow=1, byrow=F)
    idx.int<-c(idx.int, temp_idx)
  }
  
  idx.int = idx.int[-1]
  X.int<-X.int[,-1] # remove first empty column (needed for initialization)
  
  T<-matrix(data = NA, nrow = dim(X)[1], byrow = FALSE)
  idx.T<-matrix(data = NA, nrow = 1, byrow = FALSE)
  
  for (i in 1:(K-1)){ # for each SNPs id except the last one
    temp<-X[,-(1:i)] # assign to temp the k-1 snp's column
    T<-cbind(T,temp) # save it
    temp_idx<-idx[-(1:i)]
    idx.T<-c(idx.T, temp_idx)
  }
  T<-T[,-1]
  idx.T = idx.T[-1]
  
  X.quad<-cbind(X,(T*X.int)*2)
  colnames(X.quad)<-c(paste0(idx,'-', idx),paste0(idx.int, '-', idx.T))
  
  return(X.quad)
}


h2 <- 0.4
s2 <- 0.95 * h2
theta_range = 1
p = 0.5 # ratio 0/1 (binom rule)
s <- 10 # Number of signal
m <- 100 # number of SNPs
n <- 1000 # Number of individuals


# Generate X  
X<-generate.quad.design(m,n)

# Generate theta
theta<-sparse.param.quad.model.many(s,m, theta_range=theta_range)
var.theta <- var(X%*%theta)
sup<-which(theta!=0)


# Simple effect model
# p fixed
noise.fixe=(1-s2)*var.theta/s2

Y.fixe<-X%*%theta+matrix(rnorm(nrow(X), 0, sqrt(noise.fixe)), nrow(X), 1)


# Random effect model
# p follows a normal distribution with small variance
noise.al= noise.fixe
trK = sum(X^2)
var.poly.effect = (var.theta * ( -1 + (h2/(1-h2))*((1-s2)/s2))) / trK
g <- sqrt(var.poly.effect[[1]]) * X %*% rnorm(dim(theta)[1], mean=0)

Y.al <- (X%*%theta) + g + matrix(rnorm(nrow(X), 0, sqrt(noise.al)), nrow(X), 1)

# Show some result var
var(Y.al)
var(Y.fixe)

# Save results
write.table(x=colnames(X)[sup], file="SUPPORT.txt" , quote = FALSE, row.names = FALSE, col.names = F)
write.table(x=X, file="X_generated.txt", sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(x=Y.al, file="Y_generated_random_effect.txt", row.names = FALSE, col.names = FALSE)
write.table(x=Y.fixe, file="Y_generated_fixed_effect.txt", row.names = FALSE, col.names = FALSE)
