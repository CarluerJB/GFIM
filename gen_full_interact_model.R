##2 functions :

# [generate.quad.design]  generates a "quadratic" design matrix"
## with simple, interaction and pure quadra effects
## by default the simple effect matrix is binomial

#[sparse.param.quad.model] constructs a sparse parameter model by 
# breaking down the non-zero coordinates into simple, interaction and quadratic effects

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
theta_range=1
p=0.5 # ratio 0/1 (binom rule)
s<-5 # Number of signal
m<-100 # number of SNPs

for (n in c(50, 250, 500, 1000, 2500, 5000, 7500, 10000)){ # number of individuals
  # Generate X  
  X<-generate.quad.design(m,n)
  
  # Generate theta
  theta<-sparse.param.quad.model(s,m, theta_range=theta_range)
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
  write.table(x=colnames(X)[sup], file=paste0("SUPPORT_",n,".txt") , quote = FALSE, row.names = FALSE, col.names = F)
  write.table(x=X, file=paste0("X_generated_",n,".txt"), sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(x=Y.al, file=paste0("Y_generated_random_effect_",n,".txt"), row.names = FALSE, col.names = FALSE)
  write.table(x=Y.fixe, file=paste0("Y_generated_fixed_effect_",n,".txt"), row.names = FALSE, col.names = FALSE)
}
