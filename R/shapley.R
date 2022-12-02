#This function takes in a response vector, a design matrix, and calculates the
#shapley value for each column of the design matrix using either AIC or R^2.

#Test data
X <- cbind(rnorm(100),runif(100,1,5),1-rexp(100),rnorm(100))
colnames(X)<-c("X1","X2","X3","X4")
Y <- 1 + X%*%cbind(c(1,1/2,-1,0)) + rnorm(100)

shapley <- function(Y,X){
  #First, we need to calculate AIC for all possible models. We set up a logical
  #matrix which indicates inclusion (TRUE or FALSE) of each variable for all 2^p
  #possible coalitions. Since this is equivalent to counting from 0 to (2^p)-1
  #in binary, we leverage this to create this table, called subtab
  p <- ncol(X)
  nmodels <- 2^(p)
  subtab <- matrix(data=NA,ncol=p,nrow=2^p)
  colnames(subtab)<-colnames(X)
  for(i in 1:p){
    subtab[,i]<-rep(
      c(rep(FALSE,2^(i-1)),
        rep(TRUE, 2^(i-1))),
      2^(p-i))
  }

  #Next, we need to iterate over this table, fit a linear model for every
  #possible model, and calculate both AIC and Rsquared for each model
  aicvals  <- numeric(nmodels)
  rsqrvals <- numeric(nmodels)
  for(i in 1:nmodels){
    #subtab is constructed so the first row will always be the empty model. So
    #we use an intercept only model in this case:
    if(i==1){
      lmod <- lm(Y~1)
    }
    else{ #otherwise use subtab for design matrix subsetting
      lmod <- lm(Y ~ X[,subtab[i,]] )
    }
    aicvals[i]  <- AIC(lmod)
    rsqrvals[i] <- summary(lmod)$r.squared
  }
}
