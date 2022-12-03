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
  X <- as.matrix(X)
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
  #add this to subtab:
  subtab <- cbind(subtab,aicvals,rsqrvals)

  #we need the score for the zero case to be worth 0, and we would like a higher
  #score to indicate a better model. For this reason, we adjust the aic by
  #subtracting v({}) from it, and multiplying by negative one. This is a linear
  #transformation of aic.
  subtab[,p+1] <- aicvals[1]-aicvals

  #initialize an array for the shapley values
  shapvals <- matrix(data=NA,nrow=p,ncol=2)
  colnames(shapvals)<-c("AIC","Rsq")
  rownames(shapvals)<-colnames(X)

  #loop iterates once for each variable to calculate the shapley values
  for(i in 1:p){
    aicshap <- 0
    rsqshap <- 0

    #loop iterates over the subsets to calculate the sum
    for(j in 1:nmodels){
      #look for coalitions where current variable is missing.
      if(subtab[j,i]==0){
        exc_model <- subtab[j,]
        exc_subset <- subtab[j,1:p]
        inc_subset <- exc_subset
        inc_subset[i] <- 1
        inc_model <- NULL
        #search the subest table for the correct model, grab its AIC and Rsq
        for(k in 1:nmodels){
          flag <- TRUE
          for(l in 1:p){
            if(subtab[k,l]!=inc_subset[l]){
              flag <- FALSE
            }
          }
          if(flag){inc_model<-subtab[k,]}
        }
        #print models for diagnostics
        #cat("excluding model: \n",exc_model,"\n")
        #cat("including model: \n",inc_model,"\n")
        #cat("-----------------------------------\n")

        #add the correct value to the shapley values
        aicshap <- aicshap + (factorial(sum(exc_model[1:p]))*factorial(p-sum(exc_model[1:p])-1)/factorial(p))*(inc_model[p+1]-exc_model[p+1])
        rsqshap <- rsqshap + (factorial(sum(exc_model[1:p]))*factorial(p-sum(exc_model[1:p])-1)/factorial(p))*(inc_model[p+2]-exc_model[p+2])
      }

    }
    shapvals[i,1]<-aicshap
    shapvals[i,2]<-rsqshap
  }

  #okay, we've calculated everything. Now it's time to wrap it up in a nice
  #object to output. We'll output the shapley value vectors and the subset table
  shap_out <- list()
  shap_out$AIC_shapleys <- shapvals[,1]
  shap_out$rsq_shapleys <- shapvals[,2]
  shap_out$subset_table <- as.data.frame(subtab)
  shap_out$subset_table$aicvals <- aicvals
  return(shap_out)
}
