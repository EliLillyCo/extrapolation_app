#!/lrlhps/apps/R/qualified/R-3.4.0/bin/Rscript 
### post.R --- 
## 
## Filename: post.R
## Description: 
## Author: Michael D Sonksen
## Maintainer: 
## License: This is owned by Eli Lilly and Company
## Created: Mon Jan 21 20:09:17 2019 (-0500)
## Version: 
## Package-Requires: ()
## Last-Updated: 
##           By: 
##     Update #: 0
## URL: 
## Doc URL: 
## Keywords: 
## Compatibility: 
## 
######################################################################
## 
### Commentary: 
## 
## 
## 
######################################################################
## 
### Change Log:
## 
## 
######################################################################
## 
## 
## 
######################################################################
## 
### Libraries:


## 
### Global Functions: 


## 
### Code: 


#' Simulate sufficient statistics from two independent normal distributions.  
#'
#' @param n_sims number of sets of samples to return.
#' @param mu vector of two means for the simulations.
#' @param sigma vector of two standard deviations.
#' @param n vector of two samples sizes. 
#'
#' @export
simulate_studies <- function(n_sims,mu,sigma,n){
    sims_ybar <- t(matrix(rnorm(n_sims*2,mu,sigma/sqrt(n)),
                          nrow=2))
    sims_s2 <- t(matrix(rchisq(n_sims*2,n-1),nrow=2)*sigma*sigma/(n-1))
    list(sims_ybar=sims_ybar, sims_s2=sims_s2)
}


#' Fit a Bayesian model with a flat prior. Returns the average P(d[2]>eoi) and
#' the percent of the time P(d[2]>eoi)>prth.
#'
#' @param sims a list created with simulate_studies
#' @param n sample sizes for the samples in simulate_studies
#' @param eoi effect of interest.
#' @param prth probability threshold.
#' @param lower.tail direction for the probability (default=FALSE). 
#' 
#' @export
fit_flat <- function(sims,n,eoi,prth,lower.tail=FALSE){
    mn <- sims$sims_ybar[,2]-sims$sims_ybar[,1]
    se <- sqrt((2.0/n[2])*sims$sims_s2[,2])
    c(mean(1-pnorm(eoi,mn,se,lower.tail=lower.tail)),
      mean((1-pnorm(eoi,mn,se,lower.tail=lower.tail))> prth))
}

#' Fit a Bayesian model with an informative prior. Returns the average
#' P(d[2]>eoi) and the percent of the time P(d[2]>eoi)>prth.
#'
#' @param sims a list created with simulate_studies
#' @param n sample sizes for the samples in simulate_studies
#' @param eoi effect of interest.
#' @param prth probability threshold.
#' @param mu0 the prior mean of d[2].
#' @param sigma0 the prior sd of d[2]. 
#' @param lower.tail direction for the probability (default=FALSE). 
#' 
#' @export
fit_inf <- function(sims,n,eoi,prth,mu0,sigma0,lower.tail=FALSE){
    mn <- sims$sims_ybar[,2]-sims$sims_ybar[,1]
    se <- sqrt(sims$sims_s2[,2])
    post_var <- 1.0/(1.0/sigma0^2+n[2]/(2.0*se*se))
    post_mn <- post_var*(
        (1.0/sigma0^2)*mu0 + (n[2]/(2.0*se*se))*mn
    )
    post_se <- sqrt(post_var)
    c(mean(1-pnorm(eoi,post_mn,post_se,lower.tail=lower.tail)),
      mean((1-pnorm(eoi,post_mn,post_se,lower.tail=lower.tail))> prth))
}

#' Fit a Bayesian model with a power prior. Returns the average
#' P(d[2]>eoi) and the percent of the time P(d[2]>eoi)>prth.
#'
#' @param a0 the power prior parameter. 
#' @param sims a list created with simulate_studies
#' @param n sample sizes for the samples in simulate_studies
#' @param eoi effect of interest.
#' @param prth probability threshold.
#' @param mu0 the prior mean of d[2].
#' @param sigma0 the prior sd of d[2]. 
#' @param lower.tail direction for the probability (default=FALSE). 
#' 
#' @export
fit_power <- function(a0,sims,n,eoi,prth,mu0,sigma0,lower.tail=FALSE){
    if(a0==0){
        fit_inf(sims,n,eoi,prth,mu0,sigma0,lower.tail=FALSE)
    }
    else{
        sigma0 <- sigma0/sqrt(a0)
        mn <- sims$sims_ybar[,2]-sims$sims_ybar[,1]
        se <- sqrt(sims$sims_s2[,2])
        post_var <- 1.0/(1.0/sigma0^2+n[2]/(2.0*se*se))
        post_mn <- post_var*(
            (1.0/sigma0^2)*mu0 + (n[2]/(2.0*se*se))*mn
        )
        post_se <- sqrt(post_var)
        c(1-mean(pnorm(eoi,post_mn,post_se,lower.tail=lower.tail)),
          mean((1-pnorm(eoi,post_mn,post_se,lower.tail=lower.tail))> prth))
    }
}

marginal <- function(ybar,s2,mu0,sigma0,n,sigma,log=FALSE){
    Sigma <- diag(c(rep(sigma[1]^2,n[1]),rep(sigma[2]^2,n[2])))
    Sigmai <- diag(c(rep(1.0/sigma[1]^2,n[1]),rep(1.0/sigma[2]^2,n[2])))
    Sigma0 <- diag(sigma0*sigma0)
    Sigma0i <- diag(1/(sigma0*sigma0))
    X <- cbind(1,c(rep(0,n[1]),rep(1,n[2])))
    
    betahat <- c(ybar[1],ybar[2]-ybar[1])
    Xbetahat <- c(rep(ybar[1],n[1]),rep(ybar[2],n[2]))
    ##Xbetahat2 <- X%*%betahat
    ##mean(Xbetahat==Xbetahat2)
    ##XtSigmaiX <- t(X)%*%Sigmai%*%X
    XtSigmaiX <- matrix(c(
        n[2]/sigma[2]^2+n[1]/sigma[1]^2,
        rep(n[2]/sigma[2]^2,3)),2,2)
        
    log_ll <- -(sum(n)/2.0)*log(2*pi) +
        0.5*log(det(XtSigmaiX+Sigma0i))-
        0.5*(n[1]*2*log(sigma[1])+n[2]*2*log(sigma[2]))-
        0.5*(log(sigma0[1])+log(sigma0[2]))
    
    log_ll <- log_ll-0.5*(n[1]-1)*s2[1]/sigma[1]^2+(n[2]-1)*s2[2]/sigma[2]^2
    log_ll <- log_ll-0.5*t(mu0)%*%Sigma0i%*%mu0
    log_ll <- log_ll-0.5*t(Xbetahat)%*%Sigmai%*%Xbetahat
    log_ll <- log_ll+0.5*t(Sigma0i%*%mu0+t(X)%*%Sigmai%*%Xbetahat
                           )%*%solve(XtSigmaiX+Sigma0i
                                     )%*%(Sigma0i%*%mu0+t(X)%*%Sigmai%*%Xbetahat)
    if(log)
        return(log_ll)
    else
        return(exp(log_ll))
}



marginal2 <- function(ybar,s2,mu0,sigma0,n,sigma){
    Sigma <- diag(c(rep(sigma[1]^2,n[1]),rep(sigma[2]^2,n[2])))
    Sigma0 <- diag(sigma0*sigma0)
    X <- cbind(1,c(rep(0,n[1]),rep(1,n[2])))
    
    betahat <- c(ybar[1],ybar[2]-ybar[1])
    log_ll <- -(sum(n)/2.0)*log(2*pi) +
        0.5*log(det(t(X)%*%solve(Sigma)%*%X+solve(Sigma0)))-
        0.5*log(det(Sigma))-0.5*log(det(Sigma0))

    log_ll <- log_ll-0.5*(n[1]-1)*s2[1]/sigma[1]+(n[2]-1)*s2[2]/sigma[2]
    log_ll <- log_ll-0.5*t(mu0)%*%solve(Sigma0)%*%mu0
    log_ll <- log_ll-0.5*t(X%*%betahat)%*%solve(Sigma)%*%X%*%betahat
    log_ll <- log_ll+0.5*t(solve(Sigma0)%*%mu0+t(X)%*%solve(Sigma
        )%*%X%*%betahat)%*%solve(t(X)%*%solve(Sigma)%*%X+solve(Sigma0)
        )%*%(solve(Sigma0)%*%mu0+t(X)%*%solve(Sigma)%*%X%*%betahat)
    exp(log_ll)
}

#' Fit a Bayesian model with a mixture prior. Returns the average
#' P(d[2]>eoi) and the percent of the time P(d[2]>eoi)>prth.
#'
#' @param p0 the prior mixing prob. 
#' @param sims a list created with simulate_studies
#' @param n sample sizes for the samples in simulate_studies
#' @param eoi effect of interest.
#' @param prth probability threshold.
#' @param mu0 vector of  prior mean of d[1] and d[2] for the first mixture piece.
#' @param sigma0 vector of prior sd of d[1] and d[2] for the first mixture piece. 
#' @param mu1 vector of  prior mean of d[1] and d[2] for the 2nd mixture piece.
#' @param sigma1 vector of prior sd of d[1] and d[2] for the 2nd mixture piece. 
#' @param lower.tail direction for the probability (default=FALSE). 
#' 
#' @export
fit_mix <- function(p0,sims,n,eoi,prth,
                    mu0,sigma0,mu1,sigma1,
                    lower.tail=FALSE){
    if(length(p0)==1)
        p0 <- c(p0,1-p0)
    
    ybar <- cbind(sims$sims_ybar[,1],sims$sims_ybar[,2])
    s2 <- cbind(sims$sims_s2[,1],sims$sims_s2[,2])

    
    m0 <- sapply(1:dim(ybar)[1],
                 function(x) marginal(ybar[x,],s2[x,],mu0,sigma0,n,sigma=sqrt(s2[x,]),log=TRUE))
    m1 <- sapply(1:dim(ybar)[1],
                 function(x) marginal(ybar[x,],s2[x,],mu1,sigma1,n,sigma=sqrt(s2[x,]),log=TRUE))
    ##need to add largest for each.
    max_m <- apply(cbind(m0,m1),1,max)
    m0 <- exp(m0-max_m)
    m1 <- exp(m1-max_m)
        
    ybar <- cbind(sims$sims_ybar[,1],sims$sims_ybar[,2]-sims$sims_ybar[,1]) ##marginal does the conversion
    p <- p0[1]*m0/(p0[1]*m0+p0[2]*m1)
    p <- cbind(p,1-p)

    
    post_var <- cbind(1.0/(1.0/sigma0[2]^2+n[2]/(2.0*s2[,2])),
                      1.0/(1.0/sigma1[2]^2+n[2]/(2.0*s2[,2])))

    post_mn <- cbind(post_var[,1]*((1.0/sigma0[2]^2)*mu0[2] +
                                   (n[2]/(2.0*s2[,2]))*ybar[,2]),
                     post_var[,2]*((1.0/sigma1[2]^2)*mu1[2] +
                                   (n[2]/(2.0*s2[,2]))*ybar[,2]))

    prb1 <- 1-pnorm(eoi,post_mn[,1],sqrt(post_var[,1]),
                    lower.tail=lower.tail)
    prb2 <- 1-pnorm(eoi,post_mn[,2],sqrt(post_var[,2]),
                          lower.tail=lower.tail)
    prb <- p[,1]*(prb1)+
           p[,2]*(prb2)
    ##c(mean(prb1>prth),mean(prb2>prth),mean(prb>prth))
    c(mean(prb),mean(prb>prth))
}



#' Simple function for running the model fitting functions alone. 
#'
#' 
main <- function(){
    n_sims <- 1000
    mu <- c(0.0,-1.0)
    sigma <- 6.0
    n <- c(150,150)
    p0 <- .5
    sims <- simulate_studies(n_sims,mu,sigma,n)

    eoi <- 0
    prth <- 0.95
    fit_inf(sims,n,eoi,prth,-1.97,1.0)
    fit_inf(sims,n,eoi,prth,0.0,1000.0)
    fit_power(0.001,sims,n,eoi,prth,0.0,1.0)
    fit_power(1.0,sims,n,eoi,prth,0.0,1.0)
    mu0 <- c(0,0)
    mu1 <- c(0,-1.97)
    sigma0 <- c(100,1)
    sigma1 <- c(100,.33)

    fit_flat(sims,n,eoi,prth)
    fit_inf(sims,n,eoi,prth,mu0[2],sigma0[2])
    fit_inf(sims,n,eoi,prth,mu1[2],sigma1[2])

    fit_mix(0.25,sims,n,eoi,prth,mu0,sigma0,mu1,sigma1)
    fit_mix(0.5,sims,n,eoi,prth,mu0,sigma0,mu1,sigma1)
    system.time(fit_mix(0.75,sims,n,eoi,prth,mu0,sigma0,mu1,sigma1))
##   user  system elapsed 
##  2.099   0.284   2.383


    ##jondavid example
    

    sims <- simulate_studies(n_sims=1000,mu=c(0,-0.7),sigma=c(3.5,3.5),n=c(150,150))
    fit_mix(0.5,
            sims,n=c(150,150),
            eoi=0,prth=0.975,
            mu0=c(0,0),
            sigma0=c(10,3),
            mu1=c(0,-1.97),
            sigma1=c(10,0.275))
    
    fit_flat(sims,n=c(150,150),eoi=0,prth=0.975)
    fit_inf(sims,n,eoi=0.0,prth=0.975,0.0,1000.0)

    sims1 <- list(sims_ybar=cbind(c(0,0),c(-.7,-.7)),
                  sims_s2=cbind(c(3.5,3.5)^2,c(3.5,3.5)^2))

    fit_mix(0.5,
            sims1,n=c(150,150),
            eoi=0,prth=0.975,
            mu0=c(0,0),
            sigma0=c(100,1),
            mu1=c(0,-1.97),
            sigma1=c(100,0.275))
    
    fit_inf(sims,n,eoi,prth,0.0,10.0)
    fit_inf(sims,n,eoi,prth,-4.0,1.0)
    
    
    
    ybar <- c(1,10)
    s2 <- c(4,4)
    mu0 <- c(0,0)
    sigma0 <- c(5,5)
    mu1 <- c(0,3)
    sigma1 <- c(5,5)
    n <- c(100,120)
    sigma <- c(4,4)
    system.time(marginal(ybar,s2,mu0,sigma0,n,sigma))
    system.time(marginal2(ybar,s2,mu0,sigma0,n,sigma))
    ybar <- c(0,-4)
    s2 <- c(3,3)
    marginal(ybar,s2,mu0,sigma0,n,sigma,log=TRUE)
    marginal(ybar,s2,mu0,sigma0,n,sigma,log=TRUE)
    
    n <- c(100,120)
    X <- cbind(1,c(rep(0,n[1]),rep(1,n[2])))
    Y <- rnorm(sum(n),c(1,5)[X[,2]+1])
    betahat <- solve(t(X)%*%X)%*%t(X)%*%Y
    Sigma <- matrix(c(1,0,0,3),2,2)
    Sigma <- diag(c(rep(1,n[1]),rep(3,n[2])))
    t(Y-X%*%betahat)%*%solve(Sigma)%*%(Y-X%*%betahat)
    var(Y[X[,2]==0])+var(Y[X[,2]==1])
    
    mean(Y[X[,2]==0])
    mean(Y[X[,2]==1])-mean(Y[X[,2]==0])
}


######################################################################
### post.R ends here
