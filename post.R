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


#' Simulate sufficient statistics. 
#'
#'
#' @export
simulate_studies <- function(n_sims,mu,sigma,n){
    sims_ybar <- t(matrix(rnorm(n_sims*2,mu,sigma/sqrt(n)),
                          nrow=2))
    sims_s2 <- t(matrix(rchisq(n_sims*2,n-1),nrow=2)*sigma/(n-1))
    list(sims_ybar=sims_ybar, sims_s2=sims_s2)
}

fit_flat <- function(sims,n,eoi,prth,lower.tail=FALSE){
    mn <- sims$sims_ybar[,2]-sims$sims_ybar[,1]
    se <- sqrt((2.0/n[2])*sims$sims_s2[,2])
    c(mean(pnorm(eoi,mn,se,lower.tail=lower.tail)),
      mean(pnorm(eoi,mn,se,lower.tail=lower.tail)> prth))
}

fit_inf <- function(sims,n,eoi,prth,mu0,sigma0,lower.tail=FALSE){
    mn <- sims$sims_ybar[,2]-sims$sims_ybar[,1]
    se <- sqrt(sims$sims_s2[,2])
    post_var <- 1.0/(1.0/sigma0^2+n[2]/(2.0*se*se))
    post_mn <- post_var*(
        (1.0/sigma0^2)*mu0 + (n[2]/(2.0*se*se))*mn
    )
    post_se <- sqrt(post_var)
    c(mean(pnorm(eoi,post_mn,post_se,lower.tail=lower.tail)),
      mean(pnorm(eoi,post_mn,post_se,lower.tail=lower.tail)> prth))
}


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
        c(mean(pnorm(eoi,post_mn,post_se,lower.tail=lower.tail)),
          mean(pnorm(eoi,post_mn,post_se,lower.tail=lower.tail)> prth))
    }
}

marginal <- function(p0,sims,mu0,sigma0){
    mn <- sims$sims_ybar[,2]-sims$sims_ybar[,1]
    se <- sqrt(sims$sims_s2[,2])
    post_var <- 1.0/(1.0/sigma0^2+n[2]/(2.0*se*se))
    post_mn <- post_var*(
        (1.0/sigma0^2)*mu0 + (n[2]/(2.0*se*se))*mn
    )
    post_se <- sqrt(post_var)

    
    c(mean(pnorm(eoi,post_mn,post_se,lower.tail=lower.tail)),
      mean(pnorm(eoi,post_mn,post_se,lower.tail=lower.tail)> prth))
}

fit_mix <- function(p0,sims,n,eoi,prth,mu0,sigma0,lower.tail=FALSE){
    mn <- sims$sims_ybar[,2]-sims$sims_ybar[,1]
    se <- sqrt(sims$sims_s2[,2])
    
    post_var1 <- 1.0/(1.0/sigma0[1]^2+n[2]/(2.0*se*se))
    post_mn1 <- post_var1*(
        (1.0/sigma0[1]^2)*mu0[1] + (n[2]/(2.0*se*se))*mn
    )
    post_se1 <- sqrt(post_var1)

    
    post_var2 <- 1.0/(1.0/sigma0[2]^2+n[2]/(2.0*se*se))
    post_mn2 <- post_var2*(
        (1.0/sigma0[2]^2)*mu0[2] + (n[2]/(2.0*se*se))*mn
    )
    post_se2 <- sqrt(post_var2)

    
    c(mean(pnorm(eoi,post_mn,post_se,lower.tail=lower.tail)),
      mean(pnorm(eoi,post_mn,post_se,lower.tail=lower.tail)>prth))
}



main <- function(){
    n_sims <- 10000
    mu <- c(-1.0,1.0)
    sigma <- 2.0
    n <- c(100,100)
    sims <- simulate_studies(n_sims,mu,sigma,n)

    eoi <- 2.0
    prth <- 0.9
    fit_flat(sims,n,eoi,prth)
    fit_inf(sims,n,eoi,prth,0.0,1.0)
    fit_inf(sims,n,eoi,prth,0.0,1000.0)
    fit_power(0.001,sims,n,eoi,prth,0.0,1.0)
    fit_power(1.0,sims,n,eoi,prth,0.0,1.0)
    
}


######################################################################
### post.R ends here
