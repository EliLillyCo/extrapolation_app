#!/lrlhps/apps/R/qualified/R-3.4.0/bin/Rscript 
### server.R --- 
## 
## Filename: server.R
## Description: 
## Author: Michael D Sonksen
## Maintainer: 
## License: This is owned by Eli Lilly and Company
## Created: Mon Jan 14 11:20:12 2019 (-0500)
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
library(ggplot2)
library(plotly)
library(htmlwidgets)
## 
### Global Functions: 
source("post.R")

### global variables
n_cores <- 6


## 
### Code: 
shinyServer(function(input,output,clientData, session){

    sim_results <- list()
    sims <- NULL
    
    output$hot_scenarios <- renderRHandsontable({
        df <- data.frame("ScenarioName"=c("Null", "Null SD=2", "Expected",
                                  "Expected SD=2"),
                         "ControlMean"=c(0.0,0.0,0.0,1.0),
                         "ControlSD"=c(1.0,2.0,1.0,2.0),
                         "TRTMean"=c(0.0,0.0,-2.0,-1.0),
                         "TRTSD"=c(1.0,2.0,1.0,2.0))
        rhandsontable(df, useTypes=FALSE,height=300*2,
                    search=TRUE) %>%
        hot_table(highlightCol=TRUE,highlightRow=TRUE) %>%
            hot_context_menu(allowRowEdit=TRUE, allowColEdit=TRUE,
                             allowComments=TRUE,
                             customOpts = list(
                             csv = list(name = "Download to CSV",
                                      callback = htmlwidgets::JS(
                                        "function (key, options) {
                                        var csv = csvString(this, sep=',', dec='.');
                                        
                                        var link = document.createElement('a');
                                        link.setAttribute('href', 'data:text/plain;charset=utf-8,' +
                                        encodeURIComponent(csv));
                                        link.setAttribute('download', 'data.csv');
                                        
                                        document.body.appendChild(link);
                                        link.click();
                                        document.body.removeChild(link);
  }")))) %>%
        hot_cols(columnSorting=TRUE,manualColumnMove=TRUE,
                 manualColumnResize=TRUE)

    })

    output$prior2post_mix <- renderPlotly({
        n <- c(input$design_n_control,input$design_n_trt)

        mn1 <- input$prior_mix_mean1
        sd1 <- input$prior_mix_sd1
        mn2 <- input$prior_mix_mean2
        sd2 <- input$prior_mix_sd2
        mix <- input$prior_mix_weight1
        
        dnorm_c <- function(x) mix*dnorm(x,mn1,sd1)+(1-mix)*dnorm(x,mn2,sd2)
        xvalues <- data.frame(x=c(min(mn1-3*sd1,mn2-3*sd2),
                                  max(mn1+3*sd1,mn2+3*sd2)))
        ##fake data
        
        obs_mn1 <- input$obs_control_mean
        obs_sd1 <- input$obs_control_sd
        obs_mn2 <- input$obs_trt_mean
        obs_sd2 <- input$obs_trt_sd
        
        s2 <- c(obs_sd1^2,obs_sd2^2)
        ybar <-c(obs_mn1,obs_mn2)
        
        m0 <- marginal(ybar,s2,c(0,mn1),
                       c(100,sd1),n,sigma=sqrt(s2),log=TRUE)
        m1 <- marginal(ybar,s2,c(0,mn2),
                       c(100,sd2),n,sigma=sqrt(s2),log=TRUE)
        max_m <- max(m0,m1)
        m0 <- exp(m0-max_m)
        m1 <- exp(m1-max_m)
        ybar <- c(ybar[1], ybar[2]-ybar[1])
        p <- mix*m0/(mix*m0+(1-mix)*m1)
        p <- c(p,1-p)
        
        post_var <- c(1.0/(1.0/sd1^2+n[2]/(2.0*s2[2])),
                      1.0/(1.0/sd2^2+n[2]/(2.0*s2[2])))
        
        post_mn <- c(post_var[1]*((1.0/sd1^2)*mn1 +
                                  (n[2]/(2.0*s2[2]))*ybar[2]),
                     post_var[2]*((1.0/sd2^2)*mn2 +
                                  (n[2]/(2.0*s2[2]))*ybar[2]))
        
        dnorm_pc <- function(x) p[1]*dnorm(x,post_mn[1],sqrt(post_var[1]))+
                                    p[2]*dnorm(x,post_mn[2],sqrt(post_var[2]))
        pp_a0 <- input$prior_power_a0
        pp_mn <- input$prior_power_mean
        pp_sd <- input$prior_power_sd/sqrt(pp_a0)
        pp_post_var <- 1.0/(1.0/pp_sd^2+n[2]/(2.0*s2[2]))
        pp_post_mn <- pp_post_var[1]*((1.0/pp_sd^2)*pp_mn +
                                      (n[2]/(2.0*s2[2]))*ybar[2])
        dnorm_pwrp <- function(x) dnorm(x,pp_mn,pp_sd/sqrt(pp_a0))
        dnorm_pwr <- function(x) dnorm(x,pp_post_mn,sqrt(pp_post_var/pp_a0))

        txt <- paste0("Mix: P(TRT-Control<0|data)=",
                      round(p[1]*pnorm(input$design_csf_eoi,post_mn[1],
                                       sqrt(post_var[1]))+
                            p[2]*pnorm(input$design_csf_eoi,post_mn[2],
                                       sqrt(post_var[2])),3),"\n",
                      "Power: P(TRT-Control<0|data)=",
                      round(pnorm(input$design_csf_eoi,pp_post_mn,
                                       sqrt(pp_post_var)),3))
        grob <- grobTree(textGrob(txt, x=0.5,y=0.5, hjust=0,
                                  gp=gpar(col="black", fontsize=13,
                                          fontface="italic")))
        
        ggplot(xvalues,aes(x=x))+
            stat_function(fun=dnorm_c, geom="area",
                          aes(fill="MixturePrior"),alpha=.3)+
            stat_function(fun=dnorm_pc, geom="area",
                          aes(fill="MixturePosterior"),alpha=.3)+
            stat_function(fun=dnorm_pwrp, geom="area",
                          aes(fill="PowerPrior"),alpha=.3)+
            stat_function(fun=dnorm_pwr, geom="area",
                          aes(fill="PowerPosterior"),alpha=.3)+
            xlim(c(min(mn1-4*sd1,mn2-4*sd2,
                       post_mn[1]-4*sqrt(post_var[1]),
                       post_mn[2]-4*sqrt(post_var[2])),
                   max(mn2+4*sd2,mn1+4*sd1,
                       post_mn[1]+4*sqrt(post_var[1]),
                       post_mn[2]+4*sqrt(post_var[2]))))+
            geom_vline(xintercept=ybar[2],colour="red")+
            labs(x="Treatment Difference", y="Density")+
            scale_fill_manual("Distribution",
                              values=c("yellow","turquoise","purple","orange"))+
##            annotation_custom(grob)
        annotate("text",x=0,y=3,label=txt)
        
    })

    output$prior_mix_plot <- renderPlot({
        mn1 <- input$prior_mix_mean1
        sd1 <- input$prior_mix_sd1
        mn2 <- input$prior_mix_mean2
        sd2 <- input$prior_mix_sd2
        mix <- input$prior_mix_weight1
        dnorm_c <- function(x) mix*dnorm(x,mn1,sd1)+(1-mix)*dnorm(x,mn2,sd2)
        xvalues <- data.frame(x=c(min(mn1-3*sd1,mn2-3*sd2),
                                  max(mn1+3*sd1,mn2+3*sd2)))
        txt <- paste0("P(TRT-Control<0)=",round(mix*pnorm(input$design_csf_eoi,mn1,sd1)+(1-mix)*pnorm(input$design_csf_eoi,mn2,sd2),3))
        grob <- grobTree(textGrob(txt, x=0.6,y=0.95, hjust=0,
  gp=gpar(col="black", fontsize=13, fontface="italic")))
        
        ggplot(xvalues,aes(x=x))+
            stat_function(fun=dnorm_c, geom="area",fill="turquoise",alpha=.3)+
            xlim(c(min(mn1-4*sd1,mn2-4*sd2),max(mn2+4*sd2,mn1+4*sd1)))+
            ##geom_vline(xintercept=input$eoi1,colour="red")+
            labs(x="Treatment Difference", y="Density")+
            annotation_custom(grob)
        
    })

    output$prior_power_plot <- renderPlot({
        mn <- input$prior_power_mean
        sd <- input$prior_power_sd
        a0 <- input$prior_power_a0
        dnorm_c <- function(x) dnorm(x,mn,sd)^(a0)
        xvalues <- data.frame(x=c(mn-10*sd,mn-10*sd,
                                  mn+10*sd,mn+10*sd))
        txt <- paste0("P(TRT-Control<0)=",round(pnorm(input$design_csf_eoi,mn,sd/sqrt(a0)),3))
        grob <- grobTree(textGrob(txt, x=0.6,y=0.95, hjust=0,
  gp=gpar(col="black", fontsize=13, fontface="italic")))
        
        ggplot(xvalues,aes(x=x))+
            stat_function(fun=dnorm_c, geom="area",fill="orange",alpha=.3)+
            xlim(c(min(mn-11*sd,mn-11*sd),max(mn+11*sd,mn+11*sd)))+
            ##geom_vline(xintercept=input$eoi1,colour="red")+
            labs(x="Treatment Difference", y="Density")+
            annotation_custom(grob)
    })
    output$design_csf <- renderUI({
        withMathJax(paste0("P(TRT-Control < ",input$design_csf_eoi,")>",input$design_csf_prth))
    })

    ##Simulations
    sims <- reactive({
        input$sim_go
        set.seed(input$sim_seed)
        n_sims <- input$sim_n
        n <- c(input$design_n_control, input$design_n_trt)
        hot <- isolate(input$hot_scenarios)
        if(!is.null(hot)){
            scenarios <- hot_to_r(input$hot_scenarios)
            
            parallel::mclapply(1:dim(scenarios)[1],function(x)
                simulate_studies(n_sims,
                                 mu=c(scenarios$ControlMean[x],
                                             scenarios$TRTMean[x]),
                                 sigma=c(scenarios$ControlSD[x],
                                         scenarios$TRTSD[x]),n),
                mc.cores = n_cores)
        }else{
            NULL
        }
        
    })

    sims_grid <- reactive({
        input$sim_go
        set.seed(input$sim_seed)
        n_sims <- input$sim_n
        n <- c(input$design_n_control, input$design_n_trt)
        sd0 <- input$base_sd_control
        mn0 <- input$base_mn_control
        sd1 <- input$base_sd_trt
        mn1 <- seq(input$base_mn_trt_lb,
                   input$base_mn_trt_ub,
                   input$base_mn_trt_by)
        parallel::mclapply(mn1, function(x)
            simulate_studies(n_sims,
                             mu=c(mn0,x),
                             sigma=c(sd0,sd1),n),
            mc.cores = n_cores
            )
    })
    
    
    prbs <- reactive({
        sms <- sims()
        eoi <- input$design_csf_eoi
        prth <- input$design_csf_prth
        n <- c(input$design_n_control, input$design_n_trt)
        mu0 <- input$prior_power_mean
        sigma0 <- input$prior_power_sd
        a0 <- input$prior_power_a0
        p0 <- input$prior_mix_weight1
        mu0_mix <- c(0,input$prior_mix_mean1)
        mu1_mix <- c(0,input$prior_mix_mean2)
        sigma0_mix <- c(100,input$prior_mix_sd1)
        sigma1_mix <- c(100,input$prior_mix_sd2)
        parallel::mclapply(sms,function(x)
            list("flat"=fit_flat(x,n,eoi,prth,lower.tail=FALSE)[2],
                 "inf"=fit_inf(x,n,eoi,prth,mu0,sigma0,lower.tail=FALSE)[2],
                 "power"=fit_power(a0,x,n,eoi,prth,mu0,sigma0,lower.tail=FALSE)[2],
                 "mix"=fit_mix(p0,x,n,eoi,prth,mu0_mix,sigma0_mix,mu1_mix,
                               sigma1_mix,lower.tail=FALSE)[2]
                 ),
            mc.cores = n_cores
            )
        
    })


    prbs_grid <- reactive({
        sms <- sims_grid()
        eoi <- input$design_csf_eoi
        prth <- input$design_csf_prth
        n <- c(input$design_n_control, input$design_n_trt)
        mu0 <- input$prior_power_mean
        sigma0 <- input$prior_power_sd
        a0 <- input$prior_power_a0
        p0 <- input$prior_mix_weight1
        mu0_mix <- c(0,input$prior_mix_mean1)
        mu1_mix <- c(0,input$prior_mix_mean2)
        sigma0_mix <- c(100,input$prior_mix_sd1)
        sigma1_mix <- c(100,input$prior_mix_sd2)
        parallel::mclapply(sms,function(x)
            list("flat"=fit_flat(x,n,eoi,prth,lower.tail=FALSE)[2],
                 "inf"=fit_inf(x,n,eoi,prth,mu0,sigma0,lower.tail=FALSE)[2],
                 "power"=fit_power(a0,x,n,eoi,prth,mu0,sigma0,lower.tail=FALSE)[2],
                 "mix"=fit_mix(p0,x,n,eoi,prth,mu0_mix,sigma0_mix,mu1_mix,sigma1_mix,lower.tail=FALSE)[2]
                 ),
            mc.cores = n_cores
            )
        
    })

    prbs_df <- reactive({
        prb <- prbs()
        scenarios <- hot_to_r(input$hot_scenarios)
        scen_names <- scenarios$ScenarioName        
        df <- data.frame("prcsf"=as.numeric(unlist(prb)),
                         "scenario"=as.character(rep(scen_names,each=4)),
                         "prior"=as.character(rep(c("flat","informative","power","mix"),
                                                  length(scen_names))))
        df
    })


    
    prbs_df_grid<- reactive({
        prb <- prbs_grid()
        scen_names <- seq(input$base_mn_trt_lb,
                   input$base_mn_trt_ub,
                   input$base_mn_trt_by)
        df <- data.frame("prcsf"=as.numeric(unlist(prb)),
                         "scenario"=as.numeric(rep(scen_names,each=4)),
                         "prior"=as.character(rep(c("flat","informative","power", "mix"),
                                                  length(scen_names))))
        df
    })
    
    output$power_plot <- renderPlotly({
        df <-  prbs_df()
        ggplot(data = df, aes(x = scenario, y = prcsf,colour=prior)) +
            geom_point() +
            geom_line()
        
    })

    output$power_plot_grid <- renderPlotly({
        df <-  prbs_df_grid()
        ggplot(data = df, aes(x = scenario, y = prcsf,colour=prior)) +
            geom_point()+geom_line()        
    })
    
})


######################################################################
### server.R ends here
