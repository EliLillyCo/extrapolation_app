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


## 
### Code: 
shinyServer(function(input,output,clientData, session){

    sim_results <- list()
    sims <- NULL
    
    output$hot_scenarios <- renderRHandsontable({
        df <- data.frame("ScenarioName"=c("Null", "Null SD=2", "Expected",
                                  "Expected SD=2"),
                         "ControlMean"=c(0.0,0.0,1.0,1.0),
                         "ControlSD"=c(1.0,2.0,1.0,2.0),
                         "TRTMean"=c(2.0,3.0,2.0,3.0),
                         "TRTSD"=c(1.0,1.0,2.0,2.0))
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


    output$prior_mix_plot <- renderPlot({
        mn1 <- input$prior_mix_mean1
        sd1 <- input$prior_mix_sd1
        mn2 <- input$prior_mix_mean2
        sd2 <- input$prior_mix_sd2
        mix <- input$prior_mix_weight1
        dnorm_c <- function(x) mix*dnorm(x,mn1,sd1)+(1-mix)*dnorm(x,mn2,sd2)
        xvalues <- data.frame(x=c(min(mn1-3*sd1,mn2-3*sd2),
                                  max(mn1+3*sd1,mn2+3*sd2)))
        ggplot(xvalues,aes(x=xvalues))+
            stat_function(fun=dnorm_c, geom="area",fill="turquoise",alpha=.3)+
            xlim(c(min(mn1-4*sd1,mn2-4*sd2),max(mn2+4*sd2,mn1+4*sd1)))+
            ##geom_vline(xintercept=input$eoi1,colour="red")+
            labs(x="Variable", y="Density")
    })

    output$prior_power_plot <- renderPlot({
        mn <- input$prior_power_mean
        sd <- input$prior_power_sd
        a0 <- input$prior_power_a0
        dnorm_c <- function(x) dnorm(x,mn,sd)^(a0)
        xvalues <- data.frame(x=c(mn-10*sd,mn-10*sd,
                                  mn+10*sd,mn+10*sd))
        ggplot(xvalues,aes(x=xvalues))+
            stat_function(fun=dnorm_c, geom="area",fill="orange",alpha=.3)+
            xlim(c(min(mn-11*sd,mn-11*sd),max(mn+11*sd,mn+11*sd)))+
            ##geom_vline(xintercept=input$eoi1,colour="red")+
            labs(x="Variable", y="Density")
    })
    output$design_csf <- renderUI({
        withMathJax(paste0("P(TRT-Control > ",input$design_csf_eoi,")>",input$design_csf_prth))
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
            
            lapply(1:dim(scenarios)[1],function(x)
                simulate_studies(n_sims,
                                 mu=c(scenarios$ControlMean[x],
                                             scenarios$TRTMean[x]),
                                 sigma=c(scenarios$ControlSD[x],
                                         scenarios$TRTSD[x]),n))
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
        lapply(mn1, function(x)
            simulate_studies(n_sims,
                             mu=c(mn0,x),
                             sigma=c(sd0,sd1),n)
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
        lapply(sms,function(x)
            list("flat"=fit_flat(x,n,eoi,prth,lower.tail=FALSE)[2],
                 "inf"=fit_inf(x,n,eoi,prth,mu0,sigma0,lower.tail=FALSE)[2],
                 "power"=fit_power(a0,x,n,eoi,prth,mu0,sigma0,lower.tail=FALSE)[2]
                 )
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
        lapply(sms,function(x)
            list("flat"=fit_flat(x,n,eoi,prth,lower.tail=FALSE)[2],
                 "inf"=fit_inf(x,n,eoi,prth,mu0,sigma0,lower.tail=FALSE)[2],
                 "power"=fit_power(a0,x,n,eoi,prth,mu0,sigma0,lower.tail=FALSE)[2]
                 )
            )
        
    })

    prbs_df <- reactive({
        prb <- prbs()
        scenarios <- hot_to_r(input$hot_scenarios)
        scen_names <- scenarios$ScenarioName        
        df <- data.frame("prcsf"=as.numeric(unlist(prb)),
                         "scenario"=as.character(rep(scen_names,each=3)),
                         "prior"=as.character(rep(c("flat","informative","power"),length(scen_names))))
        df
    })


    
    prbs_df <- reactive({
        prb <- prbs_grid()
        scen_names <- seq(input$base_mn_trt_lb,
                   input$base_mn_trt_ub,
                   input$base_mn_trt_by)
        df <- data.frame("prcsf"=as.numeric(unlist(prb)),
                         "scenario"=as.numeric(rep(scen_names,each=3)),
                         "prior"=as.character(rep(c("flat","informative","power"),
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
        df <-  prbs_df()
        ggplot(data = df, aes(x = scenario, y = prcsf,colour=prior)) +
            geom_point()+geom_line()        
    })
    
})


######################################################################
### server.R ends here
