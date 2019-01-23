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
        df <- data.frame("Scenario Name"=c("Null", "Null SD=2", "Expected",
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
    observeEvent(input$sim_go,
    {
        set.seed(input$sim_seed)
        n_sims <- input$sim_n

        ##generate sims
        
        
    })
})


######################################################################
### server.R ends here
