#!/lrlhps/apps/R/qualified/R-3.4.0/bin/Rscript 
### ui.R --- 
## 
## Filename: ui.R
## Description: 
## Author: Michael D Sonksen
## Maintainer: 
## License: This is owned by Eli Lilly and Company
## Created: Mon Jan 14 11:19:47 2019 (-0500)
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
library(shiny)
library(ggplot2)
library(plotly)
library(shinydashboard)
library(rhandsontable)

## 
### Global Functions: 


rm(list=ls())
options(shiny.sanitize.errors=FALSE) ##informative error messages

## 
### Code: 

dashboardPage(
    dashboardHeader(
        title="Trial Simulator"
    ),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Scenarios", tabName="tab_scenarios"),
            menuItem("Priors", tabName="tab_priors"),
            menuItem("Design", tabName="tab_design"),
            menuItem("Simulations", tabName="tab_simulations")
        )
    ),
    dashboardBody(
        tabItems(
            tabItem(tabName="tab_scenarios",
                fluidRow(
                    box(title="Scenarios",rHandsontableOutput("hot_scenarios"))
                )
            ),
            tabItem(tabName="tab_priors",
                fluidRow(
                    box(title="Two Component Mixture Prior",
                        fluidRow(
                            column(6,numericInput("prior_mix_weight1",
                                                  "Mixture Prior Weight",
                                                  value=0.5,min=0.0,max=1.0,
                                                  step=0.1))
                        ),
                        fluidRow(
                        column(6,numericInput("prior_mix_mean1",
                                              "Mixture Prior Mean of Skeptical",
                                              value=0.0,step=0.1)),
                        column(6,numericInput("prior_mix_mean2",
                                              "Mixture Prior Mean of Informative",
                                              value=1.5,step=0.1))
                        ),
                        fluidRow(
                        column(6,numericInput("prior_mix_sd1",
                                              "Mixture Prior SD of Skeptical",
                                              value=1.0,min=0.0,step=0.1)),
                        column(6,numericInput("prior_mix_sd2",
                                              "Mixture Prior SD of Informative",
                                              value=3.0,min=0.0,step=0.1))
                        )
                    ),
                    box(plotOutput("prior_mix_plot"))
                ),
                fluidRow(
                    box(tabName="Power Priors",
                        numericInput("prior_power_a0",
                                     "Power Prior Discounting Parameter (a0)",
                                     value=0.5,min=0.0,max=1.0,step=0.1),
                        numericInput("prior_power_mean","Power Prior Mean",
                                     value=0.5,step=0.1),
                        numericInput("prior_power_sd", ## add plot without discount
                                     "Power Prior Standard Deviation",
                                     value=0.5,min=0.0,step=0.1)
                    ),
                    box(plotOutput("prior_power_plot"))
                )
            ),
            tabItem(tabName="tab_design", ## Add Case Examples
                fluidRow(
                    box(tabName="Design",
                        numericInput("design_n_control","Sample Size for Control",
                                     value=100,min=0,step=1),
                        numericInput("design_n_trt","Sample Size for TRT",
                                     value=100,min=0,step=1)
                    ),
                    box(tabName="Success Criterion",
                        "P(TRT-Control>EOI)>PRTH)",
                        numericInput("design_csf_eoi", "EOI",
                                     value=0.0,step=0.1),
                        numericInput("design_csf_prth", "PRTH",
                                     value=0.8,min=0.0,max=1.0,step=0.01),
                        uiOutput("design_csf")
                    )
                )
            ),
            tabItem(tabName="tab_simulations",
                fluidRow(
                    box(tabName="Simulation Settings",
                        numericInput("sim_seed",label="R random number seed",value=4212,step=1),
                        numericInput("sim_n",label="Number of Simulations",value=100,step=1),
                        actionButton("sim_go",label="Run Simulations!")
                    )
                ),
                fluidRow(
                    box(plotlyOutput("power_plot"))
                ),
                fluidRow(
                    box(numericInput("base_mn_control",label="Mean Control Group",value=0.0,step=0.1),
                        numericInput("base_sd_control",label="SD Control Group", value=1.0,step=0.1),
                        numericInput("base_sd_trt",label="SD Treatment Group", value=1.0,step=0.1),
                    numericInput("base_mn_trt_lb",label="Lower Bound of Mean Treatment Group",
                                 value=0.0,step=0.1),
                    numericInput("base_mn_trt_ub",label="Upper Bound of Mean Treatment Group",
                                 value=3.0,step=0.1),
                    numericInput("base_mn_trt_by",label="Step Size of Mean Treatment Group",
                                 value=0.1,step=0.1)),
                    box(plotlyOutput("power_plot_grid"))
                )
                
            )
        )
    ),title=''
)


######################################################################
### ui.R ends here
