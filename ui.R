# load packages ----
library(shiny)
library(shinycssloaders)
library(ggplot2)
library(plotly)
library(shinydashboard)
library(rhandsontable)
library(grid)

# define options ----
options(shiny.sanitize.errors=FALSE) ##informative error messages

# initialize user interface ----
dashboardPage(
    dashboardHeader(
        title="TrialSimulator"
    ),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Scenarios and Design", tabName="tab_scenarios"),
            menuItem("Priors", tabName="tab_priors"),
            menuItem("Simulations", tabName="tab_simulations")
        )
    ),
    dashboardBody(
        tabItems(
          tabItem(
            tabName="tab_scenarios",
              fluidRow(
                  box(tabName="Design",
                      title="Sample Sizes",
                      solidHeader=TRUE,
                      collapsible=TRUE,
                      status="primary",
                      numericInput("design_n_control","Sample Size for Control",
                                   value=150,min=0,max=10000,step=1),
                      numericInput("design_n_trt","Sample Size for TRT",
                                   value=150,min=0,max=10000,step=1),
                      uiOutput("model")
                      ),
                box(tabName="Success Criterion",title="Success Criterion",
                    solidHeader=TRUE,
                    collapsible=TRUE,status="primary",
                    "P(TRT-Control<EOI)>PrTH) \n",
                    br(),
                    "EOI = Effect of Interest\n",
                    br(),
                    "PrTH = Probability Threshold",
                    numericInput("design_csf_eoi", "EOI",
                                 value=0.0,step=0.1),
                    numericInput("design_csf_prth", "PRTH",
                                 value=0.99,min=0.0,max=1.0,step=0.01),
                    uiOutput("design_csf")
                )
              ),
              fluidRow(
                box(title="Scenarios",rHandsontableOutput("hot_scenarios"),
                    solidHeader=TRUE,
                    collapsible=TRUE,status="primary")
              )
            ),
            tabItem(
              tabName="tab_priors",
                fluidRow(
                  box(title="Two Component Mixture Prior",solidHeader=TRUE,
                    collapsible=TRUE,status="primary",
                    uiOutput("mixture_prior"),
                    fluidRow(
                      column(6,numericInput("prior_mix_weight1",
                                            "Mixture Prior Weight (p)",
                                            value=0.5,min=0.0,max=1.0,
                                            step=0.1))
                    ),
                    fluidRow(
                    column(6,numericInput("prior_mix_mean1",
                                          "Mixture Prior Mean of Skeptical (m0)",
                                          value=0.0,step=0.1,min=-100,max=100)),
                    column(6,numericInput("prior_mix_mean2",
                                          "Mixture Prior Mean of Informative (m1)",
                                          value=-2.0,step=0.1,min=-100,max=100))
                    ),
                    fluidRow(
                    column(6,numericInput("prior_mix_sd1",
                                          "Mixture Prior SD of Skeptical (s0)",
                                          value=3.0,min=0.0,step=0.1,max=100)),
                    column(6,numericInput("prior_mix_sd2",
                                          "Mixture Prior SD of Informative (s1)",
                                          value=0.275,min=0.0,step=0.1,max=100))
                    )
                  ),
                  box(withSpinner(plotOutput("prior_mix_plot")),
                      title="Mixture Prior Density Plot",solidHeader=TRUE,
                      collapsible=TRUE,status="primary")
                ),
                fluidRow(
                  box(tabName="Power Priors",title="Power Prior",
                      solidHeader=TRUE,
                      collapsible=TRUE,status="primary",
                      uiOutput("power_prior"),
                      numericInput("prior_power_a0",
                                   "Power Prior Discounting Parameter (a0)",
                                   value=0.04,min=0.0,max=1.0,step=0.1),
                      numericInput("prior_power_mean","Power Prior Mean (m)",
                                   value=-1.97,step=0.1,min=-1000,max=1000),
                      numericInput("prior_power_sd", ## add plot without discount
                                   "Power Prior Standard Deviation (s)",
                                   value=0.275,min=0.0,step=0.1,max=1000)
                  ),
                  box(withSpinner(plotOutput("prior_power_plot")),
                      title="Power Prior Density Plot",solidHeader=TRUE,
                      collapsible=TRUE,status="primary")
                ),
                fluidRow(
                  box(tabName="Hypothetical Data",title="Example Observed Data",
                      solidHeader=TRUE,
                      collapsible=TRUE,status="primary",
                      numericInput("obs_control_mean",
                                   label="Observed Control Mean",value=0.0,
                                   step=0.1,min=-100,max=100),
                      numericInput("obs_control_sd",
                                   label="Observed Control SD",value=3.5,
                                   step=0.1,
                                   min=0, max=1000),
                      numericInput("obs_trt_mean",
                                   label="Observed TRT Mean",
                                   value=-2.0,step=0.1,min=-100,max=100),
                      numericInput("obs_trt_sd",label="Observed TRT SD",
                                   value=3.5,step=0.1,min=0,max=1000),
                      numericInput("pwr_y",label="Power Y-position",
                                   value=1.5, step=0.1,
                                   min=-1000,max=1000),
                      numericInput("pwr_x",label="Power X-position",
                                   value=0.0, step=0.1,
                                   min=-1000,max=1000)
                      
                    ),
                    box(
                      htmlOutput("prior2post_text"),
                        withSpinner(plotlyOutput("prior2post_mix")),
                        title="Prior to Posterior Plots",solidHeader=TRUE,
                        collapsible=TRUE,status="primary"
                    )
                )
            ),
            tabItem(tabName="tab_simulations",
              fluidRow(
                box(tabName="Simulation Settings",solidHeader=TRUE,
                    collapsible=TRUE,status="primary",title="Simulation Setting",
                    numericInput("sim_seed",label="R random number seed",value=4212,step=1),
                    numericInput("sim_n",label="Number of Simulations",value=100,step=1),
                    actionButton("sim_go",label="Run Simulations!")
                ),
                box(withSpinner(plotlyOutput("power_plot")),title="Power By Scenario",solidHeader=TRUE,
                    collapsible=TRUE,status="primary")
              ),
              fluidRow(
                box(numericInput("base_mn_control",
                                 label="Mean Control Group",value=0.0,step=0.1,
                                 min=-1000,max=1000),
                    numericInput("base_sd_control",label="SD Control Group",
                                 value=6,step=0.1,min=-1000,max=1000),
                    numericInput("base_sd_trt",label="SD Treatment Group",
                                 value=6,step=0.1,min=-1000,max=1000),
                    numericInput("base_mn_trt_lb",
                                 label="Lower Bound of Mean Treatment Group",
                                 value=-3.0,step=0.1,min=-1000,max=1000),
                    numericInput("base_mn_trt_ub",
                                 label="Upper Bound of Mean Treatment Group",
                             value=0.0,step=0.1,min=-1000,max=1000),
                    numericInput("base_mn_trt_by",
                                 label="Step Size of Mean Treatment Group",
                                 value=0.1,step=0.1,min=-1000,max=1000),
                    actionButton("sim_grid_go",label="Run Simulations!"),
                    title="Power Curve Settings",solidHeader=TRUE,
                    collapsible=TRUE,status="primary"),
                box(withSpinner(plotlyOutput("power_plot_grid")),title="Power Curve",solidHeader=TRUE,
                    collapsible=TRUE,status="primary")
              )
            )
        )
    ),title=''
)
