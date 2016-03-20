################################################################################
# User Interface for Baker (UI panel file)
#
# Designed for PERCH data in the current version. Please do not hesitate to
# contact the author(s) below for software design and modifications.
#
# Zhenke Wu
# zhenkewu@gmail.com
# 1st version: March 17, 2016
################################################################################

library(shinyFiles)
library(shinydashboard)

mywellPanel <- function(...){
  wellPanel(...)
}

sidewellPanel <- function(...){
  wellPanel(style="background-color: #222D32;",...)
}

header <- dashboardHeader(
  title = "baker: Bayesian Analytic Kit for Etiology Research",
  titleWidth = 450,
  dropdownMenuOutput("messageMenu"),
  
  # dropdownMenu(type = "notifications",
  #              notificationItem(
  #                text = "5 new users today",
  #                icon("users")
  #              ),
  #              notificationItem(
  #                text = "12 items delivered",
  #                icon("truck"),
  #                status = "success"
  #              ),
  #              notificationItem(
  #                text = "Server load at 86%",
  #                icon = icon("exclamation-triangle"),
  #                status = "warning"
  #              )
  # ),
  dropdownMenu(type = "tasks", badgeStatus = "success",
               taskItem(value = 90, color = "green",
                        "Documentation"
               ),
               taskItem(value = 17, color = "aqua",
                        "Project X"
               ),
               taskItem(value = 75, color = "yellow",
                        "Server deployment"
               ),
               taskItem(value = 80, color = "red",
                        "Overall project"
               )
  )
)

sidebar <- dashboardSidebar(
  width = 350,
  # sidebarUserPanel(
  #   "Zhenke Wu","baker",
  #   "PERCH_logo.svg"
  #   ),
  
  sidebarMenu(id="main_method",
              menuItem("Input Data",tabName="Input",icon=icon("table")),
              menuItem("Explore Data",tabName="Explore",icon=icon("compass")),
              menuItem( "Input Prior",tabName="Prior",icon=icon("area-chart")),
              menuItem("Specify Model Likelihood",tabName="Likelihood",icon=icon("key")),
              menuItem("Fit Model",tabName="Analyze",icon=icon("pie-chart")),
              menuItem("Check Model",tabName="Check",icon=icon("stethoscope")),
              menuItem("Visualizations",tabName="Visualize",icon=icon("eye")),
              menuItem("Save Results",tabName="Download",icon=icon("download")),
              menuItem("Messages",tabName="Messages",icon=icon("history")),
              menuItem("About" , tabName="About",icon=icon("info-circle"))#,
              #menuItem("Email Developers", tabName="Email",icon=icon("reply"))
  ),
  
  conditionalPanel(condition = "input.main_method=='Input'",
                   tagList(
                     sidewellPanel(
                       h5("Input Lab Test Data"),
                       tags$p(HTML('Rows: Individuals')),
                       tags$p(HTML('Columns: Lab tests - {Pathogen}_{Specimen}{Test}, e.g., PNEU_NPPCR,
                                         Other info: Case/control dfinition variables and covariates')),
                       fileInput('input_lab_data', 'Choose File (.csv)',accept=".csv")
                     ),
                     
                     
                     sidewellPanel(
                       sidewellPanel(uiOutput("siteUI")),
                       sidewellPanel(
                         selectInput("casedef", "Case Definition Variable",
                                     list("Chest X-Ray Positive (CXRFINCAT_5)"="CXRFINCAT_5")),
                         numericInput("casedef_val", "Value",
                                      1,min=1,max=2)
                       ),
                       sidewellPanel(
                         selectInput("ctrldef", "Control Definition Variable",
                                     list("All Controls (CASECONT)"="CASECONT")),
                         numericInput("ctrldef_val", "Value",
                                      2,min=1,max=2)
                       ),
                       
                       sidewellPanel(
                         h5("Causes and Measurements"),
                         actionButton("cause_button","Choose",
                                      icon("bug"), 
                                      style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                       ),
                       
                       sidewellPanel(
                         uiOutput("X_extra_UI"),
                         actionButton("selectall_X_extra","Un/Select All")
                       ),
                       
                       actionButton("refresh_button","Initialize",
                                    icon("sliders"), 
                                    style="color: #fff; background-color: #337ab7; border-color: #2e6da4;
                                    display: inline-block")
                     )
                   )
                   
  ),
  
  #
  # Exploratory data analysis:
  #
  conditionalPanel(condition="input.main_method=='Explore'",
                   sidewellPanel(
                     radioButtons("eda_options","EDA Options",choices=c(
                       "Data Summary"="data_structure",
                       "Test Correlations"="pairLOR",
                       "Seasonal Trend"="season"),
                       selected="data_structure")
                   ),
                   conditionalPanel(condition="input.eda_options=='pairLOR'",
                                    actionButton("plot_pairLOR_button", "Plot")
                   ),
                   conditionalPanel(condition="input.eda_options=='season'",
                                    sidewellPanel(
                                      uiOutput("choose_patho_for_season_UI"),
                                      actionButton("selectall_choose_patho_for_season","Select All"),
                                      actionButton("visualize_season_button","Plot")
                                    )
                   )
  ),
  
  
  #
  # Specify prior; sidebar:
  #
  conditionalPanel(condition="input.main_method=='Prior'",
                   sidewellPanel(
                     sidewellPanel(
                       uiOutput("Eti_prior_UI")
                     ),
                     sidewellPanel(
                       uiOutput("TPR_prior_BrS_info_sidebarUI"),
                       uiOutput("TPR_prior_SS_info_sidebarUI")
                     )
                   )),
  
  
  
  #
  # specify measurement likelihoods:
  #
  
  conditionalPanel(condition="input.main_method=='Likelihood'",
                   tagList(
                     sidewellPanel(
                       uiOutput("use_measurement_UI")
                     ),
                     
                     uiOutput("subclass_UI"),
                     
                     sidewellPanel(
                       h5("Formula for Etiology Regression"),
                       uiOutput("Eti_formula_UI")
                     ),
                     
                     sidewellPanel(
                       h5("Formula for False Positive Rate Regression"),
                       uiOutput("FPR_formula_UI")
                     )
                   )
  ),
  
  
  conditionalPanel(condition = "input.main_method=='Analyze'",
                   sidewellPanel(
                     tags$h4('The file path for writing and reading fitted results'),
                     tags$p(HTML('JAGS requires a working directory to write the data, initial values and model files. 
                                     The <code> baker </code> package will read this file path for posterior sample processing.' )),
                     #shinyDirButton('res_directory', 'Select Working Directory', 'Please select a folder'),
                     
                     tagList(singleton(tags$head(tags$script(src = "sF/shinyFiles.js"), 
                                                 tags$link(rel = "stylesheet", type = "text/css", href = "sF/styles.css"), 
                                                 tags$link(rel = "stylesheet", type = "text/css", href = "sF/fileIcons.css"))), 
                             tags$button(id = 'res_directory', type = "button", 
                                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                                         class = paste(c("shinyDirectories btn",paste0("btn-", "default"), class), collapse = " "), 
                                         `data-title` = 'Please select a folder', list(icon("folder-open-o"),as.character('Choose')))
                     ),
                     tags$p(),
                     
                     verbatimTextOutput('display_res_directory_path'),
                     textInput("default_res_directory","Default working directory (end with '/')",
                               value="~/Downloads/test_baker/")
                   ),
                   sidewellPanel(
                     tags$h4('Specify options for MCMC sampling'),
                     uiOutput("specify_mcmc_UI"),
                     actionButton("fit_button","Fit Model", icon("play-circle-o"), 
                                  style="color: #fff; background-color: #337ab7; border-color: #2e6da4;
                                    display: block")
                   )
  ),
  
  
  
  conditionalPanel(condition="input.main_method=='Visualize'",
                   sidewellPanel(
                     radioButtons("visualize_options","Figures",choices = c("Pie Standardization"="pie_standardization",
                                                                            "Individual Prediction"="individual_prediction",
                                                                            "Grouped Etiologies"="group_etiology"))
                   ),
                   conditionalPanel(condition="input.visualize_options=='group_etiology'",
                                    sidewellPanel(
                                      h5("Input Pathogen Taxomony Data"),
                                      tags$p(HTML('Rows: Pathogens')),
                                      tags$p(HTML('Columns: <code> pathogen </code> (Abbrevations, e.g., PNEU) and <code> pathogen_type </code> (V or B, or F)')),
                                      fileInput('input_taxo_data', 'Choose File (.csv)',accept=".csv"),
                                      tags$p(),
                                      actionButton("plot_group_etiology_button","Plot",
                                                   style="float:center")
                                    )
                   )
  )
  
)


body <- dashboardBody(
  tags$head(tags$style(HTML('
        .skin-blue .main-header .logo {
          background-color: #3c8dbc;
        }
        .skin-blue .main-header .logo:hover {
          background-color: #3c8dbc;
        }
      '))),
  
  uiOutput("baker_statusUI"),
  
  tabItems(
    #
    # data input:
    #
    tabItem(tabName = "Input",
            
            tabsetPanel(
              
              tabPanel("Specify Causes / Add Measurements", 
                       verbatimTextOutput("frontpage_upload_lab_data_warning"),
                       conditionalPanel(condition="input.cause_button > 0",
                                        splitLayout(
                                          box(width='100%',
                                              title = "Potential Causes of Diseases", status = "primary", 
                                              solidHeader = TRUE,collapsible = TRUE,
                                              uiOutput("cause_listUI"),
                                              actionButton("selectall_cause_list","Un/Select All")
                                          ),
                                          box(width='100%',
                                              title = "Bronze-Standard Measurements", status = "info", 
                                              solidHeader = TRUE,collapsible = TRUE,
                                              tags$p(HTML('Imperfect sensitivity and imperfect specificity')),
                                              tags$p(HTML('Available for both cases and controls')),
                                              mywellPanel(
                                                uiOutput("BrS_object_NPPCR_UI"),
                                                actionButton("selectall_patho_BrS_NPPCR","Un/Select All")
                                              ),
                                              mywellPanel(
                                                uiOutput("BrS_object_WBPCR_UI"),
                                                actionButton("selectall_patho_BrS_WBPCR","Un/Select All")
                                              ),
                                              mywellPanel(
                                                uiOutput("BrS_object_NPCX_VT13_UI"),
                                                actionButton("selectall_patho_BrS_NPCX_VT13","Un/Select All")
                                              )
                                              
                                          ),
                                          box(width='100%',
                                              title = "Silver-Standard Measurements", status = "warning", 
                                              solidHeader = TRUE,collapsible = TRUE,
                                              tags$p(HTML('Perfect Specificity; Imperfect Sensitivity')),
                                              tags$p(HTML('Available for cases only')),
                                              mywellPanel(
                                                uiOutput("SS_object_BCX_UI"),
                                                actionButton("selectall_patho_SS_BCX","Un/Select All")
                                              ),
                                              mywellPanel(
                                                uiOutput("SS_object_LA_ADJ_UI"),
                                                actionButton("selectall_patho_SS_LA_ADJ","Un/Select All")
                                              ),
                                              mywellPanel(
                                                uiOutput("SS_object_PF_ADJ_UI"),
                                                actionButton("selectall_patho_SS_PF_ADJ","Un/Select All")
                                              )
                                          )
                                        )
                       )
                       
              ) 
            )
            
    ),
    
    #
    # Exploratory data analysis:
    #
    tabItem(tabName="Explore",
            uiOutput("eda_UI")
    ),
    
    #
    # Specify prior:
    #
    tabItem(tabName="Prior",
            tabsetPanel(
              tabPanel("True Positive Rate Priors",
                       splitLayout(
                         uiOutput("TPR_prior_BrS_UI"),
                         uiOutput("TPR_prior_SS_UI")
                       )
              )
            )
    ),
    
    #
    # Specify likelihood:
    #
    tabItem(tabName="Likelihood",
            tabsetPanel(
              tabPanel("Summary of Model Options",
                       mywellPanel(
                         verbatimTextOutput("print_model_options")
                       )
              )
            )
    ),
    
    # Fit model:
    tabItem(tabName = "Analyze",
            tabsetPanel(
              tabPanel("MCMC Options",
                       mywellPanel(
                         h5("Summary of Markov Chain Monte Carlo Options:"),
                         verbatimTextOutput("print_mcmc_options")
                       )
              ),
              tabPanel("Panels Plot",
                       uiOutput("download_panels_UI"),
                       uiOutput("panels_UI")
              )
            )),
    
    
    # visualize:
    tabItem(tabName="Visualize",
            uiOutput("visualize_UI")
    ),
    
    
    #########################################################################
    #########################################################################
    #########################################################################
    # about page:
    tabItem(tabName="About",
            mywellPanel(
              p('Current Version: 0.3.0'),
              p('Release Date: 2016-03-14'),
              p('Author: Zhenke Wu, Scott Zeger'),
              p('Maintainer: Zhenke Wu <zhwu@jhu.edu>'),
              p(a("Get source code",href="https://github.com/zhenkewu/baker",target="_blank")),
              p(a("Visit my research page",href="http://zhenkewu.com",target="_blank"),
                p(a("Email 'baker' developlers",href="mailto:zhenkewu@gmail.com",target="_blank"))
              )
              
            ),
            tagList(
              tags$p(HTML('Research reported in this work was partially funded 
                                     through a Patient-Centered Outcomes Research Institute (PCORI) Award ME-1408-20318 and a Gates Foundation Grant 48968.')),
              tags$p(img(src="inhealth_logo.jpg",style="display: inline-block; max-width: 25%; width: 25%; height= auto; margin-right: 10px; 
                                    margin-top: 5px; margin-bottom: 20px"),
                     img(src="pcori_logo.jpg",style="display: inline-block; max-width: 25%; width: 25%; height= auto; margin-right: 10px; 
                                    margin-top: 5px; margin-bottom: 20px"),
                     img(src="gates_logo.jpg",style="display: inline-block; max-width: 25%; width: 25%; height= auto;margin-right: 10px; 
                                    margin-top: 5px; margin-bottom: 20px"),
                     img(src="PERCH_logo.jpg",style="display: inline-block; max-width: 25%; width: 25%; height= auto; margin-right: 10px; 
                                    margin-top: 5px; margin-bottom: 20px")
              )
            )
    ),
    
    tabItem(tabName="Messages",
            mywellPanel(
              dataTableOutput("message_df")
              
            )
    )
    
    
    
    
  ) # end of tabItems.
)

dashboardPage(
  header,
  sidebar,
  body
)
