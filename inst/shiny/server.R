################################################################################
# User Interface for Baker (server file)
#
# Designed for PERCH data in the current version. Please do not hesitate to
# contact the author(s) below for software design and modifications.
#
# Zhenke Wu
# zhenkewu@gmail.com
# 1st version: March 17, 2016
################################################################################

# adjust the maximum allowed file size from csv:
options(shiny.maxRequestSize=30*1024^2) 

#library(baker)
library(shinyFiles)
library(shinydashboard)


# Define global values; can modify if other options are required, e.g., more
# causes, more stratifying variables and their values, etc:
sitename_choices <- c("Kenya"="01KEN",
                      "Gambia"="02GAM",
                      "Mali"="03MAL",
                      "Zambia"="04ZAM",
                      "South Africa"="05SAF",
                      "Thailand"="06THA",
                      "Bangladesh"="07BAN"
)
cause_list_choices <- c("BOPE"="BOPE",
                        "C_PNEU"="C_PNEU",
                        "M_PNEU"="M_PNEU",
                        "PCP"="PCP",
                        "ADENO"="ADENO",
                        "CMV"="CMV",
                        "COR"="COR",
                        "FLU_C"="FLU_C",
                        "HBOV"="HBOV",
                        "HMPV_A_B"="HMPV_A_B",
                        "FLU_A"="FLU_A",
                        "FLU_B"="FLU_B",
                        "PARA_1"="PARA_1",
                        "PARA_2"="PARA_2",
                        "PARA_3"="PARA_3",
                        "PARA_4"="PARA_4",
                        "PV_EV"= "PV_EV",
                        "RHINO"="RHINO",
                        "RSV"="RSV",
                        "HINF"="HINF",
                        "MCAT"="MCAT",
                        "PNEU_VT13"="PNEU_VT13",
                        "PNEU_NOVT13"="PNEU_NOVT13",
                        "SASP"="SASP",
                        "SAUR"="SAUR"
)

patho_BrS_NPPCR_choices <- c("BOPE"="BOPE",
                             "C_PNEU"="C_PNEU",
                             "M_PNEU"="M_PNEU",
                             "PCP"="PCP",
                             "ADENO"="ADENO",
                             "CMV"="CMV",
                             "COR"="COR",
                             "FLU_C"="FLU_C",
                             "HBOV"="HBOV",
                             "HMPV_A_B"="HMPV_A_B",
                             "FLU_A"="FLU_A",
                             "FLU_B"="FLU_B",
                             "PARA_1"="PARA_1",
                             "PARA_2"="PARA_2",
                             "PARA_3"="PARA_3",
                             "PARA_4"="PARA_4",
                             "PV_EV"= "PV_EV",
                             "RHINO"="RHINO",
                             "RSV"="RSV",
                             "HINF"="HINF",
                             "MCAT"="MCAT",
                             "PNEU"="PNEU",
                             "SASP"="SASP",
                             "SAUR"="SAUR"
)
patho_BrS_WBPCR_choices         <- c("PNEU"="PNEU")
patho_BrS_NPCX_VT13_choices     <- c("PNEU_VT13"="PNEU_VT13",
                                     "PNEU_NOVT13"="PNEU_NOVT13")

patho_SS_BCX_choices <- c("HINF"="HINF",
                          "MCAT"="MCAT",
                          "PNEU_VT13"="PNEU_VT13",
                          "PNEU_NOVT13"="PNEU_NOVT13",
                          "SASP"="SASP",
                          "SAUR"="SAUR")
patho_SS_LA_ADJ_choices <- patho_SS_BCX_choices
patho_SS_PF_ADJ_choices <- patho_SS_BCX_choices

X_extra_choices <- c("Enrollment Date"= "ENRLDATE",
                     "Patient ID"     = "patid",
                     "Age Category"   = "AGECAT",
                     "HIV Status"     = "HIV")

X_strat <- "newSITE"


# prior choices:
Eti_prior_choices <- c("Symmetric Dirichlet"="sym_dirichlet")



################################################################################

# utility functions:

clean_list <- function(A){
  res <- A[!sapply(A, is.null)]
  if (length(res)==0){
    return(NULL)
  }
  res
}

sidewellPanel <- function(...){
  wellPanel(style="background-color: #222D32;",...)
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

shinyServer(function(input,output,session) {
  
  
  output$baker_statusUI <- renderUI(
    bootstrapPage(
      # Add custom CSS & Javascript;
      tagList(
        tags$head(
          tags$link(rel="stylesheet", type="text/css",href="style.css"),
          tags$script(type="text/javascript", src = "busy.js")
        )
      ),
      div(class = "busy",  
          p("Calculation in progress.."), 
          img(src="ajaxloaderq.gif")
      )#,
      #div(class = "span4", uiOutput("obs")),
      #div(class = "span8", verbatimTextOutput("data_nplcm_structure"))
    )
  )

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
  
  #
  # site name:
  #
  output$siteUI <- renderUI({
    selectInput("sitename","Site (Split Data By Site)",sitename_choices,selected=c("Gambia"="02GAM"))
  })
  
  # covariate strata:
  output$X_extra_UI <- renderUI({
    if(input$selectall_X_extra%%2 == 1) {
      checkboxGroupInput("X_extra","Covariates", X_extra_choices,inline=!TRUE)
    } else{
      checkboxGroupInput("X_extra","Covariates", X_extra_choices,
                         selected = X_extra_choices, inline=FALSE)
    }
  })
  
  
  #
  # Specify causes and measurements:
  #
  
  # causes of disease:
  output$cause_listUI <- renderUI({

    if(input$selectall_cause_list%%2 == 1) {
      checkboxGroupInput("cause_list","Causes",cause_list_choices,inline=!TRUE)
    } else{
      checkboxGroupInput("cause_list","Causes",cause_list_choices,
                         selected=cause_list_choices,inline=!TRUE) 
    }
  })
  
  # bronze-standard data:
  output$BrS_object_NPPCR_UI <- renderUI({
    if(input$selectall_patho_BrS_NPPCR%%2 == 1) {
      checkboxGroupInput("patho_BrS_NPPCR","Nasal-Pharyngeal PCR Measurements", patho_BrS_NPPCR_choices,inline=!TRUE)
    } else{
      checkboxGroupInput("patho_BrS_NPPCR","Nasal-Pharyngeal PCR Measurements",patho_BrS_NPPCR_choices,
                         selected=patho_BrS_NPPCR_choices,inline=!TRUE) 
    }
  })
  
  output$BrS_object_WBPCR_UI <- renderUI({
    if(input$selectall_patho_BrS_WBPCR%%2 == 1) {
      checkboxGroupInput("patho_BrS_WBPCR","Whole Blood PCR Measurements", patho_BrS_WBPCR_choices,inline=!TRUE)
    } else{
      checkboxGroupInput("patho_BrS_WBPCR","Whole Blood PCR Measurements",patho_BrS_WBPCR_choices,
                         selected=patho_BrS_WBPCR_choices,inline=!TRUE) 
    }
  })
  
  output$BrS_object_NPCX_VT13_UI <- renderUI({
    if(input$selectall_patho_BrS_NPCX_VT13%%2 == 1) {
      checkboxGroupInput("patho_BrS_NPCX_VT13","NPCX_VT13 Measurements", patho_BrS_NPCX_VT13_choices,inline=!TRUE)
    } else{
      checkboxGroupInput("patho_BrS_NPCX_VT13","NPCX_VT13 Measurements",patho_BrS_NPCX_VT13_choices,
                         selected=patho_BrS_NPCX_VT13_choices,inline=!TRUE) 
    }
  })
  
  # silver-standard data:
  output$SS_object_BCX_UI <- renderUI({
    if(input$selectall_patho_SS_BCX%%2 == 1) {
      checkboxGroupInput("patho_SS_BCX","Blood Culture Measurements", patho_SS_BCX_choices,inline=!TRUE)
    } else{
      checkboxGroupInput("patho_SS_BCX","Blood Culture Measurements",patho_SS_BCX_choices,
                         selected=patho_SS_BCX_choices,inline=!TRUE) 
    }
  })
  
  output$SS_object_LA_ADJ_UI <- renderUI({
    if(input$selectall_patho_SS_LA_ADJ%%2 == 1) {
      checkboxGroupInput("patho_SS_LA_ADJ","Lung Aspirate Adjudicated Measurements", patho_SS_LA_ADJ_choices,inline=!TRUE)
    } else{
      checkboxGroupInput("patho_SS_LA_ADJ","Lung Aspirate Adjudicated Measurements",patho_SS_LA_ADJ_choices,
                         selected=patho_SS_LA_ADJ_choices,inline=!TRUE) 
    }
  })
  
  output$SS_object_PF_ADJ_UI <- renderUI({
    if( input$selectall_patho_SS_PF_ADJ%%2 == 1) {
      checkboxGroupInput("patho_SS_PF_ADJ","Pleural Fluid Adjudicated Measurements",patho_SS_PF_ADJ_choices,inline=!TRUE)
    } else{
      checkboxGroupInput("patho_SS_PF_ADJ","Pleural Fluid Adjudicated Measurements",patho_SS_PF_ADJ_choices,
                         selected=patho_SS_PF_ADJ_choices,inline=!TRUE) 
    }
  })
  
  output$frontpage_upload_lab_data_warning <- renderPrint({
      if (is.null(input$input_lab_data)){
            cat("==[baker UI] Please upload lab test data!==\n")
      } else{
           cat("==[baker UI] Great! Lab test data loaded!\n
                             Please choose data cleaning options on the sidebar.\n
                             Come back to this tab to specify causes and available measurements.==\n")
        }
  })
  
  
  # clean options: reactive:
  curr_clean_options <- reactive({ 
    if (input$refresh_button >0 ){
      
      if (is.null(input$input_lab_data)){
        cat("==[baker UI] Please upload lab test data!==\n")
      } else{
        
        BrS_objects <- vector("list",3)
        if (!is.null(input$patho_BrS_NPPCR)){
          BrS_objects[[1]] <- make_meas_object(input$patho_BrS_NPPCR,"NP","PCR","BrS",input$cause_list)
        }
        if (!is.null(input$patho_BrS_WBPCR)){
          BrS_objects[[2]] <- make_meas_object(input$patho_BrS_WBPCR,"WB","PCR","BrS",input$cause_list)
        }
        if (!is.null(input$patho_BrS_NPCX_VT13)){
          BrS_objects[[3]] <- make_meas_object(input$patho_BrS_NPCX_VT13,"NP","CX","BrS",input$cause_list)
        }
        
        SS_objects <- vector("list",3)
        if (!is.null(input$patho_SS_BCX)){
          SS_objects[[1]] <- make_meas_object(input$patho_SS_BCX,"B","CX","SS",input$cause_list)
        }
        if (!is.null(input$patho_SS_LA_ADJ)){
          SS_objects[[2]] <- make_meas_object(input$patho_SS_LA_ADJ,"LA","_ADJ","SS",input$cause_list)
        }
        if (!is.null(input$patho_SS_PF_ADJ)){
          SS_objects[[3]] <- make_meas_object(input$patho_SS_PF_ADJ,"PF","_ADJ","SS",input$cause_list)  
        }
        
        BrS_objects_clean <- clean_list(BrS_objects)
        SS_objects_clean  <- clean_list(SS_objects)
        
        ## clean PERCH data:
        list(raw_meas_dir       =  input$input_lab_data$datapath,                 # <-------- where to look for raw data.
             case_def           =  input$casedef,                                # <------- case definition variable.
             case_def_val       =  input$casedef_val,                                                    # <---------- case definition variable's value.
             ctrl_def           =  input$ctrldef,                                           # <------control definition variable.
             ctrl_def_val       =  input$ctrldef_val,                                                    # <------ control definition variable's value.
             X_strat            =  X_strat,                                            # <---- focus on the stratum defined by this variable.
             X_strat_val        =  list(input$sitename),                                       # <---- the stratum definition value.
             BrS_objects        =  BrS_objects_clean,    # <---- all bronze-standard measurements.
             SS_objects         =  SS_objects_clean,
             X_extra            =  input$X_extra,                         # <----- covariates besides case/control status.
             patho_taxo_dir     =  input$input_taxo_data$datapath)        # <---- where to look for pathogen taxonomy information.
        # you can add: "date_formats".
      }
    }
  })
  
  curr_data_nplcm <- reactive({ 
    if (input$refresh_button >0 | input$fit_button > 0){
      clean_perch_data(curr_clean_options()) # <---- time-consuming; let it be reactive conductor. Only update when button is pressed.
    }
  })
  
  output$data_nplcm_structure <- renderPrint({
    if (input$refresh_button >0 && !is.null(input$input_lab_data)){
      cat("==[baker UI] Data Cleaning Summaries==\n")
      str(curr_clean_options())
      cat("==[baker UI] Data Summaries==\n")
      str(curr_data_nplcm())
    } else{
      cat("==[baker UI] Please \n
          1) upload lab test data (sidebar),\n
          2) choose causes of disease and measurements (click 'Choose' in sidebar and click the 2nd tab in the main panel), and \n
          3) press 'Set/Update Parameters' button (bottom of the sidebar). ==\n")
    }
  })
  
  output$logORmat <-renderPlot({
    if (input$plot_pairLOR_button >0){
      isolate({
        # pairwise log odds ratio plot:
        plot_logORmat(curr_data_nplcm(),input$patho_BrS_NPPCR,1)
      })
    } 
  })
  
  # get which are selected:
  subset_vars <- reactiveValues()
  observe({
    if (input$refresh_button >0 && !is.null(input$input_lab_data)){
      subset_vars$ind_for_season <- sapply(input$patho_BrS_NPPCR,function(x){which(patho_BrS_NPPCR_choices==x)})
    }
  })
  
  output$choose_patho_for_season_UI <- renderUI({
    if(input$selectall_choose_patho_for_season == 0 | input$selectall_choose_patho_for_season%%2 == 0) {
      checkboxGroupInput("patho_for_season","Select Pathogens to Visualize Seasonal Trends (NPPCR)", 
                         patho_BrS_NPPCR_choices[ subset_vars$ind_for_season],inline=TRUE)
    } else{
      checkboxGroupInput("patho_for_season","Select Pathogens to Visualize Seasonal Trends (NPPCR)", 
                         patho_BrS_NPPCR_choices[ subset_vars$ind_for_season],
                         selected=patho_BrS_NPPCR_choices[ subset_vars$ind_for_season],inline=TRUE) 
    }
    
  })
  
  output$visualize_season <- renderPlot({
    if (input$refresh_button >0 && input$visualize_season_button > 0 ){
      isolate({
        if (!("ENRLDATE" %in% input$X_extra)){
          cat("==[baker UI] Please select 'Enrollment Date' to visualize seasonal trend!==\n")
        } else if (is.null(input$patho_for_season)){
          cat("==[baker UI] Please choose one or more pathogens to visualize seasonal trend!'==\n")
        }else{
          data_nplcm0 <-  curr_data_nplcm()
          ind_notNA   <- which(rowSums(is.na(data_nplcm0$Mobs$MBS$NPPCR))>0)
          data_nplcm_season <- data_nplcm0
          
          if (length(ind_notNA)>0){
            data_nplcm_season$Mobs$MBS$NPPCR <-  data_nplcm0$Mobs$MBS$NPPCR[-ind_notNA,,drop=FALSE]
            data_nplcm_season$Mobs$MSS$BCX <-  data_nplcm0$Mobs$MSS$BCX[-ind_notNA,,drop=FALSE]
            data_nplcm_season$X <-  data_nplcm0$X[-ind_notNA,]
            data_nplcm_season$Y <-  data_nplcm0$Y[-ind_notNA]
          }
          
          # explore covariate effect: focus on seasonality.
          ind_for_season_UI <- sapply(input$patho_for_season,function(x){which(input$patho_BrS_NPPCR==x)})
          par(mfrow=c(length(ind_for_season_UI),1))
          par(mar=c(1,5,0,5),oma = c(2,2,0,0),xpd=TRUE)
          for (i in ind_for_season_UI){
            par(mar=c(1,5,3,5))
            visualize_season(data_nplcm_season,i,1)
          }
        }
      })
    }
  })
  
  plot_height <- function() {
    return(length(input$patho_for_season)*300)
  }
  
  output$visualize_season_UI <- renderUI({
    if (input$refresh_button > 0 && !is.null(input$input_lab_data) && input$visualize_season_button > 0){
      isolate({
        plotOutput("visualize_season",height = plot_height(), width = "100%")
      })
    }
  })
  
  
  #
  # EDA UI:
  #
  output$eda_UI <- renderUI({
    
    if (input$eda_options == 'data_structure'){
      tabsetPanel(
        tabPanel("Data Structure",
                 verbatimTextOutput("data_nplcm_structure")
        )
      )
    } else if (input$eda_options=='pairLOR'){
      tabsetPanel(
        tabPanel("Pairwise Log Odds Ratio (Cases and Controls)",
                 plotOutput("logORmat",height=plot_height_pairLOR(), width='100%')
        )
      )
    } else {
      tabsetPanel(
        tabPanel("Seasonal Trend",
                 uiOutput("visualize_season_UI")
        )
      )
    }
  })
  
  #
  # visualize UI:
  #
  
  output$visualize_UI <- renderUI({
    if (input$visualize_options == "pie_standardization"){
      tabsetPanel(
        tabPanel("Pie Standardizations",
                 uiOutput("plot_pie_standardization_UI")
        )
      )
    } else if (input$visualize_options == "individual_prediction"){
      tabsetPanel(
        tabPanel("Individual Prediction",
                 uiOutput("plot_individual_prediction_UI")
        )
      )
    } else { #(input$visualize_options == "group_etiology"){
      tabsetPanel(
        tabPanel("Bacterial vs Viral",
                 uiOutput("plot_group_etiology_UI")
        )
      )  
    }
  })
  
  #
  # Specify priors:
  #
  
  # etiology priors:
  output$Eti_prior_UI <- renderUI({
    tagList(
      selectInput("Eti_prior_name","Etiology Prior Distribution",choices=Eti_prior_choices, selected="sym_dirichlet"),
      conditionalPanel(condition = "input.Eti_prior_name=='sym_dirichlet'",
                       numericInput("sym_dirichlet_alpha","Symmetric Dirichlet Parameter (alpha > 0)",1,min=0.001,max=100, step=1))
    )
  })
  
  # TPR priors for bronze-standard data:
  TPR_prior_BrS_info_choices <- c("non-informative"="non-informative",
                                  "informative"="informative")
  TPR_prior_BrS_input_choices <- c("match prior range"="match_range",
                                   "direct Beta parameters"="direct_beta_param")
  
  output$status_patho_BrS_NPPCR <- reactive({
    return(!is.null(input$patho_BrS_NPPCR))
  })
  outputOptions(output, "status_patho_BrS_NPPCR", suspendWhenHidden = FALSE)
  
  output$status_patho_BrS_WBPCR <- reactive({
    return(!is.null(input$patho_BrS_WBPCR))
  })
  outputOptions(output, "status_patho_BrS_WBPCR", suspendWhenHidden = FALSE)
  
  output$status_patho_BrS_NPCX_VT13 <- reactive({
    return(!is.null(input$patho_BrS_NPCX_VT13))
  })
  outputOptions(output, "status_patho_BrS_NPCX_VT13", suspendWhenHidden = FALSE)
  
  output$TPR_prior_BrS_info_sidebarUI <- renderUI({
    if (input$refresh_button > 0 && !is.null(input$input_lab_data)){
      tagList(
        conditionalPanel(condition = "(output.status_patho_BrS_NPPCR || output.status_patho_BrS_WBPCR || output.status_patho_BrS_NPCX_VT13)",
                         sidewellPanel(
                           h5("Bronze-Standard Measurements"),
                           selectInput("TPR_prior_BrS_info", "True Positive Rate Prior",choices = TPR_prior_BrS_info_choices,selected="informative"),
                           conditionalPanel(condition = "input.TPR_prior_BrS_info=='informative'",
                                            selectInput("TPR_prior_BrS_input", "How to specify informative TPR prior?",choices = TPR_prior_BrS_input_choices,selected="match_range")
                           )
                         )
        )
      )
    }
  })
  
  output$TPR_prior_BrS_UI <- renderUI({
    if (input$refresh_button > 0 && !is.null(input$input_lab_data)){
    tagList(
      conditionalPanel(condition = "(output.status_patho_BrS_NPPCR || output.status_patho_BrS_WBPCR || output.status_patho_BrS_NPCX_VT13)",
                       box(width='100%',
                         title=c("Bronze-Standard Measurements"),status="primary",collapsible = TRUE,solidHeader = TRUE,
                         #selectInput("TPR_prior_BrS_info", "True Positive Rate Prior",choices = TPR_prior_BrS_info_choices,selected="informative"),
                         conditionalPanel(condition = "input.TPR_prior_BrS_info=='informative'",
                                          #selectInput("TPR_prior_BrS_input", "How to specify informative TPR prior?",choices = TPR_prior_BrS_input_choices,selected="match_range"),
                                          conditionalPanel(condition="input.TPR_prior_BrS_input=='match_range'",
                                                           h4("Specify Upper and Lower Ranges:"),
                                                           conditionalPanel(condition = "output.status_patho_BrS_NPPCR",
                                                                            wellPanel(
                                                                              h5("NPPCR"),
                                                                              numericInput("val_NPPCR_down","Lower Range",0.5,min=0.001,max=0.995,step=0.01),
                                                                              numericInput("val_NPPCR_up","Upper Range",0.99,min=0.001,max=0.995,step=0.01),
                                                                              renderPlot({
                                                                                beta_parms_from_quantiles(c(input$val_NPPCR_down,input$val_NPPCR_up),plot=TRUE)
                                                                              }))),
                                                           conditionalPanel(condition = "output.status_patho_BrS_WBPCR",
                                                                            wellPanel(
                                                                              h5("WBPCR"),
                                                                              numericInput("val_WBPCR_down","Lower Range",0.5,min=0.001,max=0.995,step=0.01),
                                                                              numericInput("val_WBPCR_up","Upper Range",0.99,min=0.001,max=0.995,step=0.01),
                                                                              renderPlot({
                                                                                beta_parms_from_quantiles(c(input$val_WBPCR_down,input$val_WBPCR_up),plot=TRUE)
                                                                              }))),
                                                           conditionalPanel(condition = "output.status_patho_BrS_NPCX_VT13",
                                                                            wellPanel(
                                                                              h5("NPCX"),
                                                                              numericInput("val_NPCX_down","Lower Range",0.5,min=0.001,max=0.995,step=0.01),
                                                                              numericInput("val_NPCX_up","Upper Range",0.99,min=0.001,max=0.995,step=0.01),
                                                                              renderPlot({
                                                                                beta_parms_from_quantiles(c(input$val_NPCX_down,input$val_NPCX_up),plot=TRUE)
                                                                              })
                                                                            ))),
                                          
                                          conditionalPanel(condition="input.TPR_prior_BrS_input=='direct_beta_param'",
                                                           h4("Specify Beta Parameters:"),
                                                           conditionalPanel(condition = "output.status_patho_BrS_NPPCR",
                                                                            wellPanel(
                                                                              h5("NPPCR"),
                                                                              numericInput("val_NPPCR_alpha","alpha",6,min=0.001,max=500,step=1),
                                                                              numericInput("val_NPPCR_beta","beta",2,min=0.001,max=500,step=1),
                                                                              renderPlot({
                                                                                beta_parms_from_quantiles(qbeta(c(0.025,0.975),input$val_NPPCR_alpha,input$val_NPPCR_beta),plot=TRUE)
                                                                              }))),
                                                           conditionalPanel(condition = "output.status_patho_BrS_WBPCR",
                                                                            wellPanel(
                                                                              h5("WBPCR"),
                                                                              numericInput("val_WBPCR_alpha","alpha",6,min=0.001,max=500,step=1),
                                                                              numericInput("val_WBPCR_beta","beta",2,min=0.001,max=500,step=1),
                                                                              renderPlot({
                                                                                beta_parms_from_quantiles(qbeta(c(0.025,0.975),input$val_WBPCR_alpha,input$val_WBPCR_beta),plot=TRUE)
                                                                              }))),
                                                           conditionalPanel(condition = "output.status_patho_BrS_NPCX_VT13",
                                                                            wellPanel(
                                                                              h5("NPCX"),
                                                                              numericInput("val_NPCX_alpha","alpha",6,min=0.001,max=500,step=1),
                                                                              numericInput("val_NPCX_beta","beta",2,min=0.001,max=500,step=1),
                                                                              renderPlot({
                                                                                beta_parms_from_quantiles(qbeta(c(0.025,0.975),input$val_NPCX_alpha,input$val_NPCX_beta),plot=TRUE)
                                                                              })))
                                          )
                         ),
                         # non-informative:
                         conditionalPanel(condition = "input.TPR_prior_BrS_info=='non-informative'",
                                          h4("Specify Beta Parameters:"),
                                          conditionalPanel(condition = "output.status_patho_BrS_NPPCR",
                                                           wellPanel(
                                                             h5("NPPCR"),
                                                             numericInput("val_NPPCR_alpha_flat","alpha",1,min=1,max=1,step=1),
                                                             numericInput("val_NPPCR_beta_flat","beta",1,min=1,max=1,step=1),
                                                             renderPlot({
                                                               beta_parms_from_quantiles(qbeta(c(0.025,0.975),input$val_NPPCR_alpha_flat,input$val_NPPCR_beta_flat),plot=TRUE)
                                                             }))),
                                          conditionalPanel(condition = "output.status_patho_BrS_WBPCR",
                                                           wellPanel(
                                                             h5("WBPCR"),
                                                             numericInput("val_WBPCR_alpha_flat","alpha",1,min=1,max=1,step=1),
                                                             numericInput("val_WBPCR_beta_flat","beta",1,min=1,max=1,step=1),
                                                             renderPlot({
                                                               beta_parms_from_quantiles(qbeta(c(0.025,0.975),input$val_WBPCR_alpha_flat,input$val_WBPCR_beta_flat),plot=TRUE)
                                                             }))),
                                          conditionalPanel(condition = "output.status_patho_BrS_NPCX_VT13",
                                                           wellPanel(
                                                             h5("NPCX"),
                                                             numericInput("val_NPCX_alpha_flat","alpha",1,min=1,max=1,step=1),
                                                             numericInput("val_NPCX_beta_flat","beta",1,min=1,max=1,step=1),
                                                             renderPlot({
                                                               beta_parms_from_quantiles(qbeta(c(0.025,0.975),input$val_NPCX_alpha_flat,input$val_NPCX_beta_flat),plot=TRUE)
                                                             })))
                         )
                       )
      )
    )
    }
  })
  
  #  TPR prior for silver standard data: 
  TPR_prior_SS_info_choices <- c("informative"="informative")
  TPR_prior_SS_input_choices <- c("match prior range"="match_range",
                                  "direct Beta parameters"="direct_beta_param")
  
  output$status_patho_SS_BCX <- reactive({
    return(!is.null(input$patho_SS_BCX))
  })
  outputOptions(output, "status_patho_SS_BCX", suspendWhenHidden = FALSE)
  
  output$status_patho_SS_LA_ADJ <- reactive({
    return(!is.null(input$patho_SS_LA_ADJ))
  })
  outputOptions(output, "status_patho_SS_LA_ADJ", suspendWhenHidden = FALSE)
  
  output$status_patho_SS_PF_ADJ <- reactive({
    return(!is.null(input$patho_SS_PF_ADJ))
  })
  outputOptions(output, "status_patho_SS_PF_ADJ", suspendWhenHidden = FALSE)
  
  
  output$TPR_prior_SS_info_sidebarUI <- renderUI({
    if (input$refresh_button > 0 && !is.null(input$input_lab_data)){
    tagList(
    conditionalPanel(condition = "(output.status_patho_SS_BCX || output.status_patho_SS_LA_ADJ || output.status_patho_SS_PF_ADJ)",
                     sidewellPanel(
                       h5("Silver-Standard Measurements"),
                       selectInput("TPR_prior_SS_info", "True Positive Rate Prior",choices = TPR_prior_SS_info_choices,selected="informative"),
                       conditionalPanel(condition = "input.TPR_prior_SS_info=='informative'",
                                        selectInput("TPR_prior_SS_input", "How to specify informative TPR prior?",choices = TPR_prior_SS_input_choices,selected="match_range")
                       )
                     )
    )
    )
    }
  })
  
  output$TPR_prior_SS_UI <- renderUI({
    if (input$refresh_button > 0 && !is.null(input$input_lab_data)){
    tagList(
      conditionalPanel(condition = "(output.status_patho_SS_BCX || output.status_patho_SS_LA_ADJ || output.status_patho_SS_PF_ADJ)",
                       box(width='100%',
                           title=c("Silver-Standard Measurements"),status="warning",collapsible = TRUE,solidHeader = TRUE,
                         #selectInput("TPR_prior_SS_info", "True Positive Rate Prior",choices = TPR_prior_SS_info_choices,selected="informative"),
                         conditionalPanel(condition = "input.TPR_prior_SS_info=='informative'",
                                          #selectInput("TPR_prior_SS_input", "How to specify informative TPR prior?",choices = TPR_prior_SS_input_choices,selected="match_range"),
                                          conditionalPanel(condition="input.TPR_prior_SS_input=='match_range'",
                                                           h4("Specify Upper and Lower Ranges:"),
                                                           conditionalPanel(condition = "output.status_patho_SS_BCX",
                                                                            wellPanel(
                                                                              h5("BCX"),
                                                                              numericInput("val_BCX_down","Lower Range",0.05,min=0.001,max=0.995,step=0.01),
                                                                              numericInput("val_BCX_up","Upper Range",0.15,min=0.001,max=0.995,step=0.01),
                                                                              renderPlot({
                                                                                beta_parms_from_quantiles(c(input$val_BCX_down,input$val_BCX_up),plot=TRUE)
                                                                              }))),
                                                           conditionalPanel(condition = "output.status_patho_SS_LA_ADJ",
                                                                            wellPanel(
                                                                              h5("LA_ADJ"),
                                                                              numericInput("val_LA_ADJ_down","Lower Range",0.05,min=0.001,max=0.995,step=0.01),
                                                                              numericInput("val_LA_ADJ_up","Upper Range",0.15,min=0.001,max=0.995,step=0.01),
                                                                              renderPlot({
                                                                                beta_parms_from_quantiles(c(input$val_LA_ADJ_down,input$val_LA_ADJ_up),plot=TRUE)
                                                                              }))),
                                                           conditionalPanel(condition = "output.status_patho_SS_PF_ADJ",
                                                                            wellPanel(
                                                                              h5("PF_ADJ"),
                                                                              numericInput("val_PF_ADJ_down","Lower Range",0.05,min=0.001,max=0.995,step=0.01),
                                                                              numericInput("val_PF_ADJ_up","Upper Range",0.15,min=0.001,max=0.995,step=0.01),
                                                                              renderPlot({
                                                                                beta_parms_from_quantiles(c(input$val_PF_ADJ_down,input$val_PF_ADJ_up),plot=TRUE)
                                                                              })))
                                          ),
                                          conditionalPanel(condition="input.TPR_prior_SS_input=='direct_beta_param'",
                                                           h4("Specify Beta Parameters:"),
                                                           conditionalPanel(condition = "output.status_patho_SS_BCX",
                                                                            wellPanel(
                                                                              h5("BCX"),
                                                                              numericInput("val_BCX_alpha","alpha",1,min=0.001,max=500,step=1),
                                                                              numericInput("val_BCX_beta","beta",9,min=0.001,max=500,step=1),
                                                                              renderPlot({
                                                                                beta_parms_from_quantiles(qbeta(c(0.025,0.975),input$val_BCX_alpha,input$val_BCX_beta),plot=TRUE)
                                                                              }))),
                                                           conditionalPanel(condition = "output.status_patho_SS_LA_ADJ",
                                                                            wellPanel(
                                                                              h5("LA_ADJ"),
                                                                              numericInput("val_LA_ADJ_alpha","alpha",1,min=0.001,max=500,step=1),
                                                                              numericInput("val_LA_ADJ_beta","beta",9,min=0.001,max=500,step=1),
                                                                              renderPlot({
                                                                                beta_parms_from_quantiles(qbeta(c(0.025,0.975),input$val_LA_ADJ_alpha,input$val_LA_ADJ_beta),plot=TRUE)
                                                                              }))),
                                                           conditionalPanel(condition = "output.status_patho_SS_PF_ADJ",
                                                                            wellPanel(
                                                                              h5("PF_ADJ"),
                                                                              numericInput("val_PF_ADJ_alpha","alpha",1,min=0.001,max=500,step=1),
                                                                              numericInput("val_PF_ADJ_beta","beta",9,min=0.001,max=500,step=1),
                                                                              renderPlot({
                                                                                beta_parms_from_quantiles(qbeta(c(0.025,0.975),input$val_PF_ADJ_alpha,input$val_PF_ADJ_beta),plot=TRUE)
                                                                              })))
                                          )
                         )
                       )
      )
    )
    }
  })
  
  #
  # specify likelihood:
  #
  output$Eti_formula_UI <- renderUI({
    textInput("Eti_formula_input","",value="~ 0")
  })
  
  output$FPR_formula_UI <- renderUI({
    tagList(
      conditionalPanel(condition = "output.status_patho_BrS_NPPCR",
                       textInput("FPR_formula_input_NPPCR","NPPCR",value="~ 0")
      ),
      conditionalPanel(condition = "output.status_patho_BrS_WBPCR",
                       textInput("FPR_formula_input_WBPCR","WBPCR",value="~ 0")
      ),
      conditionalPanel(condition = "output.status_patho_BrS_NPCX_VT13",
                       textInput("FPR_formula_input_NPCX","NPCX",value="~ 0")
      )
    )
  })
  
  # fit the model:
  curr_model_options <- reactive({
    if (input$refresh_button > 0 ){
      #isolate({
      
      #
      # bronze-standard:
      #
      BrS_vals <- list()
      if (input$TPR_prior_BrS_info=="informative"){
        if ( input$TPR_prior_BrS_input=="match_range"){
          if (!is.null(input$patho_BrS_NPPCR)){
            BrS_vals[["NPPCR"]] <- list(up  = rep(input$val_NPPCR_up,length(input$patho_BrS_NPPCR)),
                                        low = rep(input$val_NPPCR_down,length(input$patho_BrS_NPPCR)))
          }
          
          if (!is.null(input$patho_BrS_WBPCR)){
            BrS_vals[["WBPCR"]] <- list(up  = rep(input$val_WBPCR_up,length(input$patho_BrS_WBPCR)),
                                        low = rep(input$val_WBPCR_down,length(input$patho_BrS_WBPCR)))
          }
          
          if (!is.null(input$patho_BrS_NPCX_VT13)){
            BrS_vals[["NPCX"]] <- list(up  = rep(input$val_NPCX_up,length(input$patho_BrS_NPCX_VT13)),
                                       low = rep(input$val_NPCX_down,length(input$patho_BrS_NPCX_VT13)))
          }
        } else {
          if (!is.null(input$patho_BrS_NPPCR)){
            BrS_vals[["NPPCR"]] <- list(alpha  = rep(input$val_NPPCR_alpha,length(input$patho_BrS_NPPCR)),
                                        beta = rep(input$val_NPPCR_beta,length(input$patho_BrS_NPPCR)))
          }
          
          if (!is.null(input$patho_BrS_WBPCR)){
            BrS_vals[["WBPCR"]] <- list(alpha  = rep(input$val_WBPCR_alpha,length(input$patho_BrS_WBPCR)),
                                        beta = rep(input$val_WBPCR_beta,length(input$patho_BrS_WBPCR)))
          }
          
          if (!is.null(input$patho_BrS_NPCX_VT13)){
            BrS_vals[["NPCX"]] <- list(alpha  = rep(input$val_NPCX_alpha,length(input$patho_BrS_NPCX_VT13)),
                                       beta = rep(input$val_NPCX_beta,length(input$patho_BrS_NPCX_VT13)))
          }
        }
      } else{
        BrS_vals <- list()
        if (!is.null(input$patho_BrS_NPPCR)){
          BrS_vals[["NPPCR"]] <- list(alpha  = rep(input$val_NPPCR_alpha_flat,length(input$patho_BrS_NPPCR)),
                                      beta   = rep(input$val_NPPCR_beta_flat,length(input$patho_BrS_NPPCR)))
        }
        
        if (!is.null(input$patho_BrS_WBPCR)){
          BrS_vals[["WBPCR"]] <- list(alpha  = rep(input$val_WBPCR_alpha_flat,length(input$patho_BrS_WBPCR)),
                                      beta = rep(input$val_WBPCR_beta_flat,length(input$patho_BrS_WBPCR)))
        }
        
        if (!is.null(input$patho_BrS_NPCX_VT13)){
          BrS_vals[["NPCX"]] <- list(alpha  = rep(input$val_NPCX_alpha_flat,length(input$patho_BrS_NPCX_VT13)),
                                     beta = rep(input$val_NPCX_beta_flat,length(input$patho_BrS_NPCX_VT13)))
        }
      }
      
      #
      # silver-standard:
      #
      SS_vals <- list()
      if (input$TPR_prior_SS_info=="informative") {
        if (input$TPR_prior_SS_input=="match_range") {
          if (!is.null(input$patho_SS_BCX)){
            SS_vals[["BCX"]] <- list(list(up = rep(input$val_BCX_up,length(input$patho_SS_BCX)),
                                          low = rep(input$val_BCX_down,length(input$patho_SS_BCX))))
          }
          
          if (!is.null(input$patho_SS_LA_ADJ)){
            SS_vals[["LA_ADJ"]] <- list(list(up = rep(input$val_LA_ADJ_up,length(input$patho_SS_LA_ADJ)),
                                             low = rep(input$val_LA_ADJ_down,length(input$patho_SS_LA_ADJ))))
          }
          
          if (!is.null(input$patho_SS_PF_ADJ)){
            SS_vals[["PF_ADJ"]] <- list(list(up  = rep(input$val_PF_ADJ_up,length(input$patho_SS_PF_ADJ)),
                                             low = rep(input$val_PF_ADJ_down,length(input$patho_SS_PF_ADJ))))
          }
        } else {
          if (!is.null(input$patho_SS_BCX)){
            SS_vals[["BCX"]] <- list(list(alpha = rep(input$val_BCX_alpha,length(input$patho_SS_BCX)),
                                          beta = rep(input$val_BCX_beta,length(input$patho_SS_BCX))))
          }
          if (!is.null(input$patho_SS_LA_ADJ)){
            SS_vals[["LA_ADJ"]] <- list(list(alpha = rep(input$val_LA_ADJ_alpha,length(input$patho_SS_LA_ADJ)),
                                             beta = rep(input$val_LA_ADJ_beta,length(input$patho_SS_LA_ADJ))))
          }
          if (!is.null(input$patho_SS_PF_ADJ)){
            SS_vals[["PF_ADJ"]] <- list(list(alpha  = rep(input$val_PF_ADJ_alpha,length(input$patho_SS_PF_ADJ)),
                                             beta = rep(input$val_PF_ADJ_beta,length(input$patho_SS_PF_ADJ))))
          }
        }
      }
      
      list(
        likelihood   = list(
          cause_list = input$cause_list,               # <---- fitted causes.
          k_subclass = c(input$k_subclass_BrS_NPPCR,input$k_subclass_BrS_WBPCR,input$k_subclass_BrS_NPCX_VT13),                         # <---- no. of subclasses.
          Eti_formula = input$Eti_formula_input,                    # <---- only apply FPR formula to specified slice of measurements; if not default to the first slice.
          FPR_formula = list(NPPCR = input$FPR_formula_input_NPPCR,
                             WBPCR = input$FPR_formula_input_WBPCR,
                             NPCX  = input$FPR_formula_input_NPCX)),                     # <---- etiology regression formula.
        use_measurements = input$use_measurements,                           # <---- which measurements to use to inform etiology
        prior        = list(Eti_prior   = overall_uniform(input$sym_dirichlet_alpha, input$cause_list) ,                       # <--- etiology prior. 
                            TPR_prior   = list(
                              BrS  = list(info  = input$TPR_prior_BrS_info,
                                          input = input$TPR_prior_BrS_input,
                                          val   = BrS_vals
                              ),
                              SS   = list(info   = input$TPR_prior_SS_info,
                                          input  = input$TPR_prior_SS_input,
                                          val    = SS_vals
                              )
                            )
        )
      )
      
      #})
    }
    
  })
  
  
  curr_use_meas <- reactiveValues() 
  
  observe({
    use_meas_choices <- list()
    if (input$refresh_button >0 && !is.null(input$input_lab_data) ){
      clean_options <- curr_clean_options()
      if (length(clean_options$BrS_objects) > 0){
        use_meas_choices[["Bronze"]] <- "BrS"
      }
      if (length(clean_options$SS_objects)   > 0){
        use_meas_choices[["Silver"]] <- "SS"
      }
      if (length(use_meas_choices)==0){
        cat("==[baker UI] Please select at least one slice of measurement: 'Input Data (sidebar) ==> Specify Causes and Measurements (tab)'  ! ==\n")
      }
    }
    curr_use_meas$use_meas_choices <- unlist(use_meas_choices)
  })
  
  
  # which subset of measurement qualities to use:
  output$use_measurement_UI <- renderUI({
    checkboxGroupInput("use_measurements","Use Measurements", choices = curr_use_meas$use_meas_choices, 
                       selected=curr_use_meas$use_meas_choices)
  })
  
  observe({
    use_meas <- curr_use_meas$use_meas_choices
    updateCheckboxGroupInput(session,"use_measurements", choices = use_meas, selected=use_meas)
  })
  
  # number of subclasses:
  output$subclass_UI <- renderUI({
    sidewellPanel(
      h5("No. of Subclasses (Bronze Standard)"),
      conditionalPanel(condition = "output.status_patho_BrS_NPPCR",
                       sliderInput("k_subclass_BrS_NPPCR","NPPCR", min=1,max=20,value=1,round=TRUE,step=1)
      ),
      conditionalPanel(condition = "output.status_patho_BrS_WBPCR",
                       sliderInput("k_subclass_BrS_WBPCR","WBPCR", min=1,max=1,value=1,round=TRUE,step=1)
      ),
      conditionalPanel(condition = "output.status_patho_BrS_NPCX_VT13",
                       sliderInput("k_subclass_BrS_NPCX_VT13","NPCX_VT13", min=1,max=1,value=1,round=TRUE,step=1)
      )
    )
  })
  
  
  volumes <- getVolumes() #c('R Installation'=R.home())
  shinyDirChoose(input, 'res_directory', roots=volumes, session=session, restrictions=system.file(package='base'))
  output$display_res_directory_path <- renderPrint({
    if (is.null(input$res_directory)){
      return("==[baker UI] Please choose a directory to store model outputs (posterior samples and figures). ==")
    }
    parseDirPath(volumes, input$res_directory)
  })
  
  #output$savefile <- renderPrint({parseSavePath(volumes, input$save)})
  
  output$specify_mcmc_UI <- renderUI({
    flowLayout(
      sliderInput("n_itermcmc","No. of Iterations Per Chain",min=100,max=1000000,value=100,step=1000,round=TRUE),
      sliderInput("n_burnin","No. of Burn-ins Per Chain",min=1,max=50,value=1,step=1000,round=TRUE),
      sliderInput("n_thin","Thinning Per Chain",min=1,max=100,value =1,step=1,round=TRUE),
      selectInput("individual_pred","Individual Prediction?",choices=c("Yes","No"),selected="No"),
      selectInput("ppd","Calculate and record posterior predictive samples?",choices=c("Yes","No"),selected="No"),
      sliderInput("n_chains","No. of Chains",min=1,max=10,value=1,step=1,round=TRUE)
    )
  })
  
  observe({
    updateSliderInput(session,"n_burnin",max=floor(input$n_itermcmc/5))
  })
  
  
  # set up fitting options:
  fit_option_vars <- reactiveValues()
  result_vars     <- reactiveValues()
  
  observe({
    if (input$refresh_button > 0 && input$fit_button > 0){
      #isolate({
      clean_options <- curr_clean_options()
      
      # date stamp for analysis:
      Date          <- gsub("-", "_", Sys.Date())
      
      if (is.null(input$res_directory)){
        cat("==[baker UI] No working directory specified. Using default...==\n")
        fit_option_vars$working_dir <- paste0(input$default_res_directory,"/")
      } else {
        fit_option_vars$working_dir <- paste0(parseDirPath(volumes, input$res_directory),"/")
      }
      
      # include stratification information in file name:
      fit_option_vars$fullname <- paste0(fit_option_vars$working_dir, Date,"_",
                                                    paste(unlist(clean_options$X_strat_val),sep="_"),"_output")
      
      # for finer scenarios, e.g., different types of analysis applicable to the
      # same data set. Here we just perform one analysis:
      fit_option_vars$result_folder <- fit_option_vars$fullname
      
      # options for MCMC chains:
      fit_option_vars$mcmc_options <- list(debugstatus = TRUE,
                                           n.chains   = input$n_chains,
                                           n.itermcmc = input$n_itermcmc,
                                           n.burnin   = input$n_burnin,
                                           n.thin     = input$n_thin,
                                           individual.pred = input$individual_pred=="Yes",
                                           ppd             = input$ppd=="Yes",
                                           result.folder = fit_option_vars$result_folder,
                                           bugsmodel.dir = fit_option_vars$result_folder,#"C:\\winbugs_model_package\\",
                                           jags.dir = "",
                                           use_jags = TRUE)
      
      # })
    }
  })
  
  observe({
    if (input$fit_button > 0){
      isolate({
        if (is.null(input$input_lab_data)){
          cat("==[baker UI] Please upload lab test data before press Fit Model'. 
                Go to 'Input Data' Menu Item on the sidebar.==\n")
        } else{
          
          data_nplcm    <- curr_data_nplcm()
          model_options <- curr_model_options()
          
          # create folder
          check_dir_create(fit_option_vars$fullname)
          
          # for finer scenarios, e.g., different types of analysis applicable to the
          # same data set. Here we just perform one analysis:
          check_dir_create(fit_option_vars$result_folder)
          
          # <----because the UI allows one to fit many models; so we have to remove the contents (may ask user to save them first):
          #file.remove(file.path(fit_option_vars$mcmc_options$result.folder, list.files(fit_option_vars$mcmc_options$result.folder))) 
          
          # Record the settings of current analysis:
          #data clean options:
          dput(curr_clean_options(),file.path(fit_option_vars$mcmc_options$result.folder,"data_clean_options.txt"))
          result_vars$gs <- nplcm(data_nplcm,model_options,fit_option_vars$mcmc_options)
        }
      })
    }
  })

  output$print_mcmc_options <- renderPrint({
    if (input$refresh_button > 0){
      if (is.null(input$input_lab_data)){
       cat( "==[baker UI] Please upload lab test data before press Fit Model'. 
                1) Go to 'Input Data' Menu Item on the sidebar.
                2) Upload lab data (.csv) and press the 'Initialize' button==\n")
        } else{
          fit_option_vars$mcmc_options
      }
    } else{
      cat( "==[baker UI] Please go back to 'Input Data' on the sidebar. 
                         Come back when data are initialized.==\n")
      }
  })
  
  curr_plot_panels <- function(){
    par(mar=c(2,5,3,5))
    plot_panels(fit_option_vars$result_folder,bg_color=NULL)
  }
  
  output$panels <- renderPlot({
    if (input$fit_button > 0 ){
      isolate({
        if (is.null(input$input_lab_data)){
          plot(0,0,xlim=c(-10,10),ylim=c(-10,10),
               pch="",axes=F, xlab="", ylab="")
          text(0,9,"==[baker UI] Please upload lab test data before press Fit Model'. 
                Go to 'Input Data' Menu Item on the sidebar.==\n",cex=2)
        } else{
        curr_plot_panels()
        }
      })
    }
  })
  
  output$panels_UI <- renderUI({
    if (input$fit_button > 0){
      isolate({
        plotOutput("panels",height = plot_height_panels(), width = plot_width_panels())
      })
    }
  })
  
  output$save_panels_plot <- downloadHandler(
    filename <- function(){paste0(input$panels_plot_name,".pdf")},
    content <- function(file){
      pdf(file,width=as.numeric(input$panels_plot_width),height=as.numeric(input$panels_plot_height))
      curr_plot_panels()
      dev.off()
    }
  )
  
  textInputRow<-function (inputId, label, value = "") 
  {
    div(style="display:inline-block",
        tags$label(label, `for` = inputId), 
        tags$input(id = inputId, type = "text", value = value,class="input-small"))
  }
  
  output$download_panels_UI <- renderUI({
    tagList(
      wellPanel(
        textInputRow("panels_plot_name","File Name","panels_plot"),
        textInputRow("panels_plot_width","Plot Width (inches)",(1+length(curr_clean_options()$BrS_objects)+length(curr_clean_options()$SS_objects))*5),
        textInputRow("panels_plot_height","Plot Height (inches)",length(input$cause_list)*1.5),
        tags$p(),
        tags$p(downloadButton("save_panels_plot","Save Plot"))
      )
    )
  })
  
  output$print_model_options <- renderPrint({
    if (input$refresh_button >0 && !is.null(input$input_lab_data)){
      str(curr_model_options())
    } else{
      cat("==[baker UI] Please \n
          1) upload lab test data (sidebar),\n
          2) choose causes of disease and measurements (click 'Choose' in sidebar and click the 2nd tab in the main panel), and \n
          3) press 'Set/Update Parameters' button (bottom of the sidebar). ==\n")
    }
  })
  
  curr_plot_group_etiology <- function(){
    plot_group_etiology(fit_option_vars$result_folder, curr_clean_options()$patho_taxo_dir,1, 5)
  }
  
  output$plot_group_etiology <- renderPlot({
    if (input$refresh_button >0 && input$fit_button >0 && 
        !is.null(input$input_lab_data) && !is.null(input$input_taxo_data)){
      isolate({
        curr_plot_group_etiology()
      })
    }
  })
  
  output$plot_group_etiology_UI <- renderUI({
    if (input$plot_group_etiology_button > 0){
      isolate({
        plotOutput("plot_group_etiology",width = 1000,height =1000)
      })
    }
  })
  
  
  #############################################################################
  #############################################################################
  #############################################################################
  
  # utility functions within server:
  plot_height_pairLOR <- function() {
    return(max(length(input$patho_BrS_NPPCR)*30,1000))
  }
  
    plot_height_panels <- function() {
    return(length(input$patho_BrS_NPPCR)*80)
  }
  
  plot_width_panels <- function() {
    return(max(1250,
               (1+length(curr_clean_options()$BrS_objects)+length(curr_clean_options()$SS_objects))*250)
    )
  }
  

  #
  # administrative stuff:
  #
  messageData <- reactive({
    if (is.null(input$refresh_button) | input$refresh_button==0){
      mat <<- matrix(c(from="Zhenke Wu", message= "Welcome to baker!",
                       icon="user",time=as.character(Sys.time())),ncol=4)
    }
    
    if (input$fit_button >0){
      isolate({
        if (input$refresh_button >0 ){
          curr_n <- nrow(mat)
          curr_msg_seg <- ""
          if (curr_n==1){
            curr_msg_seg <- "-st model fit can be found in directory:  "
          } else if (curr_n==2){
            curr_msg_seg <- "-nd model fit can be found in directory:  "
          } else if (curr_n==3){
            curr_msg_seg <- "-rd model fit can be found in directory:  "
          } else{
            curr_msg_seg <- "-th model fit can be found in directory:  "
          }
          mat <<- rbind(mat, c("baker UI",
                               paste0(curr_n,curr_msg_seg, 
                                      fit_option_vars$result_folder),
                               "desktop",
                               as.character(Sys.time())))
        }
      })
    }
    
    res <- data.frame(mat)
    colnames(res) <- c("from","message","icon","time")
    return(res)
  })
  
  output$messageMenu <- renderMenu({
    # Code to generate each of the messageItems here, in a list. This assumes
    # that messageData is a data frame with two columns, 'from' and 'message'.
    msgs <- apply(messageData(), 1, function(row) {
      messageItem(from = row[["from"]], message = row[["message"]],
                  icon = icon(row[["icon"]]), 
                  time = row[["time"]])
    })
    
    # This is equivalent to calling:
    #   dropdownMenu(type="messages", msgs[[1]], msgs[[2]], ...)
    dropdownMenu(type = "messages", .list = msgs)
  })
  
  output$message_df <- renderDataTable({
      messageData()[c("time","from","message")]
    
    })
})
