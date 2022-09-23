
# load libraries
library(shiny)
library(shinymanager)
library(shinycssloaders)
library(shinyjs)
library(shinydashboard)
library(shinyWidgets)
library(shinyBS)
source("utils.R")

credentials <- define_credentials()

ui <- fluidPage(
  
  setBackgroundColor(color = "white"),
  tags$head(includeScript("www/returnClick.js")),
  HTML("<br>"),
  
  # sidebar layout
  sidebarLayout(
    sidebarPanel(
      div(img(src="cornell-logo.png",width="90%"),style="text-align: center;"),
      HTML("<br>"),
      
      selectInput("geometry", "Device geometry:",
                  c("Cylinder" = "Cylinder",
                    "Planar slab" = "Slab1",
                    "Concentric cylinder (fixed cell layer thickness of 500 μm)" = "Annulus"),
                  selected = "Slab1"),
      bsTooltip(id="geometry", 
                title="Select islet delivery device geometry", 
                placement = "bottom", 
                trigger = "hover",
                options = NULL),
      
      numericInput("rho","Islet density, ρ (%):", 6.5, step=0.1, min = 1, max = 100),
      bsTooltip(id="rho", 
                title="Enter volumetric islet density (as a percentage) in cell-containing domain", 
                placement = "bottom", 
                trigger = "hover",
                options = NULL),
      
      numericInput("tau","Characterstic Thickness, τ (mm):", 0.65, step=0.1, min = 0.5, max = 5),
      bsTooltip(id="tau", 
                title="See 'input parameters' section in the About tab", 
                placement = "bottom", 
                trigger = "hover",
                options = NULL),
      
      HTML("<br>"),
      actionButton("go", "Calculate"),
      HTML("<br><br>")
    ),
    
    # main panel layout
    mainPanel(
      h3('SHARP-ML™'),
      tabsetPanel(
        #tabPanel(title="K-values.json",
        #         verbatimTextOutput('kvalues_text')),
        tabPanel(title="Conversion factor table",
                  tableOutput('kvalues_table')),
        tabPanel(title="Chart", 
                 plotOutput("kvalues") %>%
                   withSpinner(type=1,color="#3475B6")),
        tabPanel(title="About",
                 withMathJax(includeMarkdown("about.Rmd"))
        )
      )
    )
  )
)


# Wrap your UI with secure_app
ui <- secure_app(ui,enable_admin = TRUE)

server <- function(input, output, session) {
  
  # check_credentials returns a function to authenticate users
  res_auth <- secure_server(
    check_credentials = check_credentials(
      db= "credentials/database.sqlite",
      passphrase = NULL
    )
  )
  
  output$auth_output <- renderPrint({
    reactiveValuesToList(res_auth)
  })
  
  formatted_predictions <- eventReactive(
    input$go,{
      # Create a Progress object
      progress <- shiny::Progress$new(min=0,max=4)
      progress$set(message = "", value = 0)
      
      # progress function
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value/progress$getMax()
        }
        progress$set(value = value, detail = detail)
      }
      
      # load model, make predictions, and format
      make_predictions(
        geometry = input$geometry,
        rho = input$rho,
        tau = input$tau,
        ensemble=TRUE,
        updateProgress)
      
    }
  )
  
  # plot
  output$kvalues <- renderPlot({
    formatted_predictions() %>%
      plot_predictions(geometry = input$geometry,
                       rho = input$rho,
                       tau = input$tau)
  })
  
  # table of predictions
  output$kvalues_table <- renderTable(formatted_predictions() %>%
                                        select(`Islet size group (μm)`=bucket,
                                               `K_ieq (Standard Conversion Factors)`=kieq, 
                                               κ = k_mu),
                                      digits = 4)
  
  # # table of predictions
  # output$kvalues_text <- renderText(
  #   list(model = 'SHARP-ML',
  #        run_date = Sys.Date(),
  #        data = formatted_predictions()) %>% 
  #     toJSON(auto_unbox = TRUE) %>% 
  #     prettify()
  #   )
}


shinyApp(ui, server)


#rsconnect::deployApp(appName="sharp-ml",appTitle="SHARP-ML")