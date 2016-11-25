library(shiny)
EMSaov.env<<-new.env()

shinyUI(fluidPage(
      shiny::headerPanel("Shiny Application for ANOVA with EMS"),
      shiny::fileInput("outputfile",label="File input"),
      shiny::br(),
      shiny::wellPanel(  
        shiny::fluidRow(shiny::column(3,shiny::uiOutput("choose_Yvar"))),
        shiny::fluidRow(
          shiny::column(2,shiny::uiOutput("choose_Xvar")),
          shiny::column(2,shiny::uiOutput("choose_type")),
          shiny::column(2,shiny::uiOutput("choose_level")),
          shiny::column(2,shiny::uiOutput("choose_nested")),
          shiny::column(2,shiny::uiOutput("choose_split"))
        ),
        shiny::submitButton("Submit")
      ),
      shiny::hr(),
      shiny::mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel("EDA-main effect",shiny::plotOutput("EDA1")),
          shiny::tabPanel("EDA-interaction",shiny::plotOutput("EDA2")),     
          shiny::tabPanel("ANOVA table",shiny::tableOutput("result1"),
                          shiny::p(paste("Signif. codes : <0.0001 '***'",
                                 "0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1"))),
          shiny::tabPanel("ANOVA table with Approx. F",
                          shiny::tableOutput("result2"),
                          shiny::p(paste("Signif. codes : <0.0001 '***'",
                                  "0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1"))),
          shiny::tabPanel("Pooled ANOVA",
                          shiny::uiOutput("choose_ANOVA"),
                          shiny::submitButton("Submit1"),
                          shiny::tableOutput("result3"))
        )    
      )
    )
)