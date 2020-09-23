#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(fgsea)


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("BioPathNet"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput("inupt_file", label = "Upload your DEG file: ",
                      accept = c(".csv")
                      )
        ),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(#"output_panels",
                        tabPanel("Preview", dataTableOutput("preview_tbl")),
                        tabPanel("GSEA",
                                 dataTableOutput("ranked_preview_tbl"),
                                 actionButton("gsea_btn", label = "Run GSEA"),
                                 dataTableOutput("gsea_table"),
                                 fluidRow(
                                     column(6,
                                            plotOutput("up_gsea")
                                     ),
                                     column(6,
                                            plotOutput("down_gsea")
                                     )
                                 ),
                                 
                                 selectInput("gsea_path_sel", label = "Select Pathway: ", choices = list()),
                                 plotOutput("ind_gsea_plot")
                                 
                                 
                                 )
                        )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    gsea_sel <- reactiveValues()
    
    # output$frame <- renderUI({
    #     test <- "https://maayanlab.cloud/clustergrammer/viz/5f5a506251ad21000e214918/small_38x29_clustergrammer_matrix.txt"
    #     my_test <- tags$iframe(src=test, height=600, width=700)
    #     my_test
    # })
    
    data_set0 <- reactive({
        req(input$inupt_file)
        inFile <- input$inupt_file
        x<-read.csv(inFile$datapath, header=T,stringsAsFactors = F)
        x
        
    })

    output$preview_tbl <- renderDataTable({
        data_set0() %>% head()
    })
    
    output$ranked_preview_tbl <- renderDataTable({
        data_set0() %>% distinct(ID, .keep_all = T) %>% 
            arrange(desc(FC)) %>% mutate(rank = row_number()) %>%
            select(ID, FC, rank) -> ranked_tbl
        
        gsea_sel$table <- ranked_tbl %>% select(ID, FC)
        
        ranked_tbl
    })
    
    
    observeEvent(input$gsea_btn, {
        
        gsea_pathway <- gmtPathways("./files/Rat_GO_AllPathways_with_GO_iea_April_01_2019_symbol.gmt")
        
        ranks <- setNames(gsea_sel$table$FC, gsea_sel$table$ID)
        
        print(head(ranks))
        
        print("done")
        
        fgseaRes <- fgsea(gsea_pathway, ranks, minSize=15, maxSize=500, nperm=1000)
        
        output$gsea_table <- renderDataTable({
            fgseaRes

        })
        # 
        topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
        topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
        topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

        output$up_gsea <- renderPlot({
            plotGseaTable(gsea_pathway[topPathwaysUp], ranks, fgseaRes,
                          gseaParam=0.5)
        })

        output$down_gsea <- renderPlot({
            plotGseaTable(gsea_pathway[topPathwaysDown], ranks, fgseaRes,
                          gseaParam=0.5)
        })
        
        updateSelectInput(session, "gsea_path_sel", choices = fgseaRes$pathway)
        
        output$ind_gsea_plot <- renderPlot({
            req(input$gsea_path_sel)
            
            plotEnrichment(gsea_pathway[[input$gsea_path_sel]],
                           ranks)
            
        })
        
        
        
    })
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
