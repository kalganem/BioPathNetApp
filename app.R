
library(shiny)
library(tidyverse)
library(fgsea)


ui <- fluidPage(
    navbarPage("BioPathNet",

    # Overview UI Tab ----
    tabPanel("Overview"),
    
    # Input UI Tab ----
    tabPanel("Input", 
             
             fileInput("inupt_file", label = "Upload your DGE file: ",
                       accept = c(".csv")
             ),
             radioButtons("input_rank", "Rank By: ", choices = c("FC", "p.value")),
             actionButton("ex_btn", "Example"),
             dataTableOutput("preview_tbl"),
             dataTableOutput("ranked_preview_tbl")
             
             ), 
    
    # GSEA UI Tab ----
    tabPanel("GSEA", 
             
             tabsetPanel(#"output_panels",
                 tabPanel("Options", 
                          radioButtons("gsea_sp", "Species ", choices = c("Human", "Mouse", "Rat")),
                          sliderInput("gsea_itr", label = "Permutation", min = 100,  max = 1000, value = 500, step=100),
                          sliderInput("gsea_min", label = "Min", min = 5,  max = 50, value = 15, step=5),
                          sliderInput("gsea_max", label = "Max", min = 100,  max = 1000, value = 500, step=100),
                          
                          actionButton("gsea_btn", label = "Run GSEA")
                          ),
                 # tabPanel("Rank Table", 
                 #          dataTableOutput("ranked_preview_tbl")
                 #          ),
                 tabPanel("Table", 
                          dataTableOutput("gsea_table")
                          ),
                 tabPanel("Top",
                          fluidRow(
                              shinydashboard::box(title = "Up", solidHeader = T,
                              width = 8, collapsible = T,
                              plotOutput("up_gsea")),
                              shinydashboard::box(title = "Down", solidHeader = T,
                              width = 8, collapsible = T,
                              plotOutput("down_gsea"))
                          
                          )
                          ),
                 tabPanel("LeadingEdge",
                          selectInput("gsea_path_sel", label = "Select Pathway: ", choices = list()),
                          plotOutput("ind_gsea_plot")
                          )
             )
             
             ), 
    
    # Enrichr UI Tab ----
    tabPanel("Enrichr",
             tabsetPanel(#"output_panels",
                 tabPanel("Options", 
                          radioButtons("enrichr_list", "Select Genes: ", choices = c("Up", "Down", "Both"), selected = "Both"),
                          #sliderInput("enrichr_num", label = "Permutation", min = 50,  max = 1000, value = 500, step=50),
                          
                          actionButton("enrichr_btn", label = "Run Enrichr")
                 ),
                 tabPanel("Table", 
                          dataTableOutput("enrichr_tbl")
                 )
             )
             ), 
    
    # iLINCS UI Tab ----
    tabPanel("iLINCS"), 
    
    # Final Results UI Tab ----
    tabPanel("Final Results")

    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    input_file <- reactiveValues()
    
    # output$frame <- renderUI({
    #     test <- "https://maayanlab.cloud/clustergrammer/viz/5f5a506251ad21000e214918/small_38x29_clustergrammer_matrix.txt"
    #     my_test <- tags$iframe(src=test, height=600, width=700)
    #     my_test
    # })
    
    # Input Server ----
    observe({
        req(input$inupt_file)
        inFile <- input$inupt_file
        x<-read.csv(inFile$datapath, header=T,stringsAsFactors = F)
        input_file$data <- x
    })
    
    # example input file
    observeEvent(input$ex_btn, {
        print("ex_btn")
        input_file$data <- read.csv("files/15-1_BLA_ERvEE.csv", header=T,stringsAsFactors = F)
    })
    
    output$preview_tbl <- renderDataTable({
        input_file$data %>% head()
    })
    
    output$ranked_preview_tbl <- renderDataTable({
        req(input_file$data)
        input_file$data %>% distinct(ID, .keep_all = T) %>% 
            {if (input$input_rank == "p.value") arrange(.,p) else arrange(.,desc(FC))} %>% 
            mutate(rank = row_number()) %>%
            select(ID, FC,p, rank) -> ranked_tbl
        input_file$ranked <- ranked_tbl %>% select(ID, FC)
        ranked_tbl %>% head()
    })

    # GSEA Server ----
    observeEvent(input$gsea_btn, {
        req(input_file$ranked)
        print("run_gsea")
        withProgress(message = "Running GSEA", value = 0, {
            incProgress(0.3)
            gsea_pathway <- gmtPathways("./files/gsea_genesets/c5.bp.v7.1.symbols.gmt")
            ranks <- setNames(input_file$ranked$FC, input_file$ranked$ID)
            #print(head(ranks))
            incProgress(0.6)
            fgseaRes <- fgsea(gsea_pathway, ranks, minSize=input$gsea_min, maxSize=input$gsea_max, nperm=input$gsea_itr)
            incProgress(0.8)
            #fgseaRes <- mutate(fgseaRes1, pathway = str_remove(pathway, "GO_"), pathway = str_replace_all(pathway, "_", " "))
            
        })
        
        # gsea table
        output$gsea_table <- renderDataTable({
            fgseaRes %>% select(-leadingEdge) %>% mutate_if(is.numeric, round, 3) %>% 
                head(10)
        })
        
        # pulling top pathways
        topPathwaysUp <- fgseaRes %>% filter(ES > 0) %>% 
            arrange(pval) %>% head(10) %>% 
            pull(pathway)
        
        topPathwaysDown <- fgseaRes %>% filter(ES < 0) %>% 
            arrange(pval) %>% head(10) %>% 
            pull(pathway)
        
        # gsea ranked plot
        output$up_gsea <- renderPlot({
            plotGseaTable(gsea_pathway[topPathwaysUp], ranks, fgseaRes,
                          gseaParam=0.5)
        })

        output$down_gsea <- renderPlot({
            plotGseaTable(gsea_pathway[topPathwaysDown], ranks, fgseaRes,
                          gseaParam=0.5)
        })
        
        
        # gsea leading edge
        updateSelectInput(session, "gsea_path_sel", choices = fgseaRes$pathway)
        
        output$ind_gsea_plot <- renderPlot({
            req(input$gsea_path_sel)
            plotEnrichment(gsea_pathway[[input$gsea_path_sel]],
                           ranks)
        })
 
    })
    
    # Enrichr Server ----
    # iLINCS Server ----
    # Final Results Server ----
      
}

# Run the application 
shinyApp(ui = ui, server = server)
