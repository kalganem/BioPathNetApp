
library(shiny)
library(shinyjs)
library(shinyBS)
library(reactable)


shinyUI(navbarPage(title = "BioPathNet",
                   theme = "style/style.css",
                   footer = includeHTML("footer.html"),
                   fluid = TRUE, 
                   collapsible = TRUE,
                   

                   tabPanel("Home",
                            includeHTML("home.html"),
                            tags$script(src = "plugins/scripts.js"),
                            tags$head(
                              tags$link(rel = "stylesheet", 
                                        type = "text/css", 
                                        href = "plugins/font-awesome-4.7.0/css/font-awesome.min.css"),
                              tags$link(rel = "icon", 
                                        type = "image/png", 
                                        href = "images/logo_icon.png")
                            )
                   ),
                   

                   tabPanel("Input",
                            fileInput("inupt_file", label = "Upload your DGE file: ",
                                      accept = c(".csv")
                            ),
                            radioButtons("input_rank", "Rank By: ", choices = c("FC", "p.value")),
                            actionButton("ex_btn", "Example"),
                            dataTableOutput("preview_tbl"),
                            dataTableOutput("ranked_preview_tbl")
                            
                   ),
                   

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
                   
                   tabPanel("Enrichr", 
                            
                            tabsetPanel(#"output_panels",
                              tabPanel("Options", 
                                       radioButtons("enrichr_gene_list", "Select Genes: ", choices = c("Up", "Down", "Both"), selected = "Both"),
                                       selectInput("enrichr_db_list", label = "Select Pathway: ", choices = c("GO_Biological_Process_2018"),
                                                   selected = "GO_Biological_Process_2018"
                                       ),
                                       sliderInput("enrichr_num", label = "Number of Genes: ", min = 50,  max = 1000, value = 500, step=50),
                                       actionButton("enrichr_btn", label = "Run Enrichr")
                              ),
                              tabPanel("Table", 
                                       
                                       fluidRow(
                                         shinydashboard::box(title = "Up", status = "primary", reactableOutput("enrichr_tbl_up")),
                                         shinydashboard::box(title = "Down", status = "warning", reactableOutput("enrichr_tbl_down"))
                                       ),
                                       
                                       shinydashboard::box(title= "Both", status = "warning", reactableOutput("enrichr_tbl_both"), width = 12)
                                       
                              ),
                              tabPanel("Heatmap", 
                                       uiOutput("enrichr_htmp")
                              )
                            )
                            
                            )
                   
                   
))