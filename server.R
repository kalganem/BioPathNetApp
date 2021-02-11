library(shiny)
library(tidyverse)
library(fgsea)
library(reactable)
library(htmltools)
library(sparkline)
library(enrichR)
library(httr)


red_pal <- function(x) rgb(colorRamp(c("red", "white"))(x), maxColorValue = 255)



server <- function(input, output, session) {
  
  input_file <- reactiveValues()
  
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
    input_file$data <- read.csv("data/example_ds.csv", header=T,stringsAsFactors = F)
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
      gsea_pathway <- gmtPathways("./data/gsea_genesets/c5.bp.v7.1.symbols.gmt")
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
  
  observeEvent(input$enrichr_btn, {
    req(input_file$ranked)
    withProgress(message = "Running Enrichr", value = 0, {
      input_file$ranked %>% arrange(desc(abs(FC))) %>% 
        head(input$enrichr_num) %>% pull(ID) %>% 
        toupper() -> enrichr_genes_both
      
      input_file$ranked %>% arrange(FC) %>% 
        head(input$enrichr_num) %>% pull(ID) %>% 
        toupper() -> enrichr_genes_down
      
      input_file$ranked %>% head(input$enrichr_num) %>% 
        pull(ID) %>% toupper() -> enrichr_genes_up
      
      gene_list <- list(both = enrichr_genes_both, down = enrichr_genes_down,
                        up = enrichr_genes_up
      )
      
      input_file$ranked -> FCC
      
      #listEnrichrDbs() -> dbs_en
      #enrichr(enrichr_genes_up, databases = input$enrichr_db_list) -> res
      
      incProgress(0.3)
      res_full <- map(gene_list, enrichr, databases = input$enrichr_db_list)
      
      
      res_full$up$GO_Biological_Process_2018 %>% head(15) -> table_ex_up
      res_full$down$GO_Biological_Process_2018 %>% head(15) -> table_ex_down
      res_full$both$GO_Biological_Process_2018 %>% head(15) -> table_ex_both
      
      incProgress(0.7)
      table_ex_up %>% mutate(Overlap = strsplit(Overlap, "/")) %>% 
        select(Term, Overlap, Genes, P.value) %>% 
        extract(Term, c("Term", "GOID"), "(.*)\\((.*)\\)") -> data_up
      
      table_ex_down %>% mutate(Overlap = strsplit(Overlap, "/")) %>% 
        select(Term, Overlap, Genes, P.value) %>% 
        extract(Term, c("Term", "GOID"), "(.*)\\((.*)\\)") -> data_down
      
      table_ex_both %>% mutate(Overlap = strsplit(Overlap, "/")) %>% 
        select(Term, Overlap, Genes, P.value) %>% 
        extract(Term, c("Term", "GOID"), "(.*)\\((.*)\\)") %>% 
        separate_rows(Genes, sep = ";") %>% 
        left_join(FCC, by = c("Genes" = "ID")) %>% group_by(Term) %>%  
        mutate(Up  = sum(FC > 0),
               Down = sum(FC < 0)) %>% ungroup() %>% 
        select(-FC) %>% group_by(Term) %>% mutate(Genes = paste0(Genes, collapse = ";")) %>% 
        ungroup() %>% distinct() %>% mutate(Dir = paste(Up, Down, sep = "_"),
                                            Dir = strsplit(Dir, "_")
        ) -> data_both
      
      
    })
    output$enrichr_tbl_up <- renderReactable({
      reactable(select(data_up,  -GOID, -Genes), columns = list(
        Overlap = colDef(align = "center",cell = function(value, index) {
          sparkline(data_up$Overlap[[index]], type = "pie", sliceColors=c("#c90000", "#999494"))
        }),
        
        P.value = colDef(cell = function(value){
          scaled <- (value - min(data_up$P.value)) / (max(data_up$P.value) - min(data_up$P.value))
          color <- red_pal(scaled)
          value<- signif(value, 2)
          div(class = "spi-rating", style = list(background = color), value)
        })
      ),
      
      details = function(index) {
        paste(data_up$GOID[[index]], data_up$Genes[[index]], sep = " - ")
      }
      
      )
    })
    
    output$enrichr_tbl_down <- renderReactable({
      reactable(select(data_down,  -GOID, -Genes), columns = list(
        Overlap = colDef(align = "center",cell = function(value, index) {
          sparkline(data_down$Overlap[[index]], type = "pie", sliceColors=c("#c90000", "#999494"))
        }),
        
        P.value = colDef(cell = function(value){
          scaled <- (value - min(data_down$P.value)) / (max(data_down$P.value) - min(data_down$P.value))
          color <- red_pal(scaled)
          value<- signif(value, 2)
          div(class = "spi-rating", style = list(background = color), value)
        })
      ),
      
      details = function(index) {
        paste(data_down$GOID[[index]], data_down$Genes[[index]], sep = " - ")
      }
      
      )
    })
    
    output$enrichr_tbl_both <- renderReactable({
      reactable(select(data_both,  Term, Overlap,Dir, P.value), columns = list(
        Overlap = colDef(align = "center",cell = function(value, index) {
          sparkline(data_both$Overlap[[index]], type = "pie", sliceColors=c("#c90000","#999494"))
        }),
        
        P.value = colDef(cell = function(value){
          scaled <- (value - min(data_both$P.value)) / (max(data_both$P.value) - min(data_both$P.value))
          color <- red_pal(scaled)
          value<- signif(value, 2)
          div(class = "spi-rating", style = list(background = color), value)
        }),
        
        Dir = colDef(cell = function(values) {
          sparkline(values, type = "bar", chartRangeMin = 0, colorMap=c("#008000", "#c90000"))
        })
      ),
      
      details = function(index) {
        paste(data_both$GOID[[index]], data_both$Genes[[index]], sep = " - ")
      }
      
      )
    })
    
    table_ex_up %>% select(Term, P.value, Genes) %>%
      separate_rows(Genes, sep = ";") %>%
      pivot_wider(names_from = Genes, values_from = P.value) %>%
      column_to_rownames("Term") %>% as.matrix() %>% log10() %>% `*`(-1) -> clust_mtx
    
    write.table(clust_mtx,"data/clust_mtx.txt", sep = "\t", quote = F, col.names = NA)
    hg2 <- upload_file("data/clust_mtx.txt")
    
    url <- 'http://amp.pharm.mssm.edu/clustergrammer/matrix_upload/'
    res <- POST(url, body = list(file = hg2))
    
    clust_link <- xml2::xml_text(content(res))
    content(res)
    print(clust_link)
    
    output$enrichr_htmp <- renderUI({
      req(clust_link)
      #test <- "https://maayanlab.cloud/clustergrammer/viz/5f5a506251ad21000e214918/small_38x29_clustergrammer_matrix.txt"
      my_test <- tags$iframe(src=clust_link, style='width:100vw;height:100vh;')
      my_test
    })
    
  })
  
  # iLINCS Server ----
  # Final Results Server ----
  
}