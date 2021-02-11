library(httr)


# add genes to enrichr
ENRICHR_URL <- 'http://maayanlab.cloud/Enrichr/addList'

genes_str = c(
  'PHF14', 'RBM3', 'MSL1', 'PHF21A', 'ARL10', 'INSR', 'JADE2', 'P2RX7',
  'LINC00662', 'CCDC101', 'PPM1B', 'KANSL1L', 'CRYZL1', 'ANAPC16', 'TMCC1',
  'CDH8', 'RBM11', 'CNPY2', 'HSPA1L', 'CUL2', 'PLBD2', 'LARP7', 'TECPR2', 
  'ZNF302', 'CUX1', 'MOB2', 'CYTH2', 'SEC22C', 'EIF4E3', 'ROBO2',
  'ADAMTS9-AS2', 'CXXC1', 'LINC01314', 'ATF7', 'ATP5F1')

r <- POST(ENRICHR_URL, body = list(list = paste(genes_str, collapse="\n"),
                                   description = 'Example gene list'
                                   ))

jsonlite::fromJSON(content(r, "text"))[[2]] -> req_id


# enrich
cols_values <- "Rank, Term name, P-value, Z-score, Combined score, Overlapping genes, Adjusted p-value, Old p-value, Old adjusted p-value"
cols_values <- strsplit(cols_values, ", ")[[1]]


bg_list <- c("KEGG_2015", "KEGG_2019_Human")

get_enrich <- function(x) {
  ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/enrich'
  gene_set_library = x
  glue::glue("?userListId={req_id}&backgroundType={gene_set_library}") -> req_url
  r2 <- GET(
    paste0(ENRICHR_URL, req_url)
  )
  en_res <- jsonlite::fromJSON(content(r2, "text"))
  names(en_res) -> db
  
  as.data.frame(matrix(unlist(en_res[[1]]),nrow=length(en_res[[1]]),byrow=TRUE), stringsAsFactors = FALSE)

}



purrr::map_df(bg_list, get_enrich) %>%  View()


url = 'http://maayanlab.cloud/Enrichr/enrich'
db = "KEGG_2019_Human"
glue::glue("?userListId={req_id}&backgroundType={db}") -> url_req
r3 <- GET(
  paste0(url, url_req)
)
enr_res <- jsonlite::fromJSON(content(r3, "text"))
names(enr_res) -> db2

as.data.frame(matrix(unlist(enr_res[[1]]),nrow=length(enr_res[[1]]),byrow=TRUE), stringsAsFactors = FALSE)
              
330/36

lengths(enr_res[[1]])

(9 * 36)


library(enrichR)
library(dplyr)
library(sparkline) 
library(tidyverse)
library(reactable)
library(htmltools)

listEnrichrDbs() -> dbs_en
enrichr(c("AKT3", "AKT2"), databases = "GO_Biological_Process_2018") -> res
res$GO_Biological_Process_2018 -> table_ex

table_ex %>% mutate(Overlap = strsplit(Overlap, "/")) %>% select(Term, Overlap, Genes, P.value) -> data

data %>% 
  extract(Term, c("Term", "GOID"), "(.*)\\((.*)\\)") -> data


make_color_pal <- function(colors, bias = 1) {
  get_color <- colorRamp(colors, bias = bias)
  function(x) rgb(get_color(x), maxColorValue = 255)
}

orange_pal <- function(x) rgb(colorRamp(c("#ffe4cc", "#ff9500"))(x), maxColorValue = 255)

off_rating_color <- make_color_pal(c("#ff2700", "#f8fcf8", "#44ab43"), bias = 1.3)

orange_pal <- function(x) rgb(colorRamp(c("red", "white"))(x), maxColorValue = 255)


reactable(select(data,  -GOID, -Genes), columns = list(
  Overlap = colDef(align = "center",cell = function(value, index) {
    sparkline(data$Overlap[[index]], type = "pie")
  }),
  
  P.value = colDef(cell = function(value){
    scaled <- (value - min(data$P.value)) / (max(data$P.value) - min(data$P.value))
    color <- orange_pal(scaled)
    value<- signif(value, 2)
    div(class = "spi-rating", style = list(background = color), value)
  })
  ),
  
  details = function(index) {
  paste(data$GOID[[index]], data$Genes[[index]], sep = " - ")
  }

)



url <- 'http://amp.pharm.mssm.edu/clustergrammer/matrix_upload/'
hg <- upload_file("files/small_38x29_clustergrammer_matrix.txt")

res <- POST(url, body = list(file = hg))
content(res)

colnames()

table_ex %>% select(Term, P.value, Genes) %>% 
  separate_rows(Genes, sep = ";") %>% 
  pivot_wider(names_from = Genes, values_from = P.value) %>% 
  column_to_rownames("Term") %>% as.matrix() %>% log10() %>% `*`(-1) -> clust_mtx

write.table(clust_mtx,"files/clust_mtx.txt", sep = "\t", quote = F, col.names = NA)
hg2 <- upload_file("files/clust_mtx.txt")

res <- POST(url, body = list(file = hg2))
content(res) %>% class()

library(xml2)

read_xml("<p>This is some <b>text</b>. This is more.</p>")
xml_text(content(res))



data_both %>% separate_rows(Genes, sep = ";") %>% 
  left_join(FCC, by = c("Genes" = "ID")) %>% group_by(Term) %>%  
  mutate(Up  = sum(FC > 0),
         Down = sum(FC < 0)) %>% ungroup() %>% 
  select(-FC) %>% group_by(Term) %>% mutate(Genes = paste0(Genes, collapse = ";")) %>% 
    ungroup() %>% distinct() %>% mutate(Dir = paste(Up, Down, sep = "_"),
                                        Dir = strsplit(Dir, "_")
                                        ) -> data


red_pal <- function(x) rgb(colorRamp(c("red", "white"))(x), maxColorValue = 255)

data %>% mutate(Down = ifelse(Term == "regulation of superoxide metabolic process", 4, Down)) %>% View()

reactable(select(data,  Term, Overlap,Dir, P.value), columns = list(
  Overlap = colDef(align = "center",cell = function(value, index) {
    sparkline(data$Overlap[[index]], type = "pie")
  }),
  
  P.value = colDef(cell = function(value){
    scaled <- (value - min(data$P.value)) / (max(data$P.value) - min(data$P.value))
    color <- red_pal(scaled)
    value<- signif(value, 2)
    div(class = "spi-rating", style = list(background = color), value)
  }),
  
  Dir = colDef(cell = function(values) {
    sparkline(values, type = "bar", chartRangeMin = 0, colorMap=c("#008000", "#800000"))
  })
),

details = function(index) {
  paste(data$GOID[[index]], data$Genes[[index]], sep = " - ")
}

)
