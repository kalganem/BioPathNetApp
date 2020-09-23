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




