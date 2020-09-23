
library(httr)


url <- 'http://amp.pharm.mssm.edu/clustergrammer/matrix_upload/'



hg <- upload_file("files/small_38x29_clustergrammer_matrix.txt")

POST("http://httpbin.org/post", body = list(y = citation))


res <- POST(url, body = list(file = hg))
res %>% View()

content(res)


test_mx <- matrix(c(84,58,38,63), ncol = 2, nrow = 2,
                  dimnames = list(c("gene1", "gene2"),
                                  c("sample1", "sample2")
                                  )
                  )


res <- POST(url, body = list(data = test_mx))
res %>% View()

content(res)
