
library(httr)


url <- 'http://amp.pharm.mssm.edu/clustergrammer/matrix_upload/'
hg <- upload_file("files/small_38x29_clustergrammer_matrix.txt")

res <- POST(url, body = list(file = hg))
content(res)


test_mx <- matrix(c(84,58,38,63), ncol = 2, nrow = 2,
                  dimnames = list(c("gene1", "gene2"),
                                  c("sample1", "sample2")
                                  )
                  )

write.table(test_mx,"files/test_mtx.txt", sep = "\t", quote = F, col.names = NA)
hg2 <- upload_file("files/test_mtx.txt")


res <- POST(url, body = list(file = hg2))
res %>% View()

content(res)

hg
