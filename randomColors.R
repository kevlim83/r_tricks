#methods to generate distinct colors

#method1
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

#method2
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

#method3
library(randomcoloR)
n <- 20
palette <- distinctColorPalette(n)


