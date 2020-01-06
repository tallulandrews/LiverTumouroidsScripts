# Read in go enrichments
# shorten term desc
# olap btw terms
# cluster & terms
x <- read.table("CCA1_SC3_marker_GOrich.txt")




# shorted name rules
x$term.name <- sub("regulation", "reg", x$term.name)
x$term.name <- sub("regulates", "reg", x$term.name)
x$term.name <- sub("positive", "pos", x$term.name)
x$term.name <- sub("negative", "neg", x$term.name)
x$term.name <- sub("membrane", "mmb", x$term.name)
x$term.name <- sub("cellular", "", x$term.name)
x$term.name <- sub("of ", "", x$term.name)
x$term.name <- sub("to ", "", x$term.name)
x$term.name <- sub("amino acid", "AA", x$term.name)
x$term.name <- sub("development", "devo", x$term.name)
x$term.name <- sub("Uncertain", "U", x$term.name)
x$term.name <- sub("Medium", "M", x$term.name)
x$term.name <- sub("Supported", "S", x$term.name)
x$term.name <- sub("High", "H", x$term.name)
x$term.name <- sub("Confirmed", "C", x$term.name)
x$term.name <- sub("cells", "", x$term.name)
