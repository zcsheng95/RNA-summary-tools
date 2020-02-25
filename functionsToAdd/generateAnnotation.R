#### Generate fake anotation genesymbol
library(random)
set.seed(1234)
ngene <- 2000
genesymbol <- as.character(randomStrings(n=2000, len=5, digits=TRUE, upperalpha=TRUE,
                        loweralpha=TRUE, unique=TRUE, check=TRUE))
geneid <- paste0("gene", 1:ngene,".1")
class <- "pseudo gene"

anno <- data.frame(Geneid = geneid, GeneSymbol=genesymbol, Class = class)
rootdir <- "/home/zs68@dhe.duke.edu/Documents/dcibioinformatics/dcibioinformatics/inst/extdata"
write.table(anno, file = file.path(rootdir,"annotation.txt"),quote = F, col.names = T, sep = '\t')