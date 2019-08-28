mypkg=dcibioinformatics
mypkgver=0.1
mypkggz=${mypkg}_${mypkgver}.tar.gz

all: pkg check

pkg:
	Rscript -e "Rcpp::compileAttributes(pkg = '${mypkg}')"
	Rscript -e "devtools::document(pkg = '${mypkg}')"

check:
	R CMD build ${mypkg}
	R CMD check ${mypkggz}

