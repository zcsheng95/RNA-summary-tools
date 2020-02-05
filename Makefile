mypkg=dcibioinformatics
mypkgver=1.0
mypkggz=${mypkg}_${mypkgver}.tar.gz

all: pkg check

pkg:
	# Rscript -e "Rcpp::compileAttributes(pkg = '${mypkg}')"
	Rscript -e "devtools::document(pkg = '${mypkg}')"

check:
	R CMD build ${mypkg}
	R CMD check ${mypkggz}

