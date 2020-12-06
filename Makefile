## AUTOMATICALLY GENERATED FROM TEMPLATE (Mo Jul  8 20:50:58 CEST 2019). DO NOT EDIT IT MANUALLY!
################################################################################
##
##  Makefile
##      Author: Jitao David Zhang <jitao_david.zhang@roche.com>
##	    F. Hoffmann-La Roche AG
##      Description: Makefile for building distributions etc.
##
################################################################################
R:=R

install_deps:
	@echo '====== install dependencies ======'	
	@(${R} -e 'install.packages("remotes", repos="https://cran.rstudio.com")')
	@(${R} -e 'remotes::install_cran(c("devtools", "roxygen2", "covr"), repos="https://cran.rstudio.com")')
	@(${R} -e 'deps <- remotes::dev_package_deps(dependencies = NA); remotes::install_deps(dependencies = TRUE, upgrade="never"); if (!all(deps$$package %in% installed.packages())) { message("missing: ", paste(setdiff(deps$$package, installed.packages()), collapse=", ")); q(status = 1, save = "no")}')

roxygenise:
	@echo '====== roxygenize ======'	
	@(${R} -q -e "library(devtools);load_all();document('.')")
	@echo ' '

test:
	@echo '====== test ======'
	@(${R} -q -e "library(devtools);test('.')")
	@echo 

doVignettes:
	@echo "====== vignettes ======"
	@(${R} -q -e "library(devtools); devtools::build_vignettes()")
	@echo ' '

build: roxygenise
	@echo '====== Building Distribution ======'
	@(${R} -q -e "library(devtools); devtools::build()")
	@echo '====== Building finished ======'
	@echo ' '

install: roxygenise
	@echo '====== Installing Package ======'
	@(${R} -q -e "library(devtools); devtools::install(reload=FALSE, quick=FALSE, build=TRUE, upgrade=FALSE)")
	@echo '====== Installing finished ======'
	@echo ' '

check: roxygenise
	@echo '====== Checking Package ======'
	@(${R} -q -e "library(devtools);check('.', check_dir=\"..\")")
	@echo '====== Checking finished ======'
	@echo ' '

docs: roxygenise
	@echo '====== Building documentation with pkgdown ======'
	@(${R} -e "pkgdown::build_site()") 
	@echo '====== Building documentation finished ======'
	@echo ' '

clean:
	@echo '====== Cleaning Package ======'
	@(rm -f src/*.o src/*.so src/*.dll src/*.rds)
	@(find . -type f -name "*~" -exec rm '{}' \;)
	@(find . -type f -name ".Rhistory" -exec rm '{}' \;)
	@echo ' '
