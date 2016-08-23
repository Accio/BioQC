################################################################################
##
##  Makefile
##      Author: Jitao David Zhang <jitao_david.zhang@roche.com>
##	BEDA TRS, pRED, Hoffmann-La Roche AG
##      Description: Makefile for building distributions etc.
##                   the Makefile provides the following targets:
##                   
##                   - make install  calls R CMD INSTALL
##                   - make check    calls R CMD check (with RUnit)
##                   - make dist     calls R CMD build
##
################################################################################
## conditional: choose R version depending on the BICOSN value

R=R
PKG_ROOT_DIR=`pwd`

PKG=BioQC
PKG_VERSION=`/bin/awk 'BEGIN{FS=":"}{if ($$1=="Version") {gsub(/ /, "",$$2);print $$2}}' ${PKG}/DESCRIPTION`

PKG_SRC_DIR=$(PKG_ROOT_DIR)/src

dist:	clean
	@echo '====== Building Distribution ======'
	@(cd ..; ${R} CMD build ${DISTADD} $(PKG) )
	@echo '====== Building finished ======'
	@echo ' '

install: dist
	@echo '====== Installing Package ======'
	(cd ..; ${R} CMD INSTALL ${PKG}_${PKG_VERSION}.tar.gz) 
	@echo '====== Installing finished ======'
	@echo ' '

install-test: clean
	${R} CMD INSTALL ../${PKG}; R -e "library(testthat); test_dir('./tests')"

check:	dist
	@echo '====== Checking Package ======'
	@(cd ..; ${R} CMD check ${CHECKADD} ${PKG}_${PKG_VERSION}.tar.gz)
	@echo '====== Checking finished ======'
	@echo ' '

envcheck: dist
	@echo '====== Checking Package w/o Environmental Vars ======'
	@(cd ..; env -i BIOINFOCONFDIR=${BIOINFOCONFDIR} PATH="/usr/bin/:/usr/local/bin:/bin/:/usr/bin/:/usr/sbin/:/usr/local/bin/:/usr/X11R6/bin:/opt/oracle/client/10/run_1/bin:/usr/kerberos/bin:" LD_LIBRARY_PATH="/homebasel/beda/zhangj83/libs" ${R} CMD check ${CHECKADD} ${PKG}_${PKG_VERSION}.tar.gz) 
	@echo '====== Checking finished ======'

clean:
	@echo '====== Cleaning Package ======'
	@(rm -f $(PKG_SRC_DIR)/*.o $(PKG_SRC_DIR)/*.so $(PKG_SRC_DIR)/*.dll $(PKG_SRC_DIR)/*.rds)
	@(find . -type f -name "*~" -exec rm '{}' \;)
	@(find . -type f -name ".Rhistory" -exec rm '{}' \;)
	@echo ' '
