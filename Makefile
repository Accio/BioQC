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
ifneq ($(BICOSN), bas)
	R:= /SOFT/bi/apps/R/devel/trunk/bin/R
	CHECKADD:= ${CHECKADD} --no-manual ## for envcheck
else
	R:= R
	CHECKADD:= ${CHECKADD} --no-latex
endif 

PKG          := $(shell awk 'BEGIN{FS=":"}{if ($$1=="Package") {gsub(/ /, "",$$2);print $$2}}' DESCRIPTION)
PKG_VERSION  := $(shell awk 'BEGIN{FS=":"}{if ($$1=="Version") {gsub(/ /, "",$$2);print $$2}}' DESCRIPTION)


PKG_ROOT_DIR := $(shell pwd)
PKG_SRC_DIR := $(PKG_ROOT_DIR)/src

install: 
	@echo '====== Installing Package ======'
	@(cd ..; ${R} CMD INSTALL $(PKG))
	@echo '====== Installing finished ======'
	@echo ' '

check:	dist
	@echo '====== Checking Package ======'
	@(cd ..; ${R} CMD check ${CHECKADD} ${PKG}_${PKG_VERSION}.tar.gz)
	@echo '====== Checking finished ======'
	@echo ' '

envcheck: dist
	@echo '====== Checking Package w/o Environmental Vars ======'
	@(cd ..; env -i BIOINFOCONFDIR=${BIOINFOCONFDIR} PATH="/usr/bin/:/usr/local/bin:/bin/:/usr/bin/:/usr/sbin/:/usr/local/bin/:/usr/X11R6/bin:/opt/oracle/client/10/run_1/bin:/usr/kerberos/bin:" LD_LIBRARY_PATH="/homebasel/beda/zhangj83/libs" ${R} CMD check ${CHECKADD} ${PKG}_${PKG_VERSION}.tar.gz) 
	@echo '====== Checking finished ======'

dist:	clean
	@echo '====== Building Distribution ======'
	@(cd ..; ${R} CMD build $(PKG))
	@echo '====== Building finished ======'
	@echo ' '

clean:
	@echo '====== Cleaning Package ======'
	@(rm -f $(PKG_SRC_DIR)/*.o $(PKG_SRC_DIR)/*.so)
	@(find . -type f -name "*~" -exec rm '{}' \;)
	@(find . -type f -name ".Rhistory" -exec rm '{}' \;)
	@echo ' '
