all:BioQC.so ../tests/dyn.R
	R-devel -f ../tests/dyn.R

BioQC.so: stat_rank.c wmw_test_R.c objs
	gcc -std=gnu99 -I/SOFT/bi/apps/R/R-2.16.0/include -DDEBUG `Rscript -e 'Rcpp:::CxxFlags()'` -I/usr/local/include    -fpic  -g -O2  -c stat_rank.c -o stat_rank.o
	gcc -std=gnu99 -I/SOFT/bi/apps/R/R-2.16.0/include -DDEBUG `Rscript -e 'Rcpp:::CxxFlags()'` -I/usr/local/include    -fpic  -g -O2  -c wmw_test_R.c -o wmw_test_R.o
	g++ -shared -L/usr/local/lib64 -o BioQC.so bioqc.o gini.o read_gmt.o stat_rank.o wmw_test_R.o -L/apps/bi/apps/R/R-2.16.0/library/Rcpp/lib -lRcpp -Wl,-rpath,/apps/bi/apps/R/R-2.16.0/library/Rcpp/lib


objs:bioqc.o gini.o read_gmt.o

bioqc.o:bioqc.c
	gcc -std=gnu99 -I/SOFT/bi/apps/R/R-2.16.0/include -DDEBUG `Rscript -e 'Rcpp:::CxxFlags()'` -I/usr/local/include    -fpic  -g -O2  -c bioqc.c -o bioqc.o

gini.o:gini.c
	gcc -std=gnu99 -I/SOFT/bi/apps/R/R-2.16.0/include -DDEBUG `Rscript -e 'Rcpp:::CxxFlags()'` -I/usr/local/include    -fpic  -g -O2  -c gini.c -o gini.o

read_gmt.o:read_gmt.cpp
	g++  -I/SOFT/bi/apps/R/R-2.16.0/include -DDEBUG `Rscript -e 'Rcpp:::CxxFlags()'` -I/usr/local/include    -fpic  -g -O2  -c read_gmt.cpp -o read_gmt.o
