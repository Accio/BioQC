all:BioQC.so ../tests/dyn.R
	R-devel -f ../tests/dyn.R

BioQC.so: stat_rank.c wmw_test_R.c
	gcc -std=gnu99 -I/SOFT/bi/apps/R/R-2.16.0/include -DDEBUG `Rscript -e 'Rcpp:::CxxFlags()'` -I/usr/local/include    -fpic  -g -O2  -c stat_rank.c -o stat_rank.o
	gcc -std=gnu99 -I/SOFT/bi/apps/R/R-2.16.0/include -DDEBUG `Rscript -e 'Rcpp:::CxxFlags()'` -I/usr/local/include    -fpic  -g -O2  -c wmw_test_R.c -o wmw_test_R.o
	g++ -shared -L/usr/local/lib64 -o BioQC.so bioqc.o gini.o read_gmt.o stat_rank.o wmw_test_R.o -L/apps/bi/apps/R/R-2.16.0/library/Rcpp/lib -lRcpp -Wl,-rpath,/apps/bi/apps/R/R-2.16.0/library/Rcpp/lib


