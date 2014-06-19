#include <fstream>
#include <Rcpp.h>

#include "read_gmt.h"

using std::string;
using std::vector;
using std::ifstream;
using Rcpp::List;
using Rcpp::Named;

RcppExport SEXP read_gmt(SEXP filename) {
  string fname=Rcpp::as<string>(filename);
  ifstream in(fname.c_str());
  if(in) {
    string str;
    vector<List> gmtlist;
    int line=0;
    while(std::getline(in, str)) {
      line++;
      GmtItem git=GmtItem(str);
      if(git.isValid()) {
	List gitl=List::create(Named("name")=git.name(),
			       Named("desc")=git.desc(),
			       Named("genes")=git.genes());
	gmtlist.push_back(gitl);
      } else {
	REprintf("[Warning: invalid GMT file] Skipping line %d:'%s'\n", line, str.c_str());
      }
    }
    return Rcpp::wrap(gmtlist);
  } else {
    return R_NilValue;
  }
}
