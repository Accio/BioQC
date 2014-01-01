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
    while(getline(in, str)) {
      GmtItem git=GmtItem(str);
      List gitl=List::create(Named("name")=git.name(),
			     Named("desc")=git.desc(),
			     Named("genes")=git.genes());
      gmtlist.push_back(gitl);
    }
    return Rcpp::wrap(gmtlist);
  } else {
    return R_NilValue;
  }
}
