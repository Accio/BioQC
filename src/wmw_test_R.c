#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <math.h>

#include "stat_rank.h"

/*! \brief Wilcoxon-Mann-Whitney Test
 *  
 * \param indlist: A list of integers giving the index of gene sets
 * \param matrix: an expression matrix with features in rows and samples in columns
 * \param val_type: 0=U, 1=p(left), 2=p(right), 3=p(two.sided)
 */
#define NROW(x) INTEGER(GET_DIM((x)))[0]
#define NCOL(x) INTEGER(GET_DIM((x)))[1]

double do_wmw_test(const int* indices,
		   int n1,
		   const double* stats,
		   int n,
		   int val_type /* 0=U, 1=p(less), 2=p(right) */ ) {

  int i,j;
  int n2=n-n1;
  double irsum=0; // sum of rank of index
  static double tmp;
  double U, mu, sigma2,zlt, zut, plt, pgt;
  int ulen=0;

  DRankList list=createDRankList(stats, n);
  rankDRankList(list);

  // sum of the ranks of indexed elements
  for(i=0;i<n1;++i)
    irsum+=list->list[ indices[i] ]->rank;
  ulen=list->ulen;

  U=n1*n2+n1*(n1+1.0)*0.5-irsum; 
  if(val_type==0) {
    destroyDRankList(list);
    return(U);
  }
  
  mu=(double)n1*n2*0.5;
  sigma2=n1*n2*(n+1.0)/12;
  
  if(ulen!=n) { // has ties
    int* tbl=(int*)malloc(ulen * sizeof(int));
    int ncount=0;
    double prod=0;
    sortDRankList(list);
    for(i=0;i<n;i=j+1) {
      j=i;
      while(j<n-1 && (*(list->list[j+1]->vPtr)==*(list->list[j]->vPtr))) ++j;
      tbl[ncount++]=j-i+1;
    }
    for(i=0;i<ulen;++i)
      prod+=(0.0+tbl[i])/n*(tbl[i]+1)/(n+1)*(tbl[i]-1)/(n-1);
    sigma2*=(1-prod);
    free(tbl);
  }

  destroyDRankList(list);

  if(val_type==1) {
    zut=(U-0.5-mu)/sqrt(sigma2);
    pnorm_both(zut, &tmp, &plt, 1, 0);
    return(plt);
  } else if (val_type==2) {
    zlt=(U+0.5-mu)/sqrt(sigma2);
    pnorm_both(zlt, &pgt, &tmp, 0, 0);
    return(pgt);
  } else {
    error("Should not happen");
    return(-1);
  }
}
		   
SEXP wmw_test(SEXP indlist, SEXP matrix, SEXP val_type) {
  int i,j, ilen,k;
  int *ip;
  int type=INTEGER(val_type)[0];
  const int slen=NROW(matrix);
  double *sp=REAL(matrix);
  SEXP res;
  res=PROTECT(allocMatrix(REALSXP,
			  length(indlist), 
			  NCOL(matrix)));
  double *resp=REAL(res);
  for(i=0; i<NCOL(matrix);++i) {
    for(j=0;j<length(indlist);++j) {
      ip=INTEGER(VECTOR_ELT(indlist,j));
      ilen=length(VECTOR_ELT(indlist,j));
      resp[j+i*length(indlist)]=do_wmw_test(ip, ilen, sp, slen, type); // notice the pointer should be j+i*length(indlist) [row count + column count * nrow]
    }
    sp+=NROW(matrix);
  }
  
  UNPROTECT(1);
  return(res);
}
