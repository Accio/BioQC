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

SEXP wmw_test(SEXP indlist, SEXP matrix, SEXP val_type) {
  int i,j, n1,n2, k, m;
  int *ip;
  const int type=INTEGER(val_type)[0];
  const int slen=NROW(matrix);
  double *sp=REAL(matrix);
  SEXP res;

  DRankList list;
  double irsum; // sum of rank of index
  double tmp;
  double U;
  double mu, sigma2;
  double zlt, zut, plt, pgt;

  // for ties
  double tiecoef;

  res=PROTECT(allocMatrix(REALSXP,
			  length(indlist), 
			  NCOL(matrix)));
  double *resp=REAL(res);
  
  for(i=0; i<NCOL(matrix);++i) {
    list=createDRankList(sp, slen);    
    rankDRankList(list);

    if(list->ulen!=slen) {
      int* tbl=(int*)malloc(list->ulen * sizeof(int));
      int ncount=0;
      tiecoef=0;
      sortDRankList(list);
      for(k=0;k<slen;k=m+1) {
	m=k;
	while(m<slen-1 && (*(list->list[m+1]->vPtr)==*(list->list[m]->vPtr))) ++m;
	tbl[ncount++]=m-k+1;
      }
      for(k=0;k<list->ulen;++k)
	tiecoef+=(0.0+tbl[k])/slen*(tbl[k]+1)/(slen+1)*(tbl[k]-1)/(slen-1);
      tiecoef=1-tiecoef;
      free(tbl);
      rankDRankList(list);
    } else {
      tiecoef=1;
    }

    for(j=0;j<length(indlist);++j) {
      ip=INTEGER(VECTOR_ELT(indlist,j));
      n1=length(VECTOR_ELT(indlist,j));
      n2=slen-n1;
      irsum=0.0;
      for(k=0;k<n1;++k) 
	irsum+=list->list[ip[k]]->rank;

      U=n1*n2+n1*(n1+1.0)*0.5-irsum; 
      if(type==0) {
	resp[j+i*length(indlist)]=U;
      } else {
	mu=(double)n1*n2*0.5; // mu=n1*n2*0.5
	sigma2=n1*n2*(slen+1.0)/12*tiecoef; //sigma2 = n1*n2*(n+1)/12*tiecoef

	if(type==1) {
	  zut=(U-0.5-mu)/sqrt(sigma2);
	  pnorm_both(zut, &tmp, &plt, 1, 0);
	  resp[j+i*length(indlist)]=plt;
	} else if (type==2) {
	  zlt=(U+0.5-mu)/sqrt(sigma2);
	  pnorm_both(zlt, &pgt, &tmp, 0, 0);
	  resp[j+i*length(indlist)]=pgt;
	} else {
	  Rprintf("Unrocognized val_type. Should not happen\n");
	}
      }
      //resp[j+i*length(indlist)]=do_wmw_test(ip, n1, sp, slen, type); 
    }

    destroyDRankList(list);
    sp+=NROW(matrix);
  }
  
  UNPROTECT(1);
  return(res);
}
