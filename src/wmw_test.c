#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <math.h>

#include "stat_rank.h"

#define MIN(x,y) ((x) > (y) ? (y) : (x))
#define NROW(x) INTEGER(GET_DIM((x)))[0]
#define NCOL(x) INTEGER(GET_DIM((x)))[1]
#define ABSLOG(x) fabs(log10( (x) ))

/*
 * void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p)
 
 * i_tail in {0,1,2} means: "lower", "upper", or "both" :
 * if(lower) return  *cum := P[X <= x]
 * if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
*/

/*! \brief Wilcoxon-Mann-Whitney Test
 *  
 * \param indlist: A list of integers giving the index of gene sets
 * \param matrix: an expression matrix with features in rows and samples in columns
 * \param val_type: 
 * \parblock
 * Define f(x)=abs(log10(x))
 * 0=p(greater), 1=p(less), 2=p(two.sided), 3=U,
 * 4=f(p(greater)),5=p(less),6=f(p(two.sided)), 7=p(greater)<p(less) ? f(p(two.sided)) : -f(p(two.sided))
 * \endparblock
 *
 * This implementation uses normal approximation, which works reasonably well if sample size is large (say N>=20)
 */
SEXP wmw_test(SEXP indlist, SEXP matrix, SEXP val_type) {
  int i,j, n1,n2, k, m;
  int *ip;
  const int type=INTEGER(val_type)[0];
  const int n=NROW(matrix); // total number of samples AND nrow of matrix
  double *matColPtr; // pointer to the current column of the matrix
  SEXP res;

  DRankList list;
  double irsum; // sum of rank of index
  double U;
  double mu, sigma2;
  double zval, pgt, plt, val;
  double tieCoef;
  double *resPtr;

  res=PROTECT(allocMatrix(REALSXP,
			  length(indlist), 
			  NCOL(matrix)));
  resPtr=REAL(res);
  matColPtr=REAL(matrix);
  
  for(i=0; i<NCOL(matrix);++i) {
    list=createDRankList(matColPtr, n);    
    rankDRankList(list);

    if(list->ulen!=n) {
      int* tbl=(int*)malloc(list->ulen * sizeof(int));
      int ncount=0;
      tieCoef=0;
      sortDRankList(list);
      for(k=0;k<n;k=m+1) {
	m=k;
	while(m<n-1 && (*(list->list[m+1]->vPtr)==*(list->list[m]->vPtr))) ++m;
	tbl[ncount++]=m-k+1;
      }
      for(k=0;k<list->ulen;++k)
	tieCoef+=(0.0+tbl[k])/n*(tbl[k]+1)/(n+1)*(tbl[k]-1)/(n-1);
      tieCoef=1-tieCoef;
      free(tbl);
      rankDRankList(list);
    } else {
      tieCoef=1.0;
    }

    for(j=0;j<length(indlist);++j) {
      ip=INTEGER(VECTOR_ELT(indlist,j));
      n1=length(VECTOR_ELT(indlist,j));
      n2=n-n1;
      irsum=0.0;
      for(k=0;k<n1;++k) 
	irsum+=list->list[ip[k]]->rank;

      U=n1*n2+n1*(n1+1.0)*0.5-irsum; 
      if(type==3) {
	val=U;
      } else {
	mu=(double)n1*n2*0.5; // NOT mu=n1*n2*0.5
	sigma2=n1*n2*(n+1.0)/12*tieCoef; //NOT sigma2 = n1*n2*(n+1)/12*tieCoef

	if(type==0 || type==4) { /* greater */
	  zval=(U+0.5-mu)/sqrt(sigma2); // z lower tail
	  pnorm_both(zval, &pgt, &plt, 0, 0);
	  val=type==0 ? pgt : ABSLOG(pgt);
	} else if (type==1 || type==5) { /* less */
	  zval=(U-0.5-mu)/sqrt(sigma2); // z higher tail
	  pnorm_both(zval, &pgt, &plt, 1, 0);
	  val=type==1 ? plt : log10(plt);
	} else if (type==2 || type==6 || type==7) { /* two sided*/
	  zval=(U-mu- (U>mu ? 0.5 : -0.5))/sqrt(sigma2);
	  pnorm_both(zval, &pgt, &plt, 2, 0);
	  val= mu==0.0 ? 1.0 : 2.0*MIN(pgt, plt);
	  if(type==6) {
	    val=ABSLOG(val);
	  } else if(type==7) {
	    val= pgt<=plt ? ABSLOG(val) : -ABSLOG(val);
	  } 
	} else {
	  error("Unrecognized val_type. Should not happen\n");
	}
      }
      resPtr[j+i*length(indlist)]=val;
    }

    destroyDRankList(list);
    matColPtr+=NROW(matrix);
  }
  
  UNPROTECT(1);
  return(res);
}
