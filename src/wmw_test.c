#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <math.h>

#include "omp.h"
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

double wmw_test_core (const DRankList valList,
                    const int *inds, int nInds,
                    int nTotal, int type) {
    int nBg;
    int i;
    double U, mu, sigma2, zval, pgt, plt;
    double indRankSum; // sum of index rank
    double res;
    
    nBg = nTotal-nInds;
    indRankSum = 0.0;
    for(i = 0;i<nInds;++i)
        indRankSum += valList->list[inds[i]]->rank;
    
    U = nInds*nBg+nInds*(nInds+1.0)*0.5-indRankSum;
    if(type == 3) {
        res = U;
    } else {
        mu = (double)nInds*nBg*0.5; // NOT mu=n1*n2*0.5
        sigma2 = nInds*nBg*(nTotal+1.0)/12.0*tieCoef(valList); //NOT sigma2 = n1*n2*(n+1)/12*tieCoef
        
        if(type == 0 || type == 4) { /* greater */
            zval = (U+0.5-mu)/sqrt(sigma2); // z lower tail
            pnorm_both(zval, &pgt, &plt, 0, 0);
            res = type==0 ? pgt : ABSLOG(pgt);
        } else if (type == 1 || type == 5) { /* less */
            zval = (U-0.5-mu)/sqrt(sigma2); // z higher tail
            pnorm_both(zval, &pgt, &plt, 1, 0);
            res = type==1 ? plt : log10(plt);
        } else if (type == 2 || type == 6 || type == 7) { /* two sided*/
            zval = (U-mu-(U>mu ? 0.5 : -0.5))/sqrt(sigma2);
            pnorm_both(zval, &pgt, &plt, 2, 0);
            res = mu==0.0 ? 1.0 : 2.0*MIN(pgt, plt);
            if(type == 6) {
                res = ABSLOG(res);
            } else if(type == 7) {
                res = pgt<=plt ? ABSLOG(res) : -ABSLOG(res);
            }
        } else {
            error("Unrecognized type. Should not happen\n");
        }
    }
    return(res);
}

void wmw_test_list(const double *valPtr, int n,
                   SEXP indlist,
                   double *resPtr, int type) {
    DRankList list;
    int i;
    int n1;
    int* ip;

    list=createDRankList(valPtr, n);
    prepareDRankList(list);
    
#pragma omp parallel for
    for(i=0;i<length(indlist);++i) {
        ip=INTEGER(VECTOR_ELT(indlist,i));
        n1=length(VECTOR_ELT(indlist,i));
        resPtr[i]=wmw_test_core(list,
                                ip, n1,
                                n, type);
    }
    destroyDRankList(list);
}

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
    const int type=INTEGER(val_type)[0];
    const int m=length(indlist);
    const int n=NROW(matrix); // total number of samples AND nrow of matrix
    
    int i;
    double *matColPtr; // pointer to the current column of the matrix
    SEXP res;
    double *resPtr;
    
    res=PROTECT(allocMatrix(REALSXP, m, NCOL(matrix)));
    
    resPtr=REAL(res);
    matColPtr=REAL(matrix);
    
#pragma omp parallel for
    for(i=0; i<NCOL(matrix);++i) {
        wmw_test_list(matColPtr, n,
                      indlist,
                      resPtr, type);
        resPtr+=m;
        matColPtr+=n;
    }
    
    UNPROTECT(1);
    return(res);
}
