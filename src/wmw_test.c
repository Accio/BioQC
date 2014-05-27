#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <math.h>

#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "stat_rank.h"
#include "wmw_test.h"

iArray iArrayCreate(int n) {
  iArray res=(iArray)malloc(sizeof(intArrayStruct));
  res->len=n;
  res->value=(int*)malloc(n*sizeof(int));
  return(res);
}
iArray iArrayCreateDef(int n, int val) {
  int i;
  iArray res=iArrayCreate(n);
  for(i=0;i<n;i++)
    res->value[i]=val;
  return(res);
}

void iArrayDestroy(iArray array) {
  array->len=0;
  free(array->value);
  free(array);
}
void iArrayPrint(const iArray array) {
  int i=0;
  for(i=0;i<array->len;++i) {
    printf("%d", array->value[i]);
    if(i!=array->len) {
      printf(" ");
    }
  }
  puts("");
}

dArray dArrayCreate(int n) {
  dArray res=(dArray)malloc(sizeof(doubleArrayStruct));
  res->len=n;
  res->value=(double*)malloc(n*sizeof(double));
  return(res);
}
dArray dArrayCopy(const double* array, int len) {
  int i;
  dArray res=dArrayCreate(len);
  for(i=0;i<len;i++) 
    *(res->value++)=*(array++); 
  return(res);
}

void dArrayDestroy(dArray array) {
  array->len=0;
  free(array->value);
  free(array);
}
void dArrayPrint(const dArray array) {
  int i=0;
  for(i=0;i<array->len;++i) {
    printf("%g", array->value[i]);
    if(i!=array->len) {
      printf(" ");
    }
  }
  puts("");
}

wmwRes wmwResCreate(double u, double ltp, double utp) {
  wmwRes res=(wmwRes)malloc(sizeof(wmwResStruct));
  res->U=u;
  res->ltP=ltp;
  res->gtP=utp;
  return(res);
}

void wmwResDestroy(wmwRes res) {
  free(res);
}

wmwRes wmwTest(const iArray index,
	       const dArray stat,
	       const double cor,
	       const double df) {
  int i,j;
  int n=stat->len;
  int n1=index->len;
  int n2=n-n1;
  double irsum=0; // sum of rank of index
  double U, mu, sigma2;
  double zlt, zut; // z lower/upper tail (normal distribution approximation)
  double plt, put;
  int ulen=0;
 
  DRankList list=createDRankList(stat->value, n);

  rankDRankList(list);

  // sum of the ranks of indexed elements
  for(i=0;i<n1;++i)
    irsum+=list->list[ index->value[i] ]->rank;
  ulen=list->ulen;
  
  U=n1*n2+n1*(n1+1.0)/2-irsum; 
  mu=(double)n1*n2/2;
  if(cor==0 || n1==1) {
    sigma2=n1*n2*(n+1.0)/12;
  } else {
    sigma2=(asin(1)*n1*n2+asin(0.5)*n1*n2*(n2-1.0)+\
	    asin(cor/2)*n1*(n1-1.0)*n2*(n2-1.0)+ \
	    asin((cor+1)/2)*n1*(n1-1)*n2)/2/M_PI;
  }
  
  if(ulen!=n) { // has ties
    int* tbl=(int*)malloc(ulen * sizeof(int));
    int ncount=0;
    double prod=0;
    sortDRankList(list);
    for(i=0;i<n;i=j+1) {
      j=i;
      while(j<n-1 && (list->list[j+1]->value==list->list[j]->value)) ++j;
      tbl[ncount++]=j-i+1;
    }
    for(i=0;i<ulen;++i)
      prod+=(0.0+tbl[i])/n*(tbl[i]+1)/(n+1)*(tbl[i]-1)/(n-1);
    sigma2*=(1-prod);
    free(tbl);
  }
  zlt=(U+0.5-mu)/sqrt(sigma2);
  zut=(U-0.5-mu)/sqrt(sigma2);
  
  plt = pt(zut, df, 0, 0);
  put = pt(zlt, df, 1, 0);

  wmwRes res=wmwResCreate(U, plt, put);
  return(res);
}
