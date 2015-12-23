/*! \file rank.c
  \brief statistical ranking

  Functions for statistical (fractional) ranking
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "stat_rank.h"

/*! \brief Create a DRank object
 * 
 * A item object holds a double value, its original index, and its rank.
 *
 * The rank is initialized with -1, and changed to a positive one (starting from 1) by sortRankDRankList. This is used to check whether that function has been run or not, so please do not change the initial value of rank.
 */

DRank createDRank(const double* ptr, int index) {
  DRank it=(DRank)malloc(sizeof(DRankStruct)); // TAKE CARE: use sizeof(DRank) here will produce bugs that are very difficult to debug
  it->index=index;
  it->vPtr=ptr;
  it->rank=-1.0;
  return(it);
}

/*! \brief destroy an item object */
void destroyDRank(DRank it) {
  it->vPtr=NULL;
  free(it);
}

/*! \brief clear an DRank object */
void clearDRank(DRank it) {
  it->vPtr=NULL;
  it->rank=-1.0;
  it->index=-1;
  free(it);
}

/*! \brief compare item objects by value */
int compareDRank (const void* a, const void* b) // dereference void pointer: *((T*)ptr)
{
  if(*(*(DRank*)a)->vPtr>*(*(DRank*)b)->vPtr) {
    return 1;
  } else if(*(*(DRank*)a)->vPtr==*(*(DRank*)b)->vPtr) {
    return 0;
  } else 
    return -1; // long debug needed: note it must return an integer!
}

/*! \brief compare item objects by input index */
int compareDRankIndex (const void* a, const void* b) // dereference void pointer: *((T*)ptr)
{
  return ((*(DRank*)a)->index-(*(DRank*)b)->index);

}

/*! \brief create an DRankList object 
 \param array: a double array
 \param len: the length of the double array
*/
DRankList createDRankList(const double* array, int len) {
  int i;
  DRankList res=(DRankList)malloc(sizeof(DRankListStruct));
  res->len=len;
  res->ulen=-1;
  res->list=(DRank*)malloc(len*sizeof(DRank));
  for(i=0;i<len;++i) {
    res->list[i]=createDRank(array+i, i);
  }
  return(res);
}

void destroyDRankList(DRankList list) {
  int i;
  for(i=0;i<list->len;i++)
    destroyDRank(list->list[i]);
  free(list->list);
  free(list);
}

void clearDRankList(DRankList list) {
  int i;
  for(i=0;i<list->len;i++)
    clearDRank(list->list[i]);
  list->ulen=-1;
  list->tieCoef=1.0;
}

/*! \brief print an DRankList object
*/
#ifdef DEBUG
void printDRankList(const DRankList list) {
  int i=0;
  printf("--DRankList (Len=%d, UniqLen=%d)--\n",
	 list->len, list->ulen);
  for(i=0;i<list->len;i++)
    printf("ilist[%d]=%.2f, index=%d, rank=%.1f\n",
	   i, 
	   *(list->list[i]->vPtr), 
	   list->list[i]->index,
	   list->list[i]->rank);
}
#endif

/*! \brief test whether the DRankList has been ranked 
 *
 * if sortRankDRankList has been run, the value will be 1, otherwise 0.
 */
int isRanked(const DRankList list) {return(list->list[0]->rank>0);}

/*! \brief: sort and gives rank to a DRankList
 *
 *  sortRankDRankList sorts and gives statistical (fractional) rank to a DRankList.
 * 
 *  sortRankDRankList runs once and only once (controlled by isRanked):
 *  Once the ranks have been set (i.e. ranks>0), the function will exit
 *  without doing anything.
 */ 
void sortRankDRankList(DRankList list) {
  if(isRanked(list)) return; // make sure that this function runs only once
  DRank* ll=list->list;
  int len=list->len;
  int i, j, k;
  int ucount=0;
  double *backup=(double*)malloc(len * sizeof(double));
  for(i=0;i<len; ++i)
    backup[i]=*(ll[i]->vPtr);

  qsort(ll, len, sizeof(DRank), compareDRank);

  for(i=0; i<len;i=j+1) {
    j=i;
    while((j<len-1) && backup[ll[j]->index] == backup[ll[j+1]->index]) {
      j++;
    }
    for(k=i;k<=j;k++)
      ll[k]->rank=(i+j+2)*0.5;
    ucount++;
  }

  free(backup);
  list->ulen=ucount;
}

/*! \brief: rankDRankList
 * \param list A DRankList object
 * It calls sortRankDRankList if the DRankList has not been ranked before
 * The items in the list are sorted by input index
 */
void rankDRankList(DRankList list) {
  sortRankDRankList(list);
  DRank* ll=list->list;
  int len=list->len;
  qsort(ll, len, sizeof(DRank), compareDRankIndex);
}

/*! \brief: sortDRankList
 * \param list A DRankList object
 * It calls sortRankDRankList if the DRankList has not been ranked before
 * The items in the list are sorted by ascending order of the values.
 */
void sortDRankList(DRankList list) {
  sortRankDRankList(list);
  DRank* ll=list->list;
  int len=list->len;
  qsort(ll, len, sizeof(DRank), compareDRank);
}

/*! \brief: prepareDRankList
 * \param list A DRankList object
 * It prepares a DRankList object to be used in Wilcoxon-Mann-Whitney tests
 */

int len(DRankList list) {return(list->len);}
int ulen(DRankList list) {return(list->ulen);}
double tieCoef(DRankList list) {return(list->tieCoef);}

void prepareDRankList(DRankList list) {
  rankDRankList(list);
  
  if(len(list)==ulen(list)) {
    list->tieCoef=1.0;
    return;
  } else {
    int n=len(list);
    int un=ulen(list);
    int *tbl=(int*)malloc(ulen(list) * sizeof(int));
    int ncount=0;
    double mTieCoef=0.0;
    int i, j;
    sortDRankList(list);
    for(i=0;i<n;i=j+1) {
      j=i;
      while(j<n-1 && (*(list->list[j+1]->vPtr)==*(list->list[j]->vPtr))) ++j;
      tbl[ncount++]=j-i+1;
    }
    for(i=0;i<un;++i)
      mTieCoef+=(0.0+tbl[i])/n*(tbl[i]+1)/(n+1)*(tbl[i]-1)/(n-1);
    list->tieCoef=1-mTieCoef;
    free(tbl);
    rankDRankList(list);
    return;
  }
}
