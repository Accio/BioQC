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

DRank createDRank(double value, int index) {
  DRank it=(DRank)malloc(sizeof(DRankStruct)); // TAKE CARE: use sizeof(DRank) here will produce bugs that are very difficult to debug
  it->index=index;
  it->value=value;
  it->rank=-1.0;
  return(it);
}

/*! \brief destroy an item object */
void destroyDRank(DRank it) {
  free(it);
}

/*! \brief compare item objects by value */
int compareDRank (const void* a, const void* b) // dereference void pointer: *((T*)ptr)
{
  return ((*(DRank*)a)->value-(*(DRank*)b)->value);
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
    res->list[i]=createDRank(array[i], i);
  }
  return(res);
}

void destroyDRankList(DRankList list) {
  int i;
  for(i=0;i<list->len;i++)
    destroyDRank(list->list[i]);
  free(list);
}

/*! \brief print an DRankList object
*/
void printDRankList(const DRankList list) {
  int i=0;
  printf("--DRankList (Len=%d, UniqLen=%d)--\n",
	 list->len, list->ulen);
  for(i=0;i<list->len;i++)
    printf("ilist[%d]=%.2f, index=%d, rank=%.1f\n",
	   i, 
	   list->list[i]->value, 
	   list->list[i]->index,
	   list->list[i]->rank);
}

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
    backup[i]=ll[i]->value;

  // qsort does not work properly!
  qsort(ll, len, sizeof(DRank), compareDRank);

  for(i=0; i<len;i=j+1) {
    j=i;
#ifdef DEBUG
    if(j<len-1)
      printf("[DEBUG] i=%d, j=%d, ind1=%d, ind2=%d, val1=%.1f, val2=%.1f\n", 
	     i, j, 
	     ll[j]->index, ll[j+1]->index, 
	     backup[ll[j]->index],
	     backup[ll[j+1]->index]);
#endif

    while((j<len-1) && backup[ll[j]->index] == backup[ll[j+1]->index]) {
      j++;
    }
    for(k=i;k<=j;k++)
      ll[k]->rank=(i+j+2)/2.;
    ucount++;
  }
  //for(i=0;i<len;++i) {
  //  (ll[i]->index)++; // index starts from 1
  //}
  free(backup);
  list->ulen=ucount;
}

/*! \brief: rankDRankList
 * \param list An DRankList
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
 * \param list An DRankList
 * It calls sortRankDRankList if the DRankList has not been ranked before
 * The items in the list are sorted by ascending order of the values.
 */
void sortDRankList(DRankList list) {
  sortRankDRankList(list);
  DRank* ll=list->list;
  int len=list->len;
  qsort(ll, len, sizeof(DRank), compareDRank);
}


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

inline void dArraySetValue(dArray array, int index, double value) {array->value[index]=value;}
inline double dArrayGetValue(const dArray array, int index) {return(array->value[index]);}

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
