/*! \file stat_rank.h
    \brief Header for statistical ranking
    
    Details.
*/
#ifndef _STAT_RANK_H_
#define _STAT_RANK_H_

/*! \struct "*Dval"
 *  \struct DvalStruct
 *  \brief A class used to sort/rank double values
 *
 *  DRankStruct is synonym to *DRank
 */
typedef struct {
  int index; /*!< Input index (starting from 0) */
  double value; /*!< Value of the number */
  double rank; /*!< Fractional ranking (starting from 1) */
} DRankStruct, *DRank;

/*! \struct *DRankList
 *  \brief A dynamic array of DRanks, with length and unique length
 *
 *  DRankListStruct is synonym to DRankList
 */
typedef struct {
  DRank* list; /*!< Dynamic array of DRanks */
  int len; /*!< Length of the array */
  int ulen; /*!< Length of unique elements */
} DRankListStruct, *DRankList;

/*! \brief compare DRank by value */
int compareDRank(const void* a, const void* b);
/*! \brief compare DRank by index */
int compareDRankIndex(const void* a, const void* b);

/*! \brief initialize an DRankList 
 * \param array an array of double values
 * \param len the length of the array
 */
DRankList createDRankList(const double* array, int len);
/*! \brief print an DRankList */
void printDRankList(const DRankList list);
/*! \brief perform sorting and ranking at once */
void sortRankDRankList(DRankList list);
/*! \brief rank the DRankList and keeping the input order */
void rankDRankList(DRankList list);
/*! \brief sort the DRankList */
void sortDRankList(DRankList list);

/* private */
DRank createDRank(double, int);
void destroyDRank(DRank it);
int isRanked(const DRankList list);

#endif /* _STAT_RANK_H_ */
