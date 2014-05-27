#ifndef _WMW_TEST_H_
#define _WMW_TEST_H_

typedef struct  {
  size_t len;
  int* value;
} intArrayStruct, *iArray;

iArray iArrayCreate(int n);
iArray iArrayCreateDef(int n, int val);
void iArrayDestroy(iArray array);
void iArrayPrint(const iArray array);

typedef struct {
  size_t len;
  double* value;
} doubleArrayStruct , *dArray;

dArray dArrayCreate(int n);
dArray dArrayCopy(const double *array, int len);
inline void dArraySetValue(dArray array, int index, double value) {array->value[index]=value;}
inline double dArrayGetValue(const dArray array, int index) {return(array->value[index]);}
void dArrayDestroy(dArray array);
void dArrayPrint(const dArray array);

typedef struct {
  double U; // U statistic
  double ltP; // lower tail P value
  double gtP; // higher tail P value
} wmwResStruct, *wmwRes;
wmwRes wmwResCreate(double u, double ltp, double utp);
void wmwResDestroy(wmwRes res);
inline double wmw_U(const wmwRes res) {return res->U;}
inline double wmw_ltP(const wmwRes res) {return res->ltP;}
inline double wmw_gtP(const wmwRes res) {return res->gtP;}
wmwRes wmwTest(const iArray index, const dArray stat, const double cor, const double df);

# endif /* _WMW_TEST_H_ */
