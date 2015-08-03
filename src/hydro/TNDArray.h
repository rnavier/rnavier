#ifndef  RN_TNDArray_h
#define  RN_TNDArray_h

#include "ndarray.h"

typedef ndarray::Array<double,1,1> TNDArray1D ;
typedef ndarray::Array<int,1,1>  TNDArray1I;

typedef ndarray::Array<double,2,2> TNDArray2D ;
typedef ndarray::Array<int,2,2>   TNDArray2I ;

typedef ndarray::Array<double,3,3> TNDArray3D ;
typedef ndarray::Array<int,3,3>  TNDArray3I ;

typedef ndarray::Array<double,4,4> TNDArray4D ;
typedef ndarray::Array<int,4,4>  TNDArray4I;

typedef ndarray::Array<double,5,5> TNDArray5D ;
typedef ndarray::Array<int,5,5>  TNDArray5I ;

typedef ndarray::Array<double,6,6> TNDArray6D ;
typedef ndarray::Array<int,6,6>  TNDArray6I ;

inline ndarray::Array<double,1,1> ndarray_alloc(const size_t &n1);
inline ndarray::Array<double,2,2> ndarray_alloc(const size_t &n1, const size_t &n2);
inline ndarray::Array<double,3,3> ndarray_alloc(const size_t &n1, const size_t &n2, const size_t &n3);

inline const ndarray::Array<double,1,1> ndarray_view(double *data, const size_t &n1);
inline const ndarray::Array<double,2,2> ndarray_view(double *data, const size_t &n1, const size_t &n2);
inline const ndarray::Array<double,3,3> ndarray_view(double *data, const size_t &n1, const size_t &n2, const size_t &n3);

#endif  

#ifndef RN_TNDArray_cxx
#define RN_TNDArray_cxx

ndarray::Array<double,1,1> ndarray_alloc(const size_t &n) 
{
   return ndarray::allocate(ndarray::makeVector(n)) ;
}
ndarray::Array<double,2,2> ndarray_alloc(const size_t &n1, const size_t &n2) 
{
   return ndarray::allocate(ndarray::makeVector(n1, n2)) ;
}
ndarray::Array<double,3,3> ndarray_alloc(const size_t &n1, const size_t &n2, const size_t &n3) 
{
   return ndarray::allocate(ndarray::makeVector(n1, n2, n3)) ;
}


const ndarray::Array<double,1,1> ndarray_view(double *data, const size_t &n) 
{
   return ndarray::external(data, ndarray::makeVector(n)) ;
}
const ndarray::Array<double,2,2> ndarray_view(double *data, const size_t &n1, const size_t &n2) 
{
   return ndarray::external(data, ndarray::makeVector(n1,n2)) ;
}
const ndarray::Array<double,3,3> ndarray_view(double *data, const size_t &n1, const size_t &n2, const size_t &n3) 
{
   return ndarray::external(data, ndarray::makeVector(n1,n2,n3)) ;
}
#endif

