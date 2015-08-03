#include "TNDArray.h"
#include "ndarray.h"
#include <typeinfo>
#include <iostream>

using namespace std ;

using namespace ndarray ;

//// Acceptable function argumentes 
//double add_row(ArrayRef<double, 1,1> &row)
//{
//   double v = 0. ;
//   for (int i=0;  i <3 ; i++) {
//      v += row[i] ;
//      row[i] += 2. ;
//   }
//   return v ;
//}

// Preferred  function arguments if you have no intention of manipulating the
// pointers. Just accesing 
double add_row_v2(const Array<double, 1,1> &row)
{
   double v = 0. ;
   for (int i=0;  i <3 ; i++) {
      v += row[i] ;
      row[i] += 3. ;
   }
   return v ;
}


int main(int argc, char **argv) 
{
   cout << "=> Allocating with standard sequence <=" << endl;
   Vector<int, 3> dims = makeVector(2,3, 3) ;
   TNDArray3D y = allocate(dims);
   y.deep() = 2. ;
   cout << y << endl;

   cout << "=> Allocating with abbrieviated sequence 3d <=" << endl;
   TNDArray3D y2 = ndarray_alloc(2,3,3);
   y2.deep() = y ;
   y2(0,0,0) = 3 ;  y(0,1,2) = 3 ; y2(1, 0, 0) = 3 ; y2(1,2,1) = 3 ; 
   cout << y2 << endl;
   cout << "... but we haven't changed y " << endl;
   cout << y << endl;

   cout << "=> Ndarray allocation with abrrievated sequence 1d <=" << endl;
   TNDArray1D y1d = ndarray_alloc(3) ;
   y1d.deep() = 4. ;
   for (int k = 0 ; k < 3 ; k++)  
       cout << y1d[k] << endl;

   cout << "=> ndarray understanding when we need deep <=" << endl;
   cout << y[view(1)(2)()] << endl;
   cout << y2[view(1)(2)()] << endl;
   y2[view(1)(2)()] =  y[view(1)(2)()] ;
   cout << y[view(1)(2)()] << endl;
   cout << y2[view(1)(2)()] << endl;
   y2[view()] = y ;
   cout << y2 << endl;
   cout << y << endl;



   cout << "=> Ndarray view <=" << endl;
   double data[3]  = {2, 4, 6};
   auto data_view = ndarray_view(data, 3) ;
   for (int k = 0 ; k < 3 ; k++)  
       cout << data_view[k] << endl;


   cout << "=> Ndarray passing <=" << endl;
   TNDArray1D yr1 ;
   for (int j = 0 ; j < 3 ; j++) {
      cout << "=>Accessing the " << j << " row" << endl;
      yr1 = y[view(1)(j)()] ;
      const auto yr2 = y[view(1)(j)() ] ;
      const TNDArray1D yr3 = y[view(1)(j)() ] ;

      //cout << add_row(yr)  << endl;
      cout << add_row_v2(yr1)  << endl;
      cout << add_row_v2(y[view(1)(j)()])  << endl;
      cout << add_row_v2(yr2)  << endl;
      cout << add_row_v2(yr3)  << endl;
   }

}
