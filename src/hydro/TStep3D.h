#ifndef RN_TStep3D_h
#define RN_TStep3D_h

class TGrid3D ;

//! Base class for defining a stepping algorithm.
class TStep3D {

   public:

     virtual ~TStep3D() {; } 
     virtual int step(TGrid3D &gr, const double &dt) = 0;
     virtual double getDt(TGrid3D &gr) = 0;

} ;

#endif //define RN_TStep3D_h
