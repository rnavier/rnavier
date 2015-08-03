#ifndef RN_TLimiter_h
#define RN_TLimiter_h

#include <iostream>
#include "gsl/gsl_math.h"
#include "ndarray.h"

class TLimiter
{
    private:
        int fKmethod;
        double minmod(const double &dpm, const double &dpp) const
        { double r = dpm*dpp/GSL_MAX(dpm*dpm, GSL_DBL_MIN);
            return dpm*GSL_MAX(0.0, GSL_MIN(1.0,r)); }
        double nolim(const double &dpm, const double &dpp) const
        { return (dpm + dpp)/2.;  } 

    public:
        enum EMethod: int { kNolim=0, kMinmod=1};
        TLimiter(): fKmethod(kMinmod){};
        TLimiter(int kmethod): fKmethod(kmethod){};
        double D(const ndarray::Array<double,1> &p, const int &i) const
        {
            switch (fKmethod){
                case EMethod::kNolim:
                    return nolim( p(i)-p(i-1), p(i+1)-p(i) );
                    break;
                case EMethod::kMinmod:
                    return minmod( p(i)-p(i-1), p(i+1)-p(i) );
                    break;
                default:
                    std::cout << "** TLimiter::D ** Bad method selctor in switch!" << std::endl;
                    return 0.;
            }
        }

};
#endif
