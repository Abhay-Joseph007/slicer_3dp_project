#include"TriangleClass.h"
#include"v3.h"

    Triangle::Triangle(v3 n, v3 v0, v3 v1, v3 v2) : normal(n) {
      v[0] = v0; 
      v[1] = v1; 
      v[2] = v2;
      zMin = +99999999.9; 
      zMax = -99999999.9;
      setZMin(v0.z); setZMin(v1.z); setZMin(v2.z);
      setZMax(v0.z); setZMax(v1.z); setZMax(v2.z);
    }

    void Triangle::setZMin (float z) { 
      if (z < zMin) {
        zMin = z; 
      }
    }

    void Triangle::setZMax (float z) { 
      if (z > zMax) {
        zMax = z; 
      }
    }

    // Triangle& Triangle::operator-=(const v3 &pt) { 
    //   v[0] -= pt; 
    //   v[1] -= pt;  
    //   v[2] -= pt; 
    //   return *this;
    // }

    // bool Triangle::operator<(const Triangle &t) { 
    //    return zMin < t.zMin; 
    // }

    // //this doesnt need to be defined in teh class namespace since it is a friend of the class
    // ostream& operator<<(ostream& os, const Triangle& t) {
    //   os << "V0: (" << t.v[0] << "); V1: (" << t.v[1] << "); V2: (" << t.v[2] << ")";
    //   return os;
    // } 