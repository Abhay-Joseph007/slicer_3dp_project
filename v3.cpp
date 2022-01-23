#include"v3.h"

v3::v3( float _x,float _y, float _z) 
{
 x=_x;
 y=_y;
 z=_z;
}
    
float v3 :: distTo (const v3 &pt)const{
        return sqrt ( pow(fabs(x-pt.x),2.0) + pow(fabs(y-pt.y), 2.0) + pow(fabs(z-pt.z),2.0 ));
    }
float v3:: dotproduct (const v3 &v)const{ 
      return (x*v.x + y*v.y + z*v.z); 
    }

float v3 :: norm()const{ 
      return sqrt(x*x+y*y+z*z); 
    }
  
 
  
