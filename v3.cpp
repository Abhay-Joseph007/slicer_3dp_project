#include"v3.h";


v3::v3( float _x=0,float _y=0, float _z=0) : x(_x),y(_y),z(_z){}

// v3::v3( float _x=0,float _y=0, float _z=0) 
// {
//  x=_x;
//  y=_y;
//  z=_z;
// }
    
float v3 :: dist (const v3 &pt){
        return sqrt ( pow(fabs(x-pt.x),2.0) + pow(fabs(y-pt.y), 2.0) + pow(fabs(z-pt.z),2.0 ));
    }
float v3:: dotproduct (const v3 &v) { 
      return (x*v.x + y*v.y + z*v.z); 
    }

float v3 :: norm() { 
      return sqrt(x*x+y*y+z*z); 
    }
  
 
  
