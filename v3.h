#include<math.h>
#include<iostream>
using namespace std;
class v3 {

    public:

    v3(float _x =0,float _y=0,float _z=0);

    float distTo(const v3 &pt)const;
    float dotproduct(const v3 &v)const;
    float norm()const;
    bool operator!=(const v3 &pt) const;
    bool operator==(const v3 &pt) const;

float x;
float y;
float z; 

};