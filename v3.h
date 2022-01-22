#include<math.h>
#include<iostream>
using namespace std;
class v3 {
    public:
    v3(float _x,float _y,float _z);
float dist(const v3 &pt);
float dotproduct(const v3 &v);
float norm();

float x;
float y;
float z; 
};