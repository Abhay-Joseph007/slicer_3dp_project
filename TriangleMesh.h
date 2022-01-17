#include<iostream>
using namespace std;

class TriangleMesh {

  public:
  
    TriangleMesh();

    size_t size() const ;

    void push_back(Triangle &t) ;

    v3 meshAABBSize() const ;

    const vector<Triangle>& getvTriangle() const;

    v3 getBottomLeftVertex() ;

    v3 getUpperRightVertex() ;

  public:

    int meshSize;
    v3 bottomLeftVertex;
    v3 upperRightVertex;
    vector<Triangle> vTriangle;
};