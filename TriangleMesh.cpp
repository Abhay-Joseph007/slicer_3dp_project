#include"TriangleMesh.h"
    
    TriangleMesh::TriangleMesh() : bottomLeftVertex(999999,999999,999999), upperRightVertex(-999999,-999999,-999999) { meshSize = 0;}

    size_t TriangleMesh::size() const {
      return meshSize;
    }

    void TriangleMesh::push_back(Triangle &t) {
      meshSize++;
      vTriangle.push_back(t);
      for (size_t i = 0; i < 3; ++i) {
        if (t.v[i].x < bottomLeftVertex.x) { bottomLeftVertex.x = t.v[i].x; }
        if (t.v[i].y < bottomLeftVertex.y) { bottomLeftVertex.y = t.v[i].y; }
        if (t.v[i].z < bottomLeftVertex.z) { bottomLeftVertex.z = t.v[i].z; }
        if (t.v[i].x > upperRightVertex.x) { upperRightVertex.x = t.v[i].x; }
        if (t.v[i].y > upperRightVertex.y) { upperRightVertex.y = t.v[i].y; }
        if (t.v[i].z > upperRightVertex.z) { upperRightVertex.z = t.v[i].z; }
      }
    }

    v3 TriangleMesh::meshAABBSize() const {
        return v3 ( upperRightVertex.x - bottomLeftVertex.x, 
                    upperRightVertex.y - bottomLeftVertex.y, 
                    upperRightVertex.z - bottomLeftVertex.z );
    }

    const vector<Triangle>& TriangleMesh::getvTriangle() const { return vTriangle; }

    v3 TriangleMesh::getBottomLeftVertex() const { return bottomLeftVertex; }

    v3 TriangleMesh::getUpperRightVertex() const { return upperRightVertex; }

