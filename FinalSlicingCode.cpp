#include<iostream>
#include<fstream>
#include<math.h>

#define MAX_READ_CHUNK 80
#define FACET_READ_SIZE 12

// rounding a vertex (x,y,z) to an even multiple of eps
// v3 v3_round(float x,float y,float z,double eps){
// }

//returns a triangle given 3 vertices and normal after rounding off easch vertice to even multiple of eps
Triangle make_triangle(float v12[FACET_READ_SIZE],double eps){
    return Triangle (v3(v12[0],v12[1],v12[2]), v3_round(v12[3],v12[4],v12[5], eps), v3_round(v12[6],v12[7],v12[8], eps), v3_round(v12[9],v12[10],v12[11], eps));
}

// checks if any vertices of a triangle are coinciding
// bool degenerate (Triangle t) {
// }

void stlToMeshInMemory (const char *stlFile, TriangleMesh *mesh, double eps){

  int no_of_triangles = 0;                        //to store no of triangles
  char var[MAX_READ_CHUNK];                       //to read and store header from input file
  float v12[FACET_READ_SIZE];                    //to store the vector read from input file
  unsigned short attribute;                     //to store the attribute
  int ndegenerated = 0;

  ifstream FileRead(stlFile, ios::binary);

  FileRead.read(var,MAX_READ_CHUNK);
  FileRead.read((char*)no_of_triangles,4);

  for(int i = 0;i < no_of_triangles;i++){
    FileRead.read((char*)v12,sizeof(float)*FACET_READ_SIZE);
    FileRead.read((char*)attribute,2);

    Triangle t = make_triangle(v12,eps);

    if (!degenerate(t)) {
       mesh->push_back (t);
    }
    else {
      ndegenerated++;
    }
  }

  FileRead.close();
  cout << "number of degenerated triangles = " << ndegenerated << endl;
}

int main(int argc, char** argv)
{
  char *model;
  model = argv[argc-1];

  double eps = 0.004

  TriangleMesh mesh;

  stlToMeshInMemory (model, &mesh, eps);

  return 0;
}
