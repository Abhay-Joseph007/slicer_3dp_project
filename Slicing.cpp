#include<iostream>
#include<fstream>
#include<math.h>
#include<string.h>
#include<assert.h>
#include<unordered_map>
#include<vector>
#include<algorithm>

using namespace std;

#define MAX_READ_CHUNK 80
#define FACET_READ_SIZE 12

#define DEG_TO_RAD(x) (x*0.0174532925199f)
#define FILE_STL_BIN 0
#define FILE_STL_ASCII 1
#define FILE_AMF_ASCII 2

class v3 {

  public:
    
    v3 (float _x=0, float _y=0, float _z=0) : x(_x), y(_y), z(_z) {}

    float distTo (const v3 &pt) const { 
      return sqrt ( pow(fabs(x-pt.x), 2.0) + pow(fabs(y-pt.y), 2.0) + pow(fabs(z-pt.z), 2.0) ); 
    }

    array<float,3> getCoords() { 
      array<float,3> c = {{x, y, z}}; 
      return c; 
    }

    float dotproduct (const v3 &v) const { 
      return (x*v.x + y*v.y + z*v.z); 
    }

    // void transform (const glm::mat4 &mat) { 
    //   glm::vec4 v = glm::vec4(x, y, z, 1.0);
    //   glm::vec4 vt = mat*v; 
    //   x = (vt.x); y = (vt.y); z = (vt.z);
    // }
    
    v3& operator-=(const v3 &pt) { 
      x = (x-pt.x); 
      y = (y-pt.y); 
      z = (z-pt.z); 
      return *this;    
    }

    v3 operator-(const v3 &pt) { 
      return v3 ((x-pt.x), (y-pt.y), (z-pt.z)); 
    }
    
    v3 operator+(const v3 &pt) { 
      return v3 ((x+pt.x), (y+pt.y), (z+pt.z)); 
    }
    
    v3 operator/(float a) { 
      return v3 ((x/a), (y/a), (z/a)); 
    }
    
    v3 operator*(float a) { 
      return v3 ((x*a), (y*a), (z*a)); 
    }
    
    bool operator<(const v3 &pt) const { 
      return z < pt.z; 
    }
    
    bool operator>(const v3 &pt) const { 
      return z > pt.z; 
    }
    
    bool operator==(const v3 &pt) const {
      return distTo(pt) < 0.005; 
    }
    
    bool operator!=(const v3 &pt) const { 
      return distTo(pt) > 0.005; 
    }
    
    float normalize() const { 
      return sqrt(x*x+y*y+z*z); 
    }

    // string getLabel() const {
    //   stringstream ss;
    //   ss << x << "|" << y << "|" << z;
    //   return ss.str();
    // }
    
    // friend ostream& operator<<(ostream& os, const v3& v) {
    //   os << "x: " << v.x << "; y: " << v.y << "; z: " << v.z;
    //   return os;
    // }

  public:

    float x;
    float y;
    float z;
};

v3 operator-(const v3 &a, const v3 &b) {return v3((a.x-b.x), (a.y-b.y), (a.z-b.z)); }

v3 operator+(const v3 &a, const v3 &b) {return v3((a.x+b.x), (a.y+b.y), (a.z+b.z)); }

/**************************************************************************
 *                        CLASS  LineSegment                              * 
 **************************************************************************/
class LineSegment {

  public: 

    LineSegment (v3 p0=v3(), v3 p1=v3(), int i=0) { 
      v[0] = p0; 
      v[1] = p1;
      index = i;
      vertical = false;  
      if ((v[1].x - v[0].x) != 0) {
        a = (v[1].y - v[0].y)/(v[1].x - v[0].x);
        b = (v[0].y - (a * v[0].x));
      }
      else {
        vertical = true;
      }
    }

    bool operator==(const LineSegment &ls) const { 
      return ((v[0] == ls.v[0]) && (v[1] == ls.v[1])); 
    }

    // friend ostream& operator<<(ostream& os, const LineSegment& ls) {
    //   os << "V0: (" << ls.v[0] << "); V1: (" << ls.v[1] << ")";
    //   return os;
    // }

  public:

    v3 v[2];
    double a;
    double b;
    bool vertical;
    int index;
}; 

class Triangle    {

  public:
  
    Triangle(v3 n, v3 v0, v3 v1, v3 v2) : normal(n) {
      v[0] = v0; 
      v[1] = v1; 
      v[2] = v2;
      zMin = +99999999.9; 
      zMax = -99999999.9;
      setZMin(v0.z); setZMin(v1.z); setZMin(v2.z);
      setZMax(v0.z); setZMax(v1.z); setZMax(v2.z);
    }

    void setZMin (float z) { 
      if (z < zMin) {
        zMin = z; 
      }
    }

    void setZMax (float z) { 
      if (z > zMax) {
        zMax = z; 
      }
    }

    Triangle& operator-=(const v3 &pt) { 
      v[0] -= pt; 
      v[1] -= pt;  
      v[2] -= pt; 
      return *this;
    }

    bool operator<(const Triangle &t) { 
       return zMin < t.zMin; 
    }

    // friend ostream& operator<<(ostream& os, const Triangle& t) {
    //   os << "V0: (" << t.v[0] << "); V1: (" << t.v[1] << "); V2: (" << t.v[2] << ")";
    //   return os;
    // }
  
  public:

    v3 v[3];
    v3 normal;
    float zMin;
    float zMax;
};  

/**************************************************************************
 *                      CLASS  TriangleMesh                               * 
 **************************************************************************/
class TriangleMesh {

  public:
  
    TriangleMesh() : bottomLeftVertex(999999,999999,999999), upperRightVertex(-999999,-999999,-999999) { meshSize = 0;}

    size_t size() const {
      return meshSize;
    }

    void push_back(Triangle &t) {
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

    v3 meshAABBSize() const {
        return v3 ( upperRightVertex.x - bottomLeftVertex.x, 
                    upperRightVertex.y - bottomLeftVertex.y, 
                    upperRightVertex.z - bottomLeftVertex.z );
    }

    const vector<Triangle>& getvTriangle() const { return vTriangle; }

    v3 getBottomLeftVertex() const { return bottomLeftVertex; }

    v3 getUpperRightVertex() const { return upperRightVertex; }

  public:

    int meshSize;
    v3 bottomLeftVertex;
    v3 upperRightVertex;
    vector<Triangle> vTriangle;
};

template<typename T> inline void hash_combine(size_t &seed, const T &v) {
  hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct HashV3 {
  size_t operator() (const v3 &v) const {
    size_t h = hash<float>()(v.x);
    hash_combine (h, v.y);
    hash_combine (h, v.z);
    return h;
  }
};

typedef unordered_map<v3, vector<v3>, HashV3> PointMesh;

/*Structures*/
typedef struct node {
  Triangle t;
  struct node *next;
  struct node *prev;
} Mesh_Triangle_Node_t;

typedef struct _list {
   Mesh_Triangle_Node_t *head;
   Mesh_Triangle_Node_t *tail;
} Mesh_Triangle_List_t;

/*-----------------------------------------------------------------------*/
Mesh_Triangle_List_t* Mesh_Triangle_List_create (void) {
  Mesh_Triangle_List_t *L = (Mesh_Triangle_List_t *)malloc(sizeof(Mesh_Triangle_List_t));
  L->head = NULL;
  L->tail = NULL;
  return L;
}

/*-----------------------------------------------------------------------*/
void Mesh_Triangle_List_insert (Triangle t, Mesh_Triangle_List_t *L) {
  Mesh_Triangle_Node_t *node = (Mesh_Triangle_Node_t *)malloc(sizeof(Mesh_Triangle_Node_t));
  node->t = t;
  node->next = L->head;
  node->prev = NULL;
  if (L->head == NULL) {
    /*New head*/
    L->head = L->tail = node;
  }
  else {
    L->head->prev = node;
    L->head = node;
  }
}

/*-----------------------------------------------------------------------*/
void Mesh_Triangle_List_union (Mesh_Triangle_List_t *L1, Mesh_Triangle_List_t *L2) {
  if ( (L1->head != NULL) && (L2->head != NULL) ) {
    L1->tail->next = L2->head;
    L2->head->prev = L1->tail;
    L1->tail = L2->tail;;
  }
  else if (L2->head != NULL) {
    L1->head = L2->head;
    L1->tail = L2->tail;
  }
}

/*-----------------------------------------------------------------------*/
Mesh_Triangle_Node_t* Mesh_Triangle_List_remove (Mesh_Triangle_List_t *L, Mesh_Triangle_Node_t *node) {
  if ((node->prev == NULL) && (node->next == NULL)) {
    free (node);
    L->head = NULL;
    L->tail = NULL;
    return NULL;
  }
  else if (node->prev == NULL) {
    node->next->prev = NULL;
    L->head = node->next;
    free (node);
    return L->head;
  }
  else if (node->next == NULL) {
    node->prev->next = NULL;
    L->tail = node->prev;
    free (node);
    return NULL;
  }
  else {
    Mesh_Triangle_Node_t *next = node->next;
    node->next->prev = node->prev;
    node->prev->next = next;
    free (node);
    return next;
  }
}

v3 R3_Mesh_Side_slice (v3 vi, v3 vj, float Z) {
   double dx = vj.x - vi.x;
   double dy = vj.y - vi.y;
   double dz = vj.z - vi.z;
   assert(dz != 0);
   double frac = (Z - vi.z)/dz;
   float xint = (float)(frac*dx + (double)vi.x);
   float yint = (float)(frac*dy + (double)vi.y);
   return (v3){ xint, yint, Z };
}


typedef struct _contour {
   bool external;
   bool clockwise;
   vector<v3> P;
} contour;

//Rounds x to an integer multiple of eps
float xround (float x, double eps, int mod, int rem) {
  double y = round((double)x/(mod * eps));
  double z = (y * mod + rem) * eps;
  return (float)z;
}

// rounding a vertex (x,y,z) to an even multiple of eps
v3 v3_round(float x,float y,float z,double eps){
  v3 p;
  p.x = xround(x, eps, 2, 0);
  p.y = xround(y, eps, 2, 0);
  p.z = xround(z, eps, 2, 0);
  return p;
}

//returns a triangle given 3 vertices and normal after rounding off easch vertice to even multiple of eps
Triangle make_triangle ( 
   float n0, float n1, float n2,
   float f0, float f1, float f2,
   float f3, float f4, float f5,
   float f6, float f7, float f8,
   const char *rotate, double eps
 ) {

    if (strcmp(rotate,"true") == 0) {
      return Triangle (v3(n0, n1, n2), v3_round(f0, f2, f1, eps), v3_round(f3, f5, f4, eps), v3_round(f6, f8, f7, eps));
    }
    else {
      return Triangle (v3(n0, n1, n2), v3_round(f0, f1, f2, eps), v3_round(f3, f4, f5, eps), v3_round(f6, f7, f8, eps));
    }
}

// checks if any vertices of a triangle are coinciding
bool degenerate (Triangle t) {
  if (t.v[0].distTo(t.v[1]) < 0.000001) { return true; }
  if (t.v[1].distTo(t.v[2]) < 0.000001) { return true; }
  if (t.v[2].distTo(t.v[0]) < 0.000001) { return true; }
  return false;
}

LineSegment R3_Mesh_Triangle_slice (Mesh_Triangle_Node_t *t, float Z) {
   assert((t->t.zMin < Z) && (t->t.zMax > Z));
   int np = 0; /* Number of segment endpoints found */
   LineSegment seg;
   for (int i = 0; i < 3; i++) {
      /* Get side {i} of triangle: */
      int j = (i == 2 ? 0 : i+1);
      v3 vi = (t->t.v[i]);
      v3 vj = (t->t.v[j]);
      /* Check for intersection of plane with {vi--vj}. */
      /* Must consider segment closed at bottom and open at top in case {Z} goes through a vertex. */
      float vzMin = (vi.z < vj.z ? vi.z : vj.z);
      float vzMax = (vi.z > vj.z ? vi.z : vj.z);
      if ((vzMin <= Z) && (vzMax > Z)) {
         v3 p = R3_Mesh_Side_slice (vi, vj, Z);
         assert(np < 2);
         seg.v[np] = p;
         np++;
      }
   }
   assert(np == 2);
   return seg; 
}

/*-----------------------------------------------------------------------*/
/*Gets an arbitrary segment from {H}, removes it from {H} and returns it as a trivial chain. */
vector<v3> IncrementalStartLoop(vector<PointMesh> &H) {
   vector<v3> P;
   auto it = (H[0]).begin();
   v3 u = (*it).first;
   vector<v3> vw = (*it).second;
   v3 v = vw.at(0);
   P.push_back(u);
   P.push_back(v);
   (H[0][u]).erase(std::remove((H[0][u]).begin(), (H[0][u]).end(), v), (H[0][u]).end());
   if (H[0][u].size() == 0) { (H[0]).erase(u); }
   (H[0][v]).erase(std::remove((H[0][v]).begin(), (H[0][v]).end(), u), (H[0][v]).end());
   if (H[0][v].size() == 0) { (H[0]).erase(v); }
   return P;
}

/*-----------------------------------------------------------------------*/
/*Extends the chain {P} wih segments from {H}, removing them, while possible. */
void IncrementalExtendLoop(vector<v3> &P, vector<PointMesh> &H) { 
  int index = 0;
  int n = P.size();
  v3 first = P.front();  
  v3 current = P.back(); 
  v3 last;
         
  /* Collect other vertices: */
  while (true) {
    auto it = (H[0]).find(current);
    if (it == (H[0]).end()) { /*Vertex {current} is a dead end:*/ break; }
    v3 key1 = (*it).first; assert(key1 == current);  /*Paranoia check.*/
            
    /*Get {next}, the first unused neighbor of {current}:*/
    vector<v3> vw = (*it).second; /*Unused neighbors of {current}.*/
    assert (vw.size() != 0); 
    v3 next = vw.at(0); /*First unused neighbor of {current}.*/

    /*Append the segment {(current,next)} to {P} and delete from {H}:*/
    P.push_back(next);

    /*Remove the segment {(current,next)} from {H}:*/
    (H[0][current]).erase(std::remove((H[0][current]).begin(), (H[0][current]).end(), next), (H[0][current]).end());
    if (H[0][current].size() == 0) { (H[0]).erase(current); } 
    (H[0][next]).erase(std::remove((H[0][next]).begin(), (H[0][next]).end(), current), (H[0][next]).end());
    if (H[0][next].size() == 0) { (H[0]).erase(next); } 

    if (next == first) { /*We closed a loop:*/ break; }

    /*Move on:*/
    current = next;
  }
}

/*Reverses the chain {P}.*/
void IncrementalReverseLoop(vector<v3> &P) { 
  std::reverse(P.begin(),P.end());
}

void ContourConstruction (vector<LineSegment> segs, vector<contour> polygons[], int plane) {

   bool verbose = false;

   //clock_t contour_begin = clock();

   /*Creating the hash table.*/
   vector<PointMesh> H(1);

   /*Rounding vertices and filling the hash table.*/
   double eps = 1/128.0;
   for (std::vector<LineSegment>::iterator i = segs.begin(); i != segs.end(); i++) {
      LineSegment q = *i;
      q.v[0].x = round(q.v[0].x / eps) * eps;
      q.v[0].y = round(q.v[0].y / eps) * eps;
      q.v[0].z = plane;
      q.v[1].x = round(q.v[1].x / eps) * eps;
      q.v[1].y = round(q.v[1].y / eps) * eps;
      q.v[1].z = plane;
      if (q.v[0].distTo(q.v[1]) > 0.0001) {
         (H[0][q.v[0]]).push_back(q.v[1]);
         (H[0][q.v[1]]).push_back(q.v[0]);
      }
   }

   /* Count vertices by degree: */
   if (verbose) {
     int degmax = 10;
     int ctdeg[degmax+1];
     for (int deg = 0; deg <= degmax; deg++) { ctdeg[deg] = 0; }
     for (auto i = (H[0]).begin(); i != (H[0]).end(); i++) {
        vector<v3> L = (*i).second;
        int deg = L.size();
        if (deg > degmax) { deg = degmax; }
        ctdeg[deg]++;
     }
     assert(ctdeg[0] == 0);
     bool closedSlice = true;
     for (int deg = 1; deg <= degmax; deg++) { 
       if (((deg % 2) != 0) && (ctdeg[deg] > 0)) { closedSlice = false; }
       if ((verbose || (deg != 2)) && (ctdeg[deg] != 0))
         { cout << "there are " << ctdeg[deg] << " vertices of degree " << deg << " on plane " << plane << endl; }
     }
     if (!closedSlice) { cout << "** contours of plane " << plane << " are not closed" << endl; }
   }

   /*Contour construction.*/
   bool maximal = true;
   while (!(H[0]).empty()) {
     if (maximal) {
       vector<v3> P = IncrementalStartLoop(H);
       IncrementalExtendLoop(P,H);
       if (P.front() != P.back()) { //Chain {P} is open
         IncrementalReverseLoop(P);
         IncrementalExtendLoop(P,H);
       }
       polygons[plane].push_back({false, false, P});
     }
     else {
       vector<v3> P = IncrementalStartLoop(H);
       IncrementalExtendLoop(P,H);
       polygons[plane].push_back({false, false, P});
     }
   }
   //clock_t contour_end = clock();
   //loopclosure_time += double(contour_end - contour_begin)/CLOCKS_PER_SEC;
}

Mesh_Triangle_List_t** IncrementalSlicing_buildLists (bool srt, double delta, const TriangleMesh *mesh, vector<float> P) {

  int k = P.size(); /* Number of planes. */

  Mesh_Triangle_List_t **L = (Mesh_Triangle_List_t **)malloc((k+1) * sizeof(Mesh_Triangle_List_t *));

  for (size_t p = 0; p <= k; p++) { L[p] = Mesh_Triangle_List_create(); }

  const vector<Triangle> &T = mesh->getvTriangle();

  int n = T.size(); /* Number of triangles. */

    /* Uniform slicing - compute list index: */
    for (auto it = T.begin(), itEnd = T.end(); it != itEnd; ++it) {
      Triangle t = *it;
      int p;
      if (t.zMin < P[0]) {
        p = 0;
      }
      else if (t.zMin > P[k-1]) {
        p = k;
      }
      else {
        p = floor((t.zMin - P[0])/delta) + 1;
      }
      Mesh_Triangle_List_insert (t, L[p]);
    }

      return L;
}
  
void IncrementalSlicing (const TriangleMesh *mesh, vector<float> P, float delta, bool srt, vector<contour> polygons[], bool chaining, bool orienting) {

  /*Slicing*/

  int k = P.size();

  vector<LineSegment> segs[k];

  /* Classify triangles by the plane gaps that contain their {zMin}: */
  Mesh_Triangle_List_t **L = IncrementalSlicing_buildLists (srt, delta, mesh, P);
  /* Now perform a plane sweep from bottom to top: */

  Mesh_Triangle_List_t *A = Mesh_Triangle_List_create(); /* Active triangle list. */
  for (int p = 0; p < k; p++) {
    /* Add triangles that start between {P[p-1]} and {P[p]}: */
    Mesh_Triangle_List_union (A, L[p]);
    /* Scan the active triangles: */
    Mesh_Triangle_Node_t *aux = A->head;
    while (aux != NULL) {
      Mesh_Triangle_Node_t *next = aux->next;
      if (aux->t.zMax < P[p]) {
        /* Triangle is exhausted: */
        Mesh_Triangle_List_remove (A, aux);
      } else {
        /* Compute intersection: */
        if ((aux->t.zMin < P[p]) && (aux->t.zMax > P[p])) {
          LineSegment seg = R3_Mesh_Triangle_slice (aux, P[p]);
          segs[p].push_back(seg);
          //intersections++;
        }
      }
      aux = next;
    }
  }
  free(L);
  /*End-Slicing*/

    if (chaining) {
    /*Contour construction:*/
    for (size_t p = 0; p < k; p++) {
      if (!segs[p].empty()) {
         ContourConstruction (segs[p], polygons, p);
         #ifdef DEBUG
           char fname[256];
           sprintf(fname, "slice_%03d.txt", (int)p);
           FILE *fslice = fopen(fname, "w");
           fprintf(fslice, "----------- Segmentos ---------------\n");
           for (int ii = 0; ii < segs[p].size(); ii++) {
              LineSegment ss = segs[p].at(ii); 
              fprintf(fslice, "%f %f   %f %f\n", ss.v[0].x, ss.v[0].y, ss.v[1].x, ss.v[1].y);           
           }
           fprintf(fslice, "---------- End Segmentos --------------\n");
           record_polygons (polygons[p], fslice);
           fclose(fslice);
         #endif
        //  if (orienting) {   
        //    ray_casting (polygons[p]);
        //  }
         segs[p].clear(); 
      }
    }
    /*End construction.*/
  }
  else {
    // export_svg_no_chaining ("segments.svg", segs, k, mesh->meshAABBSize());
  }
}

/*Compute uniform and adaptive z-plane coordinates!*/
vector<float> compute_planes (const TriangleMesh *mesh, float max_thickness, char *adaptive, double eps, float *delta) {

  bool rounding = true; /*To avoid that the Z-coordinates of all planes are distinct from the Z-coordinates of all vertices.*/

  /* Vector to keep the plane coordinates: */
  vector<float> Planes;

  /* Assuming the model as a 3D axis-aligned bounding-box: */
  double model_zmax = std::max(mesh->getUpperRightVertex().z, mesh->meshAABBSize().z);

  double model_zmin = mesh->getBottomLeftVertex().z;

  if (strcmp(adaptive, "false") == 0)  { /*Uniform slicing: */

    double spacing = (rounding ? xround (max_thickness, eps, 2, 0) : max_thickness); /*Plane spacing even multiple of {eps}*/

    double P0 = xround (model_zmin - spacing, eps, 2, 1); /*First plane odd multiple of {eps}.*/

    int no_planes = 1 + (int)((model_zmax + spacing - P0)/spacing); /* Number of planes: */

    cout << "eps = " << eps << endl;
    cout << "max thickness = " << max_thickness << endl;
    cout << "rounded plane spacing spacing = " << spacing << endl;
    cout << "model zmin = " << model_zmin << ", model zmax = " << model_zmax << ", first plane Z = " << P0 << ", number of planes = " << no_planes << endl;

    for (size_t i = 0; i < no_planes; i++) {
      /* Building the vector with the slice z coordinates: */
      float Pi = (float)(P0 + i * spacing);
      if ((Pi > model_zmin) && (Pi < model_zmax)) {
          Planes.push_back ((float)(P0 + i * spacing));
      }
    }
    *delta = (float)(spacing);
  }
  else { /*Adaptive slicing z-planes: */

    float zplane = 0.0;
    float min_thickness = 0.016;
    Planes.push_back (model_zmin + zplane);

    while ((model_zmin + zplane) <= model_zmax) {
      double vrandom = min_thickness + (max_thickness - min_thickness) * (rand() / (double)RAND_MAX);
      double coordinate = xround (model_zmin + zplane + vrandom, eps, 2, 1);
      if (coordinate >= model_zmax) { break; }
      Planes.push_back (coordinate);
      zplane += vrandom;
    }
  }
  return Planes;
}


//* Read the given STL file name (ascii or binary is set using ‘isBinaryFormat’)
//and generate a Triangle Mesh object in output parameter ‘mesh’
int stlToMeshInMemory (const char *stlFile, TriangleMesh *mesh, bool isBinaryFormat, const char *rotate, double eps) {

  int ndegenerated = 0;

  if (!isBinaryFormat) {

    ifstream in(stlFile);

    if (!in.good()) {
      return 1;
    }

    std::string s0, s1;

    float n0, n1, n2;
    float f0, f1, f2;
    float f3, f4, f5;
    float f6, f7, f8;

    while (!in.eof()) {
      in >> s0; /*s0 can be facet or solid!*/
      if (s0 == "facet") {
        in >> s1 >> n0 >> n1 >> n2; /* normal x y z. */
        in >> s0 >> s1;             /* loop. */
        in >> s0 >> f0 >> f1 >> f2; /* vertex x y z. */
        in >> s0 >> f3 >> f4 >> f5; /* vertex x y z. */
        in >> s0 >> f6 >> f7 >> f8; /* vertex x y z. */
        in >> s0;                   /* endloop. */
        in >> s0;                   /* endfacet.*/
        Triangle t = make_triangle (n0, n1, n2, f0, f2, f1, f3, f5, f4, f6, f8, f7, rotate, eps);
        if (!degenerate(t)) {
           mesh->push_back (t);
        }
        else {
          ndegenerated++;
        }
      }
      else if (s0 == "endsolid") {
         break;
      }
    }
    in.close();
  }
  else {
    FILE *f = fopen (stlFile, "rb");
    if (!f) {
      return 1;
    }
    char title[80];
    int nFaces;
    int err;
    err = fread (title, 80, 1, f);
    err = fread ((void*)&nFaces, 4, 1, f);
    float v[12]; /* normal = 3, vertices = 9 (12) */
    unsigned short uint16;
    /* Every Face is 50 Bytes: Normal(3*float), Vertices(9*float), 2 Bytes Spacer */
    for (size_t i=0; i<nFaces; ++i) {
      for (size_t j=0; j<12; ++j) {
        err = fread((void*)&v[j], sizeof(float), 1, f);
      }
      err = fread((void*)&uint16, sizeof(unsigned short), 1, f); // spacer between successive faces
      Triangle t = make_triangle (v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], rotate, eps);
      if (!degenerate(t)) {
         mesh->push_back (t);
      }
      else {
        ndegenerated++;
      }
    }
    fclose(f);
  }
  cout << "number of degenerated triangles = " << ndegenerated << endl;
  return 0;
}

int checkASCIIFile (const char *fileName) {
  string line1, line2;
  ifstream input(fileName);
  if (!input.good()) {
    return -1;
  }
  getline(input, line1);
  getline(input, line2);
  if (line1.find("solid")!=string::npos && line2.find("facet")!=string::npos) {
    return FILE_STL_ASCII;
  }
  if (line1.find("xml")!=string::npos && line2.find("amf")!=string::npos) {
    return FILE_AMF_ASCII;
  }
  return FILE_STL_BIN;
}

int main(int argc, char** argv)
{
  bool chaining = true;

  double eps = 0.004;

  /*Total time:*/
  //clock_t begin = clock();

  char *model;

  if (strcmp(argv[2], "-model") == 0) {
     model = argv[3];
  }

  float max_thickness, delta;

  if (strcmp(argv[4], "-thickness") == 0) {
    max_thickness = atof(argv[5]);
  }  
  else {
    printf("Error: specify the slicing spacing in mm (thickness)!!!\n");
  }

  char *adaptive;

  if (strcmp(argv[6], "-adaptive") == 0) {
    adaptive = argv[7];
  }  

  char *write_option = argv[8];

  char *rotate;

  if (strcmp(argv[9], "-rotate") == 0) {
     rotate = argv[10];
  }

  bool orienting;

  if (strcmp(argv[11], "-orienting_contours") == 0) {
    if (strcmp(argv[12],"true") == 0) {
      orienting = true;
    }
    else {
      orienting = false;
    }
  }

  TriangleMesh mesh;
    
  switch (checkASCIIFile(model)) {
    case FILE_STL_ASCII:
      if (stlToMeshInMemory (model, &mesh, false, rotate, eps) != 0)
        return 1;
      break;
    case FILE_STL_BIN:
      if (stlToMeshInMemory (model, &mesh, true, rotate, eps) != 0)
        return 1;
      break;
    default:
      cerr << "Unexpected error" << endl;
      return 1;
  }

  std::string path = model;
    
  std::string lastFileName;
    
  size_t pos = path.find_last_of("/");
    
  if (pos != std::string::npos)
    lastFileName.assign(path.begin() + pos + 1, path.end());
  else
    lastFileName = path;
    
  vector<float> P = compute_planes (&mesh, max_thickness, adaptive, eps, &delta);
  
  int nplanes = P.size();

  vector<contour> polygons[nplanes];
 
  bool srt = false;
  
  if (strcmp(write_option, "-No") == 0) {
     chaining = false;
  }

  if (strcmp(argv[1], "-Trivial") == 0) {
    //TrivialSlicing (&mesh, P, polygons, chaining, orienting);
  } 
  else if (strcmp(argv[1], "-Park") == 0) {
    //ParkSlicing (&mesh, P, polygons, chaining);
  }
  else if (strcmp(argv[1], "-Incremental") == 0) {
    if (strcmp(adaptive, "true") == 0) {
      delta = 0.0;
    }
    IncrementalSlicing (&mesh, P, delta, srt, polygons, chaining, orienting);
  }

  return 0;
}