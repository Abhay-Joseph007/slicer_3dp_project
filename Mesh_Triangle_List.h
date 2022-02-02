#include<iostream>
#include"TriangleClass.h"

using namespace std;

typedef struct node {
  Triangle t;
  struct node *next;
  struct node *prev;
} Mesh_Triangle_Node_t;

typedef struct _list {
   Mesh_Triangle_Node_t *head;
   Mesh_Triangle_Node_t *tail;
} Mesh_Triangle_List_t;

Mesh_Triangle_List_t* Mesh_Triangle_List_create();

void Mesh_Triangle_List_insert (Triangle, Mesh_Triangle_List_t);

void Mesh_Triangle_List_union (Mesh_Triangle_List_t, Mesh_Triangle_List_t);

Mesh_Triangle_Node_t* Mesh_Triangle_List_remove (Mesh_Triangle_List_t , Mesh_Triangle_Node_t);