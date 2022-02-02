#include"Mesh_Triangle_List.h"
#include"TriangleClass.h"
#include"v3.h"

Mesh_Triangle_List_t* Mesh_Triangle_List_create(){
    Mesh_Triangle_List_t* List_t = (Mesh_Triangle_List_t *)malloc(sizeof(Mesh_Triangle_List_t));
    List_t->head = NULL;
    List_t->tail = NULL;
    return List_t;
}

void Mesh_Triangle_List_insert(Triangle triangle,Mesh_Triangle_List_t* L){

    Mesh_Triangle_Node_t* Node_t = (Mesh_Triangle_Node_t *)malloc(sizeof(Mesh_Triangle_Node_t));

    if(L->head == NULL){
        Node_t->prev = NULL;
        Node_t->next = NULL;
        L->head = Node_t;
        L->tail = Node_t;
    }
    else{
        Node_t->prev = NULL;
        Node_t->next = L->head;
        L->head->prev = Node_t;
        L->head = Node_t;        
    }    
}

void Mesh_Triangle_List_union (Mesh_Triangle_List_t *L1, Mesh_Triangle_List_t *L2){

    if(L1->tail != NULL && L2->head!=NULL){
        L1->tail->next = L2->head;
        L2->head->prev = L1->tail;
        L1->tail = L2->tail;
    }
    else if(L2->head!= NULL){
        L1->head = L2->head;
        L1->tail = L2->tail;
    }
}

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

