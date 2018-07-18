// An example project for COMP9024 (UNSW)
// The project includes a structure (adjacency list graph representation)
// and functionalility / algorithms (insertion, deletion, shortest path etc)

// Heap structure from earlier task scheduler is imported
// in order to utilise the heap based priority queue


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

////// STRUCTURES      //////////

// A vertex is a 2D point
typedef struct Vertex { 
	int x; // x-coordinate
	int y; // y-coordinate
} Vertex;

// each edge is a pair of vertices (end-points)
typedef struct Edge {
	Vertex* p1; // first end point
	Vertex* p2; // second end point
} Edge;

typedef struct AdjacentNode {
	struct VertexNode* vnode;
	struct AdjacentNode* next;
	int visited;
} AdjacentNode;


// A vertex node stores a vertex and other information, and you need to expand this type
typedef struct VertexNode {
	struct Vertex* v;
	struct AdjacentNode* adj;
	struct AdjacentNode* lastadj;
	struct VertexNode* next;
	struct HeapNode* h;
	int visited;
	double min_distance;
	struct VertexNode* prior;
} VertexNode;

typedef struct GraphRep *Graph;
typedef struct GraphRep { // graph header
	VertexNode *vertices; // pointer to head of vertices
	VertexNode *lastvert; // pointer to last vertex
	int nV; // #vertices
	int nE; // #edges
} GraphRep;


typedef struct SimpleQueue {
	struct SimpleQueueNode* head;
} SimpleQueue;

typedef struct SimpleQueueNode {
	struct VertexNode* vnode;
	struct SimpleQueueNode* next;
} SimpleQueueNode;
////// END STRUCTURES   //////////


////// HEAP STRUCTURES FROM ASSIGNMENT 3 ///////////

typedef struct HeapNode { 
	// each node stores the priority (key), name, execution time,
	//  release time and deadline of one task
	double key; //key of this item
	struct VertexNode *vnode;  // task name 
	struct HeapNode *parent; //pointer to parent
	struct HeapNode *left; //pointer to left child
	struct HeapNode *right; //pointer to right child
} HeapNode;

//data type for a priority queue (heap) 
typedef struct Heap{ //this is heap header
	int  size;      // count of items in the heap
	HeapNode *LastNode; // last node pointer 
	HeapNode *root; // pointer to the root of heap
} Heap;

// create a new heap node to store an item (task) 
HeapNode *newHeapNode(double k, VertexNode *node, HeapNode *L, HeapNode *R, HeapNode *P){
  // L, R, L: pointers to left child, right child and parent, respectively 	 
	HeapNode *new;
	new = malloc(sizeof(HeapNode));
	assert(new != NULL);
	new->key = k;
	new->vnode = node;
	new->left = L; // left child
	new->right = R; // righ child
	new->parent = P; // parent
	return new;
}

// create a new empty heap-based priority queue
 Heap *newHeap()
{ // this function creates an empty heap-based priority queue
	Heap *T;
	T = malloc(sizeof(Heap));
	assert (T != NULL);
	T->size = 0;
	T->LastNode=NULL; 
	T->root = NULL;
	return T;
}

// helper function that bubbles up the value inserted
// if the value is less than the parent it will swap the positions of the node with its parent
// recursively bubbles through parent stream

HeapNode *BubbleUp(Heap*T, HeapNode *node){
	if (!node->parent){ //base case trying to swap at root
		T->root = node;
		node->vnode->h = node;
		return node;
	}
	else{
		if (node->parent->key > node->key){
			if (node->parent->left == node){ // case of child is on left
				HeapNode *temp = newHeapNode(node->key,node->vnode,node->parent,node->parent->right,NULL);
				temp->left->left = node->left;
				temp->left->right = node->right;
				if (temp->left->left){ temp->left->left->parent = temp->left;}
				if (temp->left->right){ temp->left->right->parent = temp->left;}
				if (temp->right){temp->right->parent = temp;}
				if (node->parent->parent){
					temp->parent = node->parent->parent;
					if (node->parent->parent->left == node->parent) {node->parent->parent->left = temp;}
					else{ node->parent->parent->right = temp; }
				}
				temp->left->parent = temp;
				// free(node);
				node = temp;
				node->vnode->h = node;
				BubbleUp(T,node); 
				return node->left;
			}
			else if (node->parent->right == node){ // case of child is on right
				HeapNode *temp = newHeapNode(node->key,node->vnode,node->parent->left,node->parent,NULL);
				temp->right->right = node->right;
				temp->right->left = node->right;
				if (temp->right->right){ temp->right->right->parent = temp->right;}
				if (temp->right->left){ temp->right->left->parent = temp->left;}
				if (temp->left){temp->left->parent = temp;}
				if (node->parent->parent){
					temp->parent = node->parent->parent;
					if (node->parent->parent->left == node->parent) {node->parent->parent->left = temp;}
					else{ node->parent->parent->right = temp; }
				}
				temp->right->parent = temp;
				// free(node);
				node = temp;
				node->vnode->h = node;
				BubbleUp(T,node);
				return node->right;
			}
		}
		else{ // doesn't need to be swapped
			node->vnode->h = node;
			return node; 
		}
	}
}


HeapNode *BubbleDown(Heap*T, HeapNode *node){
	if (!node->left && !node->right){ //base case trying to swap at root
		node->vnode->h = node;
		return node;
	}
	else{
		if (node->left->key > node->key){
			HeapNode *temp = newHeapNode(node->left->key,node->left->vnode,node,node->right,NULL);
			node->right = node->left->right;
			node->left = node->left->left;
			if (node->parent){
				if (node->parent->left == node){ node->parent->left = temp;}
				else {node->parent->right = temp;}
			}
			else {T->root = temp;}
			node->parent = temp;
			temp->vnode->h = temp;
			BubbleDown(T,temp);
			return node;
		}
		else if ( node->right->key > node->key){
			HeapNode *temp = newHeapNode(node->right->key,node->right->vnode,node->left,node,NULL);
			node->left = node->right->left;
			node->right = node->right->right;
			if (node->parent){
				if (node->parent->left == node){ node->parent->left = temp;}
				else {node->parent->right = temp;}
			}
			else {T->root = temp;}
			node->parent = temp;
			temp->vnode->h = temp;
			BubbleDown(T,temp);
			return node;
		}
	}
}


// Function to find the next last node (after last node)
// Used to replace LastNode
// Time complexity worst case is log(n)
HeapNode *FindNextLastNode(Heap *T){
	HeapNode *node;
	if (T->size <= 1 ){
		return NULL;
	}
	if (T->size == 2){
		return T->root;
	}
	else{
		if (T->LastNode == T->LastNode->parent->right){	
			return T->LastNode->parent->left;
		}
		else{	//left node 
			// moves up tree until finds a right node, then goes up one
			// and down the left subtree, then all the way down right node until leaf
			node = T->LastNode->parent;
			while (node->parent != NULL){
				if (node == node->parent->right){
					node = node->parent->left; // traverse around the point
					break;
				}
				node = node->parent;

			}
			
			while (node->right != NULL){
				node = node->right;
			}
			return node;
		}
	}
	
}


// put the time complexity analysis for Insert() here  
// Time complexity is O(log n) for insertion
// This is due to bubble up 
HeapNode *Insert(Heap *T, double k, VertexNode *VertNode)
{ // k: key, n: task name, c: execution time, r: release time, d:deadline 
	HeapNode *node = newHeapNode(k,VertNode,NULL,NULL,NULL),*temp;	
	if (T->size == 0){	// Base cases for small tree
		HeapNode *node = newHeapNode(k,VertNode,NULL,NULL,NULL);	
		T->root = node;
	}
	else if (T->size == 1 ){ // case of only 1 node 
		node->parent = T->root;
		if (T->size == 1) { T->root->left = node ; }
		else { T->root->right = node ; }
	}
	else {
		if (!T->LastNode->parent->right){
			T->LastNode->parent->right = node;
			node->parent = T->LastNode->parent;
		}
		else{	//move up until find root or leaf node
			while(1){	//temporarily moves the last node to move through tree)
				if (T->LastNode == T->root){ // hits the root i.e.
					while (T->LastNode->left){
						T->LastNode = T->LastNode->left;
					}
					T->LastNode->left = node;
					node->parent = T->LastNode;
					break;
				}
				else if (T->LastNode->parent->left == T->LastNode){ 	//Left leaf
					T->LastNode = T->LastNode->parent->right;
					while (T->LastNode->left) {
						T->LastNode = T->LastNode->left;
					}
					T->LastNode->left = node;
					node->parent = T->LastNode;
					break;
				}
				T->LastNode = T->LastNode->parent;
			}
		}
	}
	T->LastNode = BubbleUp(T,node);
	T->size++;
	return node;
}


// returns 1 if right child is smaller/equal or 0 if left is smaller
int FindMinChild(HeapNode *node){ 
	int i;
	if (node->right == NULL){
		i=0;
	}
	else{
		if (node->left->key >= node->right->key){
			i=1;
		}
		else{
			i=0;
		}

	}
	return i;
}

int NumChildNodes(HeapNode *node){
	int i;
	if (!node->left && !node->right){
		i=0;
	}
	else if (!node->right){
		i=1;
	}
	else{
		i=2;
	}
	return i;
}

// recursively bubbles down tree until it hits the end
// then moves last node to the current position and updates last node
// exception cases exists for example if the last node is the current node
HeapNode *BubbleOut(Heap*T, HeapNode *node){
	int ChildNodes = NumChildNodes(node); // number of child nodes
	if (ChildNodes == 0){
		if (node->parent == NULL){ 	// root case

			T->size--;
			T->root = NULL;
			return node;
		}
		else{
			HeapNode *NextLastNode = FindNextLastNode(T);
			T->size--;
			if (T->LastNode == node){
				T->LastNode = NextLastNode;
				return NULL;
			}
			else{
				HeapNode *NewNode = newHeapNode(T->LastNode->key,T->LastNode->vnode,NULL,NULL,NULL);
				if (T->LastNode->parent->left == T->LastNode){T->LastNode->parent->left = NULL;}
				else {T->LastNode->parent->right = NULL;}
				if (NextLastNode == node){T->LastNode = NewNode;} // case when equal to NewNode
				else {T->LastNode = NextLastNode;}

				return NewNode;
			}
		}
		
	}
	else if (ChildNodes == 1){
		HeapNode *NewNode = newHeapNode(node->left->key,node->left->vnode,NULL,NULL,node->parent);
		T->LastNode = FindNextLastNode(T);
		T->size--;
		if (!NewNode->parent){T->root = NewNode;}
		return NewNode;
	}
	else if (ChildNodes == 2){
		int IndicatorMin = FindMinChild(node);
		if (IndicatorMin == 0){
			HeapNode *NewNode;
			if (T->LastNode != node->right){NewNode = newHeapNode(node->left->key,node->left->vnode,BubbleOut(T,node->left),node->right,node->parent);}
 			else {NewNode = newHeapNode(node->left->key,node->left->vnode,BubbleOut(T,node->left),NULL,node->parent);}
 			if (NewNode->right){NewNode->right->parent = NewNode;}
 			if (NewNode->left){NewNode->left->parent = NewNode;}
 			if (!NewNode->parent){T->root = NewNode;}
 			return NewNode;
 		}
		else if (IndicatorMin == 1){
 			HeapNode *temp, *NewNode = newHeapNode(node->right->key,node->right->vnode,node->left,BubbleOut(T,node->right),node->parent);
 			if (NewNode->left){NewNode->left->parent = NewNode;}
 			if (NewNode->right){NewNode->right->parent = NewNode;}
 			if (!NewNode->parent){T->root = NewNode;}
 			NewNode->vnode->h = NewNode;
 			return NewNode;
 		}

	}

}

// put your time complexity for RemoveMin() here:
// Time complexity of Remove min will be O(log n)
HeapNode *RemoveMin(Heap *T){
	HeapNode *node = T->root;
	BubbleOut(T,node);	
 	return node;
}



////// END HEAP STRUCTURES FROM ASSIGNMENT 3 ///////////


////// HELPER FUNCTIONS //////////



Vertex *newVertex(int X, int Y){
	struct Vertex* v = malloc(sizeof(struct Vertex));
	assert (v != NULL);
	v->x = X;
	v->y = Y;
	return v;
}

AdjacentNode *newAdjacentNode(VertexNode *vert){
	struct AdjacentNode* node = malloc(sizeof(struct AdjacentNode));
	assert (node != NULL);
	node->vnode = vert;
	node->next = NULL;
	node->visited = 0;
	return node;
}

Edge *newEdge(Vertex *v1,Vertex *v2){
	struct Edge* e = malloc(sizeof(struct Edge));
	assert (e != NULL);
	e->p1 = v1;
	e->p2 = v2;
	return e;
}

VertexNode *newVertexNode(Vertex *vert)
{
    struct VertexNode* newNode = malloc(sizeof(struct VertexNode));
    assert(newNode!=NULL);
    newNode->v = vert;
    newNode->next = NULL;
    newNode->adj = NULL;
    newNode->lastadj = NULL;
    newNode->visited = 0;
    newNode->h = NULL;
    newNode->min_distance = 0;
    newNode->prior = NULL;
    return newNode;
}

// Time complexity is O(1);
Graph CreateEmptyGraph() {
	Graph g = malloc(sizeof(GraphRep));
	assert(g != NULL);
	g->vertices = NULL; g->lastvert = NULL;
	g->nV = 0; g->nE = 0;
	return g;
}

SimpleQueueNode *newQueueNode(VertexNode *vert){
	struct SimpleQueueNode* qnode = malloc(sizeof(struct SimpleQueueNode));
	assert (qnode != NULL);
	qnode->vnode = vert;
	return qnode;
}


SimpleQueue *newQueue(SimpleQueueNode *node){
	struct SimpleQueue* q = malloc(sizeof(struct SimpleQueue));
	assert (q != NULL);
	q->head = node;
	return q;
}
///// END HELPER FUNCTIONS ///////

///// TRAVERSING FUNCTIONS ///////

// Checks whether the vertex is in graph g
// if not returns NULL
VertexNode *CheckVertexInGraph(Graph g,Vertex* vert){
	VertexNode* node = g->vertices;
	if (node){
		while (node!= NULL){
			if (node->v->x == vert->x && node->v->y == vert->y){
				return node;
			}
			node = node->next;
		}
	}
	return NULL;
}


// Checks whether the vertex is the adjacent list of a vertex node
// if not returns NULL
AdjacentNode *CheckVertexInAdj(VertexNode* original,Vertex* vert){
	AdjacentNode* node = original->adj;
	if (node){
		while (node!= NULL){
			if (node->vnode->v->x == vert->x && node->vnode->v->y == vert->y){
				return node;
			}
			node = node->next;
		}
	}
	return NULL;
}


//// END TRAVERSING FUNCTIONS ////



////  MAIN FUNCTIONS 	   ///////

// Time complexity
// worst case is O(m+n)
int InsertEdge(Graph g, Edge *e) {
	VertexNode* v1;		VertexNode* v2;
	AdjacentNode* adj1; AdjacentNode* adj2;

	v1 = CheckVertexInGraph(g,e->p1);	v2 = CheckVertexInGraph(g,e->p2);
	if (v1 && v2) { // case of both vertex in LL
		adj1 = CheckVertexInAdj(v1,e->p2);
		adj2 = CheckVertexInAdj(v2,e->p1);
		if (adj1 || adj2) { return 0; } // edge is in graph
		else{
			adj1 = newAdjacentNode(v2);	adj2 = newAdjacentNode(v1);
			v1->lastadj->next = adj1;	v2->lastadj->next = adj2;
			v1->lastadj = adj1;			v2->lastadj = adj2;
			g->nE += 1;
		}
	}
	else if (v1){	// case of vertex 1 in LL but not vertex 2
		v2 = newVertexNode(e->p2);
		g->lastvert->next = v2;		g->lastvert = v2;	
		adj1 = newAdjacentNode(v2);		adj2 = newAdjacentNode(v1);
		v2 -> adj = v2->lastadj = adj2;	
		v1 -> lastadj-> next = adj1;  v1 -> lastadj = adj1;
		g->nV += 1;		g->nE += 1;
	}
	else if (v2){	// case of vertex 1 in LL but not vertex 2
		v1 = newVertexNode(e->p1);
		g->lastvert->next = v1;		g->lastvert = v1;	
		adj1 = newAdjacentNode(v2);		adj2 = newAdjacentNode(v1);
		v1 -> adj = v1->lastadj = adj1;	
		v2 -> lastadj-> next = adj2;  v2 -> lastadj = adj1;
		g->nV += 1;		 g->nE += 1;
	}
	else{ // both vertex not in list
		v1 = newVertexNode(e->p1);		v2 = newVertexNode(e->p2);
		adj1 = newAdjacentNode(v2);		adj2 = newAdjacentNode(v1);
		v1->next = v2;
		v1->adj = adj1;		v1->lastadj = adj1;
		v2->adj = adj2;		v2->lastadj = adj2;
		g->nV += 2;			g->nE += 1;
		if (g->vertices == NULL){ g->vertices = v1; }
		else { g->lastvert->next = v1; }
		g->lastvert = v2;
	}
	return 1;
}



// helper function to remove adj from LL
// 0 means vertex not found on node
// 1 means vertex is found and adj is deleted from LL
// vert must be corresponding vertex
int deleteAdj(VertexNode *vertnode, Vertex *vert){ 
	AdjacentNode* adjnode = vertnode->adj;
	if (adjnode->vnode->v->x == vert->x && adjnode->vnode->v->y == vert->y){//
		if (adjnode->next == NULL){ vertnode -> adj = NULL; vertnode -> lastadj = NULL; }
		else { vertnode->adj = adjnode->next; }
		free(adjnode);
		return 1;
	}
	else{
		while(adjnode->next != NULL){
			if (adjnode->next->vnode->v->x == vert->x && adjnode->next->vnode->v->y == vert->y){ 
				if (adjnode->next == vertnode->lastadj) { vertnode->lastadj = adjnode; }
				AdjacentNode* temp = adjnode->next;
				adjnode->next = adjnode->next->next;
				free(temp);
				return 1;
			}
			adjnode= adjnode->next;
		}
	}
	return 0;
}

// helper function to delete a vertice (if all adj are gone)
// Already knows vert is in graph
// 0 if graph has no vertices left
// 1 if not to be deleted
int deleteVert(Graph g, VertexNode *goalnode){ 
	VertexNode *vertnode = g-> vertices;
	if (vertnode == goalnode){

		if (vertnode->next){ g->vertices = vertnode->next; free(vertnode); return 1;}
		else{g->vertices = NULL; free(vertnode); return 0;}
	}
	else{ 
		while (vertnode->next != NULL){
			if (vertnode->next->v->x == goalnode->v->x && vertnode->next->v->y == goalnode->v->y){
				if (g->lastvert == goalnode) { g->lastvert = vertnode; }
				VertexNode *temp = vertnode->next;
				vertnode -> next = goalnode->next;
				free(temp);
				return 1;
			}
			vertnode = vertnode->next;
		}
	}
	
}


// Time complexity
// O(m+n)
// traverses to find vertices;
// traverses edges to find adjacency list
// traverses to delete vertices;
void DeleteEdge(Graph g, Edge *e){
	VertexNode* v1;		VertexNode* v2;
	AdjacentNode* adj1; AdjacentNode* adj2;
	v1 = CheckVertexInGraph(g,e->p1);	v2 = CheckVertexInGraph(g,e->p2);
	if (v1 && v2){
		int v1_indicator, v2_indicator;
		if (deleteAdj(v1,e->p2)== 1 && deleteAdj(v2,e->p1) == 1){
			g->nE -= 1;
			if (v1->adj == NULL){
				deleteVert(g,v1);
				g->nV -= 1;
			}
			if (v2->adj == NULL){
				deleteVert(g,v2);
				g->nV -=1;
			}
		}


	}
}






// BFS alogirthm takes time complexity of O(m+n)
int BreadSearch(int arr[][2],Graph g, VertexNode *node){
	
	SimpleQueueNode *qnode = newQueueNode(node);
	SimpleQueueNode *temp;
	SimpleQueueNode *lastnode = qnode;
	SimpleQueue *queue = newQueue(qnode);
	AdjacentNode *adjtrav;
	int lastcord = 0;
	
	node->visited = 1;
	arr[lastcord][0] = qnode->vnode->v->x;
	arr[lastcord][1] = qnode->vnode->v->y;
	while(qnode != NULL){
		adjtrav = qnode->vnode->adj;
		while (adjtrav != NULL){ 
			if(adjtrav -> vnode->visited == 0){ 
				lastcord = lastcord +1;
				arr[lastcord][0] = adjtrav->vnode->v->x;
				arr[lastcord][1] = adjtrav->vnode->v->y;
				lastnode->next = newQueueNode(adjtrav->vnode); 
				lastnode = lastnode->next;
				adjtrav->vnode->visited = 1;
			}
			adjtrav = adjtrav->next;
		}
		qnode = qnode->next;
	}
	qnode = queue->head;
	return lastcord;

}


// function to compare tuple in qsort
int compare(const void *x, const void *y){
	int (*a)[2] = (int(*)[2])x;
	int (*b)[2] = (int(*)[2])y;
	if ((*a)[0]==(*b[0])){
		return (*a)[1]-(*b)[1];
	}
	
	else{
		return (*a)[0]-(*b)[0];
	}
	
}


// Time complexity:
// Reachable vertices using BFS which will take O(m+n) time complexity
// qsort is then used to order these vertices
// the worst case of qsort is O(n^2) but average is O(nlog(n))
// therefore time complexity worst case is O(n^2) and average is O(nlogn)

void ReachableVertices(Graph g, Vertex *v)
{
	VertexNode* desired=NULL; VertexNode* temp = g->vertices;
	AdjacentNode* adjtemp; 
	if (!temp) {return;}
	while(temp != NULL){ // set all visitied to zero
		temp->visited = 0;
		if (temp->v->x == v->x && temp->v->y == v->y){
			desired = temp;
		}
		temp = temp->next;
	}
	if (!desired) {printf("Vertex not in graph"); return;}
	int coords[g->nV][2]; 
	int lastdigit = BreadSearch(coords,g,desired);
	qsort(coords,lastdigit+1,sizeof(coords[0]),compare);
	for (int i=1;i<=lastdigit;i++){
		if (i < lastdigit){ printf("(%d,%d),",coords[i][0],coords[i][1]);}	
		else{printf("(%d,%d)\n",coords[i][0],coords[i][1]);}
	}
	return;
}


// Time complexity is O(m+n)
// every edge and vertex is set to free
void FreeGraph(Graph g){
	VertexNode* tempnode = g->vertices;
	VertexNode* pastnode;
	AdjacentNode* tempadj;
	AdjacentNode* pastadj;
	while(tempnode != NULL){
		tempadj = tempnode->adj;
		free(tempnode->v);
		while(tempadj != NULL){
			pastadj = tempadj;
			tempadj = tempadj->next;
			free(pastadj);
		}
		pastnode = tempnode;
		tempnode = tempnode->next;
		free(pastnode);
	}
	free(tempnode);
	free(tempadj);
	free(g);

}



// Time complexity is O(m+n);
// every edge and vertex will be expanded
// Prints in breadth first order 
// picks first vertex in graph and expands all edges setting visited values
// iterates through all edges if they're unconnected to this
// continues until all vertices have been expanded. 
void ShowGraph(Graph g)
{
	VertexNode *temp = g->vertices;
	AdjacentNode *adjtemp;
	AdjacentNode *adjtrav;
	VertexNode *node = g->vertices;
	SimpleQueueNode *qnode = newQueueNode(node);
	SimpleQueueNode *qtemp;
	SimpleQueueNode *lastnode = qnode;
	SimpleQueue *queue = newQueue(qnode);
	if (!temp) {return;}
	while(temp != NULL){ // set all visitied to zero
		temp->visited = 0;
		adjtemp = temp->adj;
		while(adjtemp != NULL){
			adjtemp->visited =0;
			adjtemp = adjtemp->next;
		}
		temp = temp->next;
	}
	while(qnode != NULL){
		if (qnode != NULL){
			if (qnode->vnode->visited == 0){ 
				adjtrav = qnode->vnode->adj;
				while (adjtrav != NULL){ 
					if (adjtrav->visited == 0){
						printf("(%d,%d),(%d,%d) ", qnode->vnode->v->x,qnode->vnode->v->y,adjtrav->vnode->v->x,adjtrav->vnode->v->y);
						adjtemp = adjtrav->vnode->adj;
						while (adjtemp->vnode!=qnode->vnode){
							adjtemp = adjtemp->next;
						}
						adjtemp->visited=1;
						adjtrav->visited=1;
						lastnode->next = newQueueNode(adjtrav->vnode); 
						lastnode = lastnode->next;
					}
					adjtrav=adjtrav->next;
				}
			}
			qnode->vnode->visited = 1;
			qnode = qnode->next;
		}
		if (qnode == NULL){
			if (node->next){
				node = node->next;
				lastnode->next = newQueueNode(node);
				lastnode = lastnode->next;
				qnode = lastnode;
			}
		}
	}
	printf("\n");
	return;
}

// Function to calculate distance between two points
double distance(VertexNode *v1, VertexNode *v2){
	double x1 = v1->v->x, x2 = v2->v->x, y1 = v1->v->y , y2= v2->v->y;
	double x_dif = (x1-x2) * (x1-x2);
	double y_dif = (y1-y2) * (y1-y2);
	double sum = x_dif + y_dif;
	double value = sqrt(sum);
	return value;
}


// Time complexity:
// O((m+n)log n) where m and n are the number of edges & vertices in graph respectively
// Due to implementation of Dijkstra's shortest path algorithm using heap based priority que
void ShortestPath(Graph g, Vertex *u, Vertex *v){
	VertexNode *temp = g->vertices,*A,*B;
	AdjacentNode *adjtemp;
	Heap *que = newHeap();
	HeapNode *heaptemp;
	double length;
	int PathStore[g->nV][2];
	int counter = 0;
	// resets all variables and finds A and B
	while(temp != NULL){ 
		temp->h = NULL;
		temp->visited = 0;
		temp->min_distance = 0;
		temp->prior = NULL;
		if (temp->v->x == u->x && temp->v->y == u->y){
			A = temp;
		}
		if (temp->v->x == v->x && temp->v->y == v->y){
			B = temp;
		}
		temp=temp->next;
	}
	// Checks for valid input
	if (!A || !B || (A->v->x == B->v->x && A->v->y == B->v->y)){ return; }
	Insert(que,0,A);
	// loops through each vertex and expands all edges if vertex hasn't been visisted
	// continually recalculates the minimum value 
	while(que->root != NULL){
		heaptemp = RemoveMin(que);
		if (heaptemp->vnode->visited == 0){
			heaptemp ->vnode->visited = 1;
			adjtemp = heaptemp->vnode->adj;
			while (adjtemp != NULL){
				length = distance(heaptemp->vnode,adjtemp->vnode);
				if (adjtemp->vnode != A){
					if (adjtemp->vnode->min_distance == 0){
						Insert(que,length+heaptemp->vnode->min_distance,adjtemp->vnode);
						adjtemp->vnode->min_distance = length+heaptemp->vnode->min_distance;
						adjtemp->vnode->prior = heaptemp->vnode;
					}
					else if ((heaptemp->vnode->min_distance + length) < adjtemp->vnode->min_distance ){
						if (!adjtemp->vnode->h && heaptemp->vnode->min_distance + length < adjtemp->vnode->min_distance){
							Insert(que,length+adjtemp->vnode->min_distance,adjtemp->vnode);
						}
						else {
							adjtemp->vnode->h->key = length + adjtemp->vnode->min_distance;
							BubbleUp(que,adjtemp->vnode->h);
							BubbleDown(que,adjtemp->vnode->h);
						adjtemp->vnode->min_distance = length+heaptemp->vnode->min_distance;
						adjtemp->vnode->prior = heaptemp->vnode;
						}
					}
				}		
				adjtemp = adjtemp->next;
			}

		}
		heaptemp->vnode->h = NULL;
	}
	temp = B;
	// small loop to back to original path and print in reverse order
	if (B->prior != NULL){
		while (temp != NULL){
			PathStore[counter][0] = temp->v->x;
			PathStore[counter][1] = temp->v->y;
			counter++;
			temp = temp->prior;
		}
		counter--;
		for (counter;counter>=0;counter--){
			if (counter != 0){printf("(%d,%d),",PathStore[counter][0],PathStore[counter][1]);}
			else{printf("(%d,%d)\n",PathStore[counter][0],PathStore[counter][1]);}
		}
	}
}


//  END MAIN FUNCTIONS   ///////


// TO DELETE //
void printgraph(Graph g){
	VertexNode* node = g->vertices;
	AdjacentNode* adjnode;
	while ( node!= NULL ){
		printf("(%d,%d) -> ",node->v->x,node->v->y);
		adjnode = node->adj;
		while(adjnode!=NULL){
			printf("[%d,%d] - ",adjnode->vnode->v->x,adjnode->vnode->v->y);
			adjnode = adjnode->next;
		}
		printf("\n");
		node = node->next;
	}
}
// END DELETE


int main() //sample main for testing 
{Graph g1;
 Edge *e_ptr; 
 Vertex *v1, *v2;
  
 // Create an empty graph g1;
 g1=CreateEmptyGraph();
  
 // Create first connected component 
 // Insert edge (0,0)-(0,10)
 e_ptr = (Edge*) malloc(sizeof(Edge));
 assert(e_ptr != NULL);
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v2=(Vertex *) malloc(sizeof(Vertex));
 assert(v2 != NULL);
 v1->x=0;
 v1->y=0;
 v2->x=0;
 v2->y=10;
 e_ptr->p1=v1;
 e_ptr->p2=v2;
 if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
 
 // Insert edge (0,0)-(5,6)
 e_ptr = (Edge*) malloc(sizeof(Edge));
 assert(e_ptr != NULL);
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v2=(Vertex *) malloc(sizeof(Vertex));
 assert(v2 != NULL);
 v1->x=0;
 v1->y=0;
 v2->x=5;
 v2->y=6;
 e_ptr->p1=v1;
 e_ptr->p2=v2;
 if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
 
 // Insert edge (0, 10)-(10, 10)
 e_ptr = (Edge*) malloc(sizeof(Edge));
 assert(e_ptr != NULL);
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v2=(Vertex *) malloc(sizeof(Vertex));
 assert(v2 != NULL);
 v1->x=0;
 v1->y=10;
 v2->x=10;
 v2->y=10;
 e_ptr->p1=v1;
 e_ptr->p2=v2;
 if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
 
 // Insert edge (0,10)-(5,6)
 e_ptr = (Edge*) malloc(sizeof(Edge));
 assert(e_ptr != NULL);
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v2=(Vertex *) malloc(sizeof(Vertex));
 assert(v2 != NULL);
 v1->x=0;
 v1->y=10;
 v2->x=5;
 v2->y=6;
 e_ptr->p1=v1;
 e_ptr->p2=v2;
 if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
 
 // Insert edge (0,0)-(5,4)
 e_ptr = (Edge*) malloc(sizeof(Edge));
 assert(e_ptr != NULL);
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v2=(Vertex *) malloc(sizeof(Vertex));
 assert(v2 != NULL);
 v1->x=0;
 v1->y=0;
 v2->x=5;
 v2->y=4;
 e_ptr->p1=v1;
 e_ptr->p2=v2;
 if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
 
 // Insert edge (5, 4)-(10, 4)
 e_ptr = (Edge*) malloc(sizeof(Edge));
 assert(e_ptr != NULL);
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v2=(Vertex *) malloc(sizeof(Vertex));
 assert(v2 != NULL);
 v1->x=5;
 v1->y=4;
 v2->x=10;
 v2->y=4;
 e_ptr->p1=v1;
 e_ptr->p2=v2;
 if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
 
 // Insert edge (5,6)-(10,6)
 e_ptr = (Edge*) malloc(sizeof(Edge));
 assert(e_ptr != NULL);
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v2=(Vertex *) malloc(sizeof(Vertex));
 assert(v2 != NULL);
 v1->x=5;
 v1->y=6;
 v2->x=10;
 v2->y=6;
 e_ptr->p1=v1;
 e_ptr->p2=v2;
 if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
 
 // Insert edge (10,10)-(10,6)
 e_ptr = (Edge*) malloc(sizeof(Edge));
 assert(e_ptr != NULL);
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v2=(Vertex *) malloc(sizeof(Vertex));
 assert(v2 != NULL);
 v1->x=10;
 v1->y=10;
 v2->x=10;
 v2->y=6;
 e_ptr->p1=v1;
 e_ptr->p2=v2;
 if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
 
 // Insert edge (10, 6)-(10, 4)
 e_ptr = (Edge*) malloc(sizeof(Edge));
 assert(e_ptr != NULL);
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v2=(Vertex *) malloc(sizeof(Vertex));
 assert(v2 != NULL);
 v1->x=10;
 v1->y=6;
 v2->x=10;
 v2->y=4;
 e_ptr->p1=v1;
 e_ptr->p2=v2;
 if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
 
 // Create second connected component
 // Insert edge (20,4)-(20,10)
 e_ptr = (Edge*) malloc(sizeof(Edge));
 assert(e_ptr != NULL);
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v2=(Vertex *) malloc(sizeof(Vertex));
 assert(v2 != NULL);
 v1->x=20;
 v1->y=4;
 v2->x=20;
 v2->y=10;
 e_ptr->p1=v1;
 e_ptr->p2=v2;
 if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
 
 // Insert edge (20,10)-(30,10)
 e_ptr = (Edge*) malloc(sizeof(Edge));
 assert(e_ptr != NULL);
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v2=(Vertex *) malloc(sizeof(Vertex));
 assert(v2 != NULL);
 v1->x=20;
 v1->y=10;
 v2->x=30;
 v2->y=10;
 e_ptr->p1=v1;
 e_ptr->p2=v2;
 if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
 
 // Insert edge (25,5)-(30,10) 	
 e_ptr = (Edge*) malloc(sizeof(Edge));
 assert(e_ptr != NULL);
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v2=(Vertex *) malloc(sizeof(Vertex));
 assert(v2 != NULL);
 v1->x=25;
 v1->y=5;
 v2->x=30;
 v2->y=10;
 e_ptr->p1=v1;
 e_ptr->p2=v2;
 if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n"); 
 
 //Display graph g1
 ShowGraph(g1);
	
 // Find the shortest path between (0,0) and (10,6) 
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v2=(Vertex *) malloc(sizeof(Vertex));
 assert(v2 != NULL);
 v1->x=0;
 v1->y=0;
 v2->x=10;
 v2->y=6;
 ShortestPath(g1, v1, v2);
 free(v1);
 free(v2);
	  
 // Delete edge (0,0)-(5, 6)
 e_ptr = (Edge*) malloc(sizeof(Edge));
 assert(e_ptr != NULL);
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v2=(Vertex *) malloc(sizeof(Vertex));
 assert(v2 != NULL);
 v1->x=0;
 v1->y=0;
 v2->x=5;
 v2->y=6;
 e_ptr->p1=v1;
 e_ptr->p2=v2; 	 
 DeleteEdge(g1, e_ptr);
 free(e_ptr);
 free(v1);
 free(v2);
 	 
 // Display graph g1
 ShowGraph(g1);
	
 // Find the shortest path between (0,0) and (10,6) 
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v2=(Vertex *) malloc(sizeof(Vertex));
 assert(v2 != NULL);
 v1->x=0;
 v1->y=0;
 v2->x=10;
 v2->y=6; 
 ShortestPath(g1, v1, v2);
 free(v1);
 free(v2);
 
 // Find the shortest path between (0,0) and (25,5)
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v2=(Vertex *) malloc(sizeof(Vertex));
 assert(v2 != NULL);
 v1->x=0;
 v1->y=0;
 v2->x=25;
 v2->y=5;
 ShortestPath(g1, v1, v2);
 free(v1);
 free(v2);	
 
 // Find reachable vertices of (0,0)
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v1->x=0;
 v1->y=0;
 ReachableVertices(g1, v1);
 free(v1);
 
 // Find reachable vertices of (20,4)
 v1=(Vertex*) malloc(sizeof(Vertex));
 assert(v1 != NULL);
 v1->x=20;
 v1->y=4;
 ReachableVertices(g1, v1);
 free(v1);
 
 // Free graph g1
 FreeGraph(g1);
 
 return 0; 
}




