/*	In this module, I define all the functions that are, direct or indirectly, useful
	for the main fibration routines needed to a correct implementation for the partition
	refinement algorithm. Also, here I define my own functions to construct directed 
	networks. All the handmade data structures used are defined in 'structforfiber.h'.

	Author: Higor da S. Monteiro - Universidade Federal do Ceará
	Email: higor.monteiro@fisica.ufc.br
	Complex System Lab - Departament of Physics/Universidade Federal do Ceará (UFC)
*/

#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structforfiber.h"

////////////////////////////////////////////////////////////////////
////// Implementation to define a directed, unweighted graph ///////
////////////////////////////////////////////////////////////////////
NodeAdj* createNode(int neigh, int type_link)
{
	NodeAdj* newnode = (NodeAdj*)malloc(sizeof(NodeAdj));
	newnode->neighbor = neigh;
	newnode->next = NULL;
    newnode->type_link = type_link;
	return newnode;
}

extern Graph* createGraph(int N, char* nodenames, int name_bool)
{
	Graph* graph = (Graph*)malloc(sizeof(Graph));
	graph->size = N;

	// Array of NodeAdj's //
	graph->array = (AdjList*)malloc(N*sizeof(AdjList));

	int j;
	for(j=0; j<N; j++)
	{ 
		graph->array[j].head_in = NULL; 
		graph->array[j].head_out = NULL; 
	}

	int nodeID;
	char tempname[60];
	/*	if 'name_bool' is one, then we read the file containing
		the names of each node to assign it to each node in our 
		constructed network.	*/
	if(name_bool==1)
	{
		FILE* NAMES = fopen(nodenames, "r");
		for(j=0; j<N; j++)
		{			
			fscanf(NAMES, "%s\t%d\n", &tempname, &nodeID);
			strcpy(graph->array[nodeID].gene_name, tempname);
		}
		fclose(NAMES);
	}
	return graph;
}

/////////////////////////////////////////////////////////////////
int findroot(int node, int* psite)
{
	int root = node;
	int cont = 0;
	int path[100];                                  
	while(psite[root]>=0)                     
	{
		path[cont++] = root;
		root = psite[root];
	}
	if(psite[node]>=0) psite[node] = root;
	while(cont) psite[path[--cont]] = root;
	return root;
}

void merge(int node1, int root1, int node2, int root2, int* psite)
{
	if(psite[root1] <= psite[root2])
	{
		psite[root1] += psite[root2];
		psite[root2] = root1;
		psite[node2] = root1;
	}
	else
	{
		psite[root2] += psite[root1];
		psite[root1] = root2;
		psite[node1] = root2;
	}
}
//////////////////////////////////////////////////////////////////

/*	Here I not just add the proper edges to the network but I
	dynamically defines its weakly connected components through
	a percolation-like process using disjoint sets operations.	*/
void addEdges(int** edges, int* components, Graph* graph, int* regulator, int nE)
{
	int i, j;
	int root1, root2;
	int num_component = graph->size;
	for(j=0; j<(graph->size); j++) components[j] = -1;
	
	int index, node1, node2, reg;
	for(j=0; j<nE; j++)
	{
		// 'node1' -> 'node2' directed link.
        node1 = edges[j][0];
		node2 = edges[j][1];
		reg = regulator[j];

		///// UNION-FIND OPERATIONS //////
		root1 = findroot(node1, components);
		root2 = findroot(node2, components);
		if(root1!=root2) { merge(node1, root1, node2, root2, components); num_component--; } 
		///////////////////////////////////////////////////

		NodeAdj* newnode1 = createNode(node1, reg);
		NodeAdj* newnode2 = createNode(node2, reg);

		newnode2->next = graph->array[node1].head_out;
		newnode1->next = graph->array[node2].head_in;
		graph->array[node1].head_out = newnode2;
		graph->array[node2].head_in = newnode1;
	}
	graph->num_component = num_component;
}

extern int** defineNetwork(int** edges, int* components, Graph* graph, char* filename)
{
	FILE *EDGE_FILE = fopen(filename, "r");
	if(EDGE_FILE==NULL) printf("ERROR in file reading");

	int i, j;
	char type[20];
	int r = 1;
	int nlink = 0; // number of links.
	while(r) // Calculates the number of lines in the file
	{
		r = fscanf(EDGE_FILE, "%d\t%d\t%s\n", &i, &j, &type);
		if(r==EOF) break;
		nlink++;
	}
	rewind(EDGE_FILE);

	edges = (int**)malloc(nlink*sizeof(int*));
    int* regulator = (int*)malloc(nlink*sizeof(int));
	for(j=0; j<nlink; j++)
	{
		edges[j] = (int*)malloc(2*sizeof(int));
		r = fscanf(EDGE_FILE, "%d\t%d\t%s\n", &edges[j][0], &edges[j][1], &type);
        if(strcmp("positive", type)==0) regulator[j] = 0;
        else if(strcmp("negative", type)==0) regulator[j] = 1;
        else if(strcmp("dual", type)==0) regulator[j] = 2;
        else regulator[j] = -1;
	}
	fclose(EDGE_FILE);
	
    // Defines the network structure and its weakly connected components.
	addEdges(edges, components, graph, regulator, nlink);
    return edges;
}

extern int GetNodeNumber(char* edgefile)
{
	FILE *EDGE_FILE = fopen(edgefile, "r");
	if(EDGE_FILE==NULL) printf("ERROR in file reading");

	int max = -1;
	int i, j;
	char type[20];
	int r = 1;
	int nlink = 0; // number of links.
	while(r) // Calculates the number of lines in the file
	{
		r = fscanf(EDGE_FILE, "%d\t%d\t%s\n", &i, &j, &type);
		if(r==EOF) break;
		if(i>max && i>j) max = i;
		else if(j>max && j>i) max = j;
	}
	return (max+1);
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/////////////////////////// GRAPH OPERATIONS ////////////////////////////////
int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

extern int* GET_INNEIGH(Graph* graph, int node)
{
	int n_in = 0;	
	NodeAdj* NODE;
	for(NODE=graph->array[node].head_in; NODE!=NULL; NODE=NODE->next) n_in++;

	int n_index = 0;
	int* in_neighbors = (int*)malloc(n_in*sizeof(int));
	for(NODE=graph->array[node].head_in; NODE!=NULL; NODE=NODE->next)
		in_neighbors[n_index++] = NODE->neighbor;
	
	qsort(in_neighbors, n_in, sizeof(int), cmpfunc);
	return in_neighbors; 
}

extern int* GET_OUTNEIGH(Graph* graph, int node)
{
	int n_in = 0;	
	NodeAdj* NODE;
	for(NODE=graph->array[node].head_out; NODE!=NULL; NODE=NODE->next) n_in++;

	int n_index = 0;
	int* in_neighbors = (int*)malloc(n_in*sizeof(int));
	for(NODE=graph->array[node].head_out; NODE!=NULL; NODE=NODE->next)
		in_neighbors[n_index++] = NODE->neighbor;
	
	qsort(in_neighbors, n_in, sizeof(int), cmpfunc);
	return in_neighbors; 
}

extern int GETNin(Graph* graph, int node)
{
	int n_in = 0;	
	NodeAdj* NODE;
	for(NODE=graph->array[node].head_in; NODE!=NULL; NODE=NODE->next) n_in++;
	return n_in;	
}

extern int GETinType(Graph* graph, int node, int type)
{
	int n_in = 0;	
	NodeAdj* NODE;
	for(NODE=graph->array[node].head_in; NODE!=NULL; NODE=NODE->next)
		if(NODE->type_link==type) n_in++;
	return n_in;	
}

extern int GETNout(Graph* graph, int node)
{
	int n_out = 0;	
	NodeAdj* NODE;
	for(NODE=graph->array[node].head_out; NODE!=NULL; NODE=NODE->next) n_out++;
	return n_out;	
}

extern int NinREG(Graph* graph, int node, int type)
{
	int n_in = 0;	
	NodeAdj* NODE;
	for(NODE=graph->array[node].head_in; NODE!=NULL; NODE=NODE->next)
		if(NODE->type_link==type) n_in++;
	return n_in;	
}


extern int CHECK_REGULATION(Graph* graph, int regulator, int regulated)
{
	NodeAdj* Node;
	for(Node=(graph->array[regulator].head_out); Node!=NULL; Node=Node->next)
		if(Node->neighbor==regulated) return 1;
	return 0;
}

extern int IDENTIFY_SOLITAIRE(Graph* graph, int node)
{
	int i;
	int n_in = GETNin(graph, node);
	if(n_in>0)
	{
		int* neigh = GET_INNEIGH(graph, node);
		for(i=0; i<n_in; i++) if(neigh[i]!=node) return -1;
		return 1; // receives information only from itself.
	}
	else return 0; // do not receive information not even from itself.
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

////////////////////////////////
////// STACK OPERATIONS ///////
// To insert an element in the stack.
extern void push(STACK** top, int node)
{
	STACK *ptr;
	ptr = (STACK*)malloc(sizeof(STACK));
	ptr->node_ID = node;
	if(*top==NULL)
	{
		ptr->next = NULL;
		*top = ptr;
	}
	else
	{
		ptr->next = *top;
		*top = ptr;
	}
}

// To delete an element in the STACK.
extern int pop(STACK** top)
{
	int v;
	STACK *ptr;
	ptr = *top;
	if((*top)!=NULL)
	{
		v = (*top)->node_ID;
		*top = (*top)->next;
		free(ptr);
	}
	return v;
}

extern int STACKSIZE(STACK* top)
{
	int size = 0;
	STACK* current;
	for(current=top; current!=NULL; current=current->next) size++;
	return size;
}
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
////////////////// DOUBLE LINKED LIST DATA STRUCTURE //////////////////
extern int doublycheck_element(NODELIST* head, int element)
{
	NODELIST* nodelist = head;
	for(nodelist=head; nodelist!=NULL; nodelist=nodelist->next)
		if(nodelist->data==element) return 1;
	return 0;
}

extern int edgesfromSet(Graph* graph, int node, BLOCK* Set, int type)
{
	int check;	
	NodeAdj* Node;

	int n_in = 0;
	for(Node=graph->array[node].head_in; Node!=NULL; Node=Node->next)
	{
		if(Node->type_link!=type) continue;		
		check = doublycheck_element(Set->head, Node->neighbor);
		if(check==1) n_in++;
	}
	return n_in;	
}

extern void push_doublylist(NODELIST** head, int insertion)
{
    NODELIST* new_node = (NODELIST*)malloc(sizeof(NODELIST));
    
    new_node->data = insertion;
    new_node->next = (*head);
    new_node->prev = NULL;

    if((*head)!=NULL) (*head)->prev = new_node;
    (*head) = new_node;
}

extern void add_to_block(BLOCK** block, int node)
{
	BLOCK* temp = *block;	
	push_doublylist(&(temp->head), node);
	(*block)->size++;
}

extern void Copy_NodeList(NODELIST** dest, NODELIST* source)
{
	NODELIST* aux = source;
	while(aux)
	{
		push_doublylist(dest, aux->data);
		aux = aux->next;
	}
}

extern void insertAfter(NODELIST* prev_node, int insertion)
{
    if (prev_node == NULL) return; 

    /* 2. allocate new node */
    NODELIST* new_node = (NODELIST*)malloc(sizeof(NODELIST)); 
  
    new_node->data = insertion; 
    new_node->next = prev_node->next; 
    prev_node->next = new_node; 
    new_node->prev = prev_node; 
  
    /* Change previous of new_node's next node */
    if (new_node->next != NULL) 
        new_node->next->prev = new_node; 
}

extern void append(NODELIST** head, int new_data)
{
    /* Allocate node */
    NODELIST* new_node = (NODELIST*)malloc(sizeof(NODELIST)); 
  
    NODELIST* last = *head; /* used in step 5*/
  
    new_node->data = new_data;
    new_node->next = NULL; 
  
    /* If the Linked List is empty, then make the new 
        node as head */
    if (*head == NULL) 
    { 
        new_node->prev = NULL; 
        *head = new_node; 
        return; 
    } 
  
    /* Else traverse till the last node */
    while (last->next != NULL)  last = last->next; 
  
    last->next = new_node; 
    new_node->prev = last; 
  
    return; 
}

extern void KOSAJARU(NODELIST** nodes_in_scc, int root, Graph* graph)
{
	int i, j, v;
	int N = graph->size;
	int* visited = (int*)malloc(N*sizeof(int));
	int* ivisited = (int*)malloc(N*sizeof(int));
	for(i=0; i<N; i++) { visited[i] = 0; ivisited[i] = 0; }

	// DFS process in the network starting from the root node.
	STACK* nodetocheck = NULL;
	push(&nodetocheck, root);
	while(nodetocheck)
	{
		v = pop(&nodetocheck);
		if(visited[v]==0)
		{
			visited[v] = 1;
			int n = GETNout(graph, v);
			int* neigh_out = GET_OUTNEIGH(graph, v);
			for(j=0; j<n; j++) push(&nodetocheck, neigh_out[j]);
			free(neigh_out);
		}
	}
	// DFS process in the transpose network from the root node.
	push(&nodetocheck, root);
	while(nodetocheck)
	{
		v = pop(&nodetocheck);
		if(visited[v]==1)
		{
			int check = doublycheck_element(*nodes_in_scc, v);
			if(check==0) push_doublylist(nodes_in_scc, v); // it is in the root SCC.
		}
		if(ivisited[v]==0)
		{
			ivisited[v] = 1;
			int n = GETNin(graph, v);
			int* neigh_in = GET_INNEIGH(graph, v);
			for(j=0; j<n; j++) push(&nodetocheck, neigh_in[j]);
			free(neigh_in);
		}
	}

	// 'nodes_in_scc' contains all the nodes that are in the same SCC than root.
}

int EqualBlocks(BLOCK* block1, BLOCK* block2)
{
	if(block1->size!=block2->size) return -1;
	
	NODELIST* head1 = block1->head;
	NODELIST* head2 = block2->head;	
	while(head1)
	{
		if(head1->data!=head2->data) return -1;
		head1 = head1->next;
		head2 = head2->next;
	}
	if(head1==NULL && head2!=NULL) return -1;
	return 1;
}

void deleteNode(NODELIST** head_ref, NODELIST* del) 
{ 
    /* base case */
    if (*head_ref == NULL || del == NULL) 
        return; 
  
    /* If node to be deleted is head node */
    if (*head_ref == del) 
        *head_ref = del->next; 
  
    /* Change next only if node to be deleted is NOT the last node */
    if (del->next != NULL) 
        del->next->prev = del->prev; 
  
    /* Change prev only if node to be deleted is NOT the first node */
    if (del->prev != NULL) 
        del->prev->next = del->next; 
  
    /* Finally, free the memory occupied by del*/
    free(del); 
    return; 
}

void DeletePart(PART** head_ref, PART* del) 
{ 
	BLOCK* tempblock = (del)->block;    
	if(tempblock!=NULL)	// First delete all the nodes of the block.
	{
		NODELIST* temp = tempblock->head;
		while(temp) { deleteNode(&temp, temp); }
	}
  
    if (*head_ref == NULL || del == NULL) 
        return; 
  
    /* If node to be deleted is head node */
    if (*head_ref == del) 
        *head_ref = (del)->next; 
  
    /* Change next only if node to be deleted is NOT the last node */
    if ((del)->next != NULL) 
        (del)->next->prev = (del)->prev; 
  
    /* Change prev only if node to be deleted is NOT the first node */
    if ((del)->prev != NULL) 
        (del)->prev->next = (del)->next; 
  
    /* Finally, free the memory occupied by del*/
    free(del); 
    return; 
}

void freePart(PART** head_ref)
{
	PART* AUX;
	while(head_ref)
	{
		AUX = *head_ref;
		*head_ref = (*head_ref)->next;
		DeletePart(&AUX, AUX);
		if(AUX==NULL) printf("okay\n");
	}
}

extern void push_block(PART** head, BLOCK* insertion)
{
    PART* new_node = (PART*)malloc(sizeof(PART));
    
    new_node->block = insertion;
    new_node->next = (*head);
    new_node->prev = NULL;
    if((*head)!=NULL) (*head)->prev = new_node;
    (*head) = new_node;
}

extern void DeleteBlockInPartition(PART** part, BLOCK* P)
{
	int equal;	
	PART* current_part = *part;
	while(current_part)
	{
		equal = EqualBlocks(P, current_part->block);
		if(equal==1) break;
		current_part = current_part->next;
	}
	
	if(current_part!=NULL) DeletePart(part, current_part);
}

int GetPartitionSize(PART* part)
{
	int i = 0;
	PART* temp = part;
	while(temp)
	{
		i++;
		temp = temp->next;
	}
	return i;
}

int GetBlockSize(BLOCK* P)
{
	return P->size;
}

int GetListSize(NODELIST* list)
{
	int size = 0;
	NODELIST* nodelist;
	for(nodelist=list; nodelist!=NULL; nodelist=nodelist->next) size++;
	return size;
}

int GetFiberNumber(PART* part)
{
	int nfibers = 0;
	PART* current_part;
	for(current_part=part; current_part!=NULL; current_part=current_part->next)
		if(current_part->block->size>1) nfibers++;
	return nfibers;
}

int GetFiberNumber1(PART* part, PART* nullpart)
{
	int nfibers = 0;
	PART* current_part;
	NODELIST* nodelist;
	for(current_part=part; current_part!=NULL; current_part=current_part->next)
	{
		int ntemp = 0;
		for(nodelist=current_part->block->head; nodelist!=NULL; nodelist=nodelist->next) ntemp++;
		if(ntemp>1) nfibers++;
	}
	return nfibers;
}


/////////////////////////////////////////////////////////////

///////// IMPLEMENTATION OF QUEUE DATA STRUCTURE ////////
/////////////////////////////////////////////////////////
BLOCK* COPYBLOCK(BLOCK* P)
{
	
	BLOCK* new = (BLOCK*)malloc(sizeof(BLOCK));

	int size = P->size;
	int index = P->index;
	new->size = size;
	new->index = index;
	new->head = NULL;

	NODELIST* current_node = P->head;
	while(current_node)
	{
		push_doublylist(&(new)->head, current_node->data);
		current_node = current_node->next;
	}
	return new;
}

extern void enqueue_block(QBLOCK** head, QBLOCK** tail, BLOCK* P)
{
	BLOCK* current_block = COPYBLOCK(P);
	QBLOCK* new_element = (QBLOCK*)malloc(sizeof(QBLOCK));
	
	new_element->next = NULL;
	new_element->block = current_block;
	
	if(*head==NULL)
	{
		(*head) = new_element;
		(*tail) = new_element;
		return;
	}
	(*tail)->next = new_element;
	(*tail) = new_element;
}

extern BLOCK* dequeue_block(QBLOCK** head, QBLOCK** tail)
{
	BLOCK* block = NULL;
	
	if(*head==NULL) return NULL;
	QBLOCK* top = *head;
	(*head) = (*head)->next;
	block = top->block;
	free(top);
	if(*head==NULL) (*tail) = NULL;
	return block;
}

BLOCK* peek_block(QBLOCK** head)
{
	QBLOCK* top = *head;
	if(top!=NULL) return top->block;
	else return NULL;
}

void printQueueSize(QBLOCK* head)
{
	int i = 0;	
	QBLOCK* temp = head;
	while(temp)
	{
		i++;
		temp = temp->next;
	}
	printf("%d\n", i);
}

void printQueue(QBLOCK* head)
{
	if(head==NULL) printf("EMPTY\n");	
	QBLOCK* temp = head;
	while(temp)
	{
		printf("Block %d with size %d: ", temp->block->index, temp->block->size);
		printBlock(temp->block);
		temp = temp->next;
	}
}
//////////////////////////////////////////////////////////////////////////////////////////

#endif
