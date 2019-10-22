/*	This module contains a compilation of utility functions of some
	type of situations. Specially, here we define the functions necessary
	to construct directed and undirected networks. 

	Author: Higor da S. Monteiro - Universidade Federal do Ceará
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include "structforfiber.h"

////////////////////////////////////////////////////////////////////
/// Structures to create and define a directed, unweighted graph ///
////////////////////////////////////////////////////////////////////

// Define a graph containing 'size' nodes
struct Graph
{
	// each node has its adjancency list
	int size;
	struct adjList* array;
};
typedef struct Graph Graph;

struct adjList
{
	struct NodeAdj* head_in;
    struct NodeAdj* head_out;
};
typedef struct adjList AdjList;

// Structure to define the adjacency list of a node
struct NodeAdj
{
	int neighbor;
    int regulator;
	struct NodeAdj* next;
    //struct NodeAdj* next_out;
};
typedef struct NodeAdj NodeAdj;

NodeAdj* createNode(int neigh, int type_link)
{
	NodeAdj* newnode = (NodeAdj*)malloc(sizeof(NodeAdj));
	newnode->neighbor = neigh;
	newnode->next = NULL;
    newnode->regulator = type_link;
	return newnode;
}

extern Graph* createGraph(int N)
{
	Graph* graph = (Graph*)malloc(sizeof(Graph));
	graph->size = N;

	// Array of NodeAdj's
	graph->array = (AdjList*)malloc(N*sizeof(AdjList));

	int j;
	for(j=0; j<N; j++) { graph->array[j].head_in = NULL; graph->array[j].head_out = NULL; }
	return graph;
}

void addEdges(int** edges, Graph* graph, int* regulator, int nE)
{
	int index, j, node1, node2, reg;
	
	for(j=0; j<nE; j++)
	{
		// 'node1' -> 'node2' directed link.
        node1 = edges[j][0];
		node2 = edges[j][1];
		reg = regulator[j];

		NodeAdj* newnode1 = createNode(node1, reg);
		NodeAdj* newnode2 = createNode(node2, reg);

		newnode2->next = graph->array[node1].head_out;
		newnode1->next = graph->array[node2].head_in;
		graph->array[node1].head_out = newnode2;
		graph->array[node2].head_in = newnode1;
	}
}

extern int** defineNetwork(int** edges, int* regulator, Graph* graph, char* filename)
{
	FILE *EDGE_FILE = fopen(filename, "r");
	if(EDGE_FILE==NULL) printf("ERROR in file reading");

	int i, j;
	char k[5];
	int r = 1;
	int nlink = 0; // number of links.
	while(r) // Calculates the number of lines in the file
	{
		r = fscanf(EDGE_FILE, "%d\t%d\t%s\n", &i, &j, &k);
		if(r==EOF) break;
		nlink++;
	}
	rewind(EDGE_FILE);

	edges = (int**)malloc(nlink*sizeof(int*));
    regulator = (int*)malloc(nlink*sizeof(int));
	for(j=0; j<nlink; j++)
	{
		edges[j] = (int*)malloc(2*sizeof(int));
		r = fscanf(EDGE_FILE, "%d\t%d\t%s\n", &edges[j][0], &edges[j][1], &k);
        if(strcmp("+", k)==0) regulator[j] = 0;
        else if(strcmp("-", k)==0) regulator[j] = 1;
        else if(strcmp("+-", k)==0) regulator[j] = 2;
        else regulator[j] = -1;
	}
	fclose(EDGE_FILE);
	
    // Defines the network structure.
	addEdges(edges, graph, regulator, nlink);
    
    return edges;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/////////////////////////// GRAPH OPERATIONS ////////////////////////////////

struct neigh_vector
{
	int neigh;
	int color;
};
typedef struct neigh_vector neigh_vector;

// To sort in ascendent order //
int cmpcolor(const void *a, const void *b)
{
    neigh_vector *a1 = (neigh_vector *)a;
    neigh_vector *a2 = (neigh_vector *)b;
    if ((*a1).color > (*a2).color)
        return 1;
    else if ((*a1).color < (*a2).color)
        return -1;
    else
        return 0;
}

extern neigh_vector* GET_inNEIGH(Graph* graph, int* nodecolor, int node, int N)
{
	int n_in = 0;	
	NodeAdj* NODE = graph->array[node].head_in;
	while(NODE)
	{
		n_in++;
		NODE = NODE->next;
	}
	neigh_vector* neighbors = (neigh_vector*)malloc((n_in)*sizeof(neigh_vector));

	int n_index = 0;
	NODE = graph->array[node].head_in;
	while(NODE)
	{
		neighbors[n_index].neigh = NODE->neighbor;
		neighbors[n_index].color = nodecolor[NODE->neighbor];
		n_index++;
		NODE = NODE->next;
	}
	qsort(neighbors, n_in, sizeof(neighbors[0]), cmpcolor);
	return neighbors;
}

extern neigh_vector* GET_outNEIGH(Graph* graph, int* nodecolor, int node, int N)
{
	int n_out = 0;	
	NodeAdj* NODE = graph->array[node].head_out;
	while(NODE)
	{
		n_out++;
		NODE = NODE->next;
	}

	neigh_vector* neighbors = (neigh_vector*)malloc(n_out*sizeof(neigh_vector));

	n_out = 0;
	NODE = graph->array[node].head_out;
	while(NODE)
	{
		neighbors[n_out].neigh = NODE->neighbor;
		neighbors[n_out].color = nodecolor[NODE->neighbor];
		n_out++;
		NODE = NODE->next;
	}
	qsort(neighbors, n_out, sizeof(neighbors[0]), cmpcolor);
	
	return neighbors;
}

extern int GETNin(Graph* graph, int node)
{
	int n_in = 0;	
	NodeAdj* NODE = graph->array[node].head_in;
	while(NODE)
	{
		n_in++;
		NODE = NODE->next;
	}
	return n_in;
}

extern void PrintInNeighbors(Graph* graph, int node, int N)
{
	NodeAdj* Inode = graph->array[node].head_in;
	printf("Node %d: ", node);
	while(Inode)
	{
		printf("%d ", Inode->neighbor);
		Inode = Inode->next;
	}
	printf("\n");
}

extern void PrintOutNeighbors(Graph* graph, int node, int N)
{
	NodeAdj* Inode = graph->array[node].head_out;
	printf("Node %d: ", node);
	while(Inode)
	{
		printf("%d ", Inode->neighbor);
		Inode = Inode->next;
	}
	printf("\n");
}
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
/////////////////// STACK DATA STRUCTURE ////////////////////
struct stack
{
	int extn;
	int node_ID;
	struct stack *next;
};
typedef struct stack Stack;
////////////////////////////////
////// STACK OPERATIONS ///////
// To insert an element in the stack.
extern Stack *push(Stack* top, int node)
{
	Stack *ptr;
	ptr = (Stack*)malloc(sizeof(Stack));
	ptr->node_ID = node;
	if(top==NULL)
	{
		ptr->next = NULL;
		top = ptr;
	}
	else
	{
		ptr->next = top;
		top = ptr;
	}
	return top;
}

// To delete an element in the stack.
extern Stack *pop(Stack* top)
{
	Stack *ptr;
	ptr = top;
	if(top!=NULL)
	{
		top = top->next;
		free(ptr);
	}
	return top;
}
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

//////////////////////// GENERAL USE ////////////////////////////

// Reads the number of effective voxels and total number of voxels of the subject map.
extern int nlines_file(char* filename, int ncolumns)
{
	FILE *MODFILE = fopen(filename, "r");
	if(MODFILE==NULL) printf("ERROR in file reading");

	char ch[5];
    int i, j, k, m, x, y, z;
	int ne = 0;
	int r = 1;
	while(r) // Calculates the number of lines in the file
	{
		if(ncolumns==3) r = fscanf(MODFILE, "%d\t%d\t%s\n", &i, &j, &ch);
		else if(ncolumns==6) r = fscanf(MODFILE, "%d\t%d\t%d\t%d\t%d\t%d\n", &i, &j, &k, &x, &y, &z);
		if(r==EOF) break;
		ne++;
	}
	fclose(MODFILE);
	return ne;
}

extern void printGraph(Graph* graph)
{
	int j;
	for(j=0; j<graph->size; j++)
	{
		printf("NODE %d\nin:", j);
        NodeAdj* SeeAux = graph->array[j].head_in;
        NodeAdj* SeeAux1 = graph->array[j].head_out;
		while(SeeAux)
		{
			printf("<-%d", SeeAux->neighbor);
			SeeAux = SeeAux->next; 
		}
		printf("\nout:");
        while(SeeAux1)
        {
            printf("->%d", SeeAux1->neighbor);
			SeeAux1 = SeeAux1->next;
        }
        printf("\n");
	}
}

extern int** arrint2d(int lines, int columns)
{
    int i, j;
    int** arr = (int**)malloc(lines*sizeof(int*));
    for(i=0; i<lines; i++)
    {
        arr[i] = (int*)malloc(columns*sizeof(int));
        for(j=0; j<columns; j++)
            arr[i][j] = 0;
    }
    return arr;
}

extern void free2d(int** arr, int rows, int columns)
{
	int i;
	for(i=0; i<rows; i++) free(arr[i]);
	free(arr);
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

////////////////// DOUBLE LINKED LIST DATA STRUCTURE //////////////////
struct DoublyLinkNode
{
    int data;
    struct DoublyLinkNode* prev;
    struct DoublyLinkNode* next;
};
typedef struct DoublyLinkNode DoublyLinkNode;

// The doubly linked list of the block containing the elements on it
struct BLOCK
{
    int size;
	int index;
    DoublyLinkNode* head;
};
typedef struct BLOCK BLOCK;

struct Partition
{
	BLOCK* block;
	struct Partition* prev;
	struct Partition* next;
};
typedef struct Partition PART;
///////////////////////////////////////////////////////////

extern int nblock(PART** head)
{
	if(*head==NULL) return 0;
	
	int i = 0;	
	PART* temp = *head;
	while(temp)
	{
		i++;
		temp = temp->next;
	}
	return i;
}

extern int doublycheck_element(DoublyLinkNode* head, int element)
{
	DoublyLinkNode* temp = head;
	while(temp)
	{
		if(temp->data==element) return 1;
		temp = temp->next;
	}
	return 0;
}

extern int intersection_edges(Graph* graph, int node, BLOCK* Set)
{
	int i, check;	
	NodeAdj* Node = graph->array[node].head_in;

	int n_in = 0;
	while(Node)
	{
		check = doublycheck_element(Set->head, Node->neighbor);
		if(check==1) n_in++;
		Node = Node->next;
	}
	return n_in;	
}

extern void push_doublylist(DoublyLinkNode** head, int insertion)
{
    DoublyLinkNode* new_node = (DoublyLinkNode*)malloc(sizeof(DoublyLinkNode));
    
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

extern void insertAfter(DoublyLinkNode* prev_node, int insertion)
{
    if (prev_node == NULL) return; 

    /* 2. allocate new node */
    DoublyLinkNode* new_node = (DoublyLinkNode*)malloc(sizeof(DoublyLinkNode)); 
  
    new_node->data = insertion; 
    new_node->next = prev_node->next; 
    prev_node->next = new_node; 
    new_node->prev = prev_node; 
  
    /* Change previous of new_node's next node */
    if (new_node->next != NULL) 
        new_node->next->prev = new_node; 
}

extern void append(DoublyLinkNode** head, int new_data)
{
    /* Allocate node */
    DoublyLinkNode* new_node = (DoublyLinkNode*)malloc(sizeof(DoublyLinkNode)); 
  
    DoublyLinkNode* last = *head; /* used in step 5*/
  
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

int EqualBlocks(BLOCK* block1, BLOCK* block2)
{
	if(block1->size!=block2->size) return -1;
	
	DoublyLinkNode* head1 = block1->head;
	DoublyLinkNode* head2 = block2->head;	
	while(head1)
	{
		if(head1->data!=head2->data) return -1;
		head1 = head1->next;
		head2 = head2->next;
	}
	if(head1==NULL && head2!=NULL) return -1;
	return 1;
}

void deleteNode(DoublyLinkNode** head_ref, DoublyLinkNode* del) 
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
		DoublyLinkNode* temp = tempblock->head;
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

void printBlock(BLOCK* P)
{
    DoublyLinkNode* List = P->head;
    while(List)
    {
        printf("%d ", List->data);
        List = List->next;
    }
	printf("\n");
}

void printBlockSize(BLOCK* P)
{
	printf("Size: %d\n", P->size);
}

int GetBlockSize(BLOCK* P)
{
	return P->size;
}

void printAllPartition(PART* head)
{
	if(head==NULL) printf("EMPTY\n");	
	PART* temp = head;
	BLOCK* aux_block = NULL;
	while(temp)
	{
		aux_block = temp->block;
		printf("Block %d with size %d: ", aux_block->index, aux_block->size);
		printBlock(aux_block);
		temp = temp->next;
	}
}

void printPartition(PART* P)
{
    PART* part = P;
	BLOCK* tempblock;
    while(part)
    {
		tempblock = part->block;       
		printf("%d ", tempblock->index);
        part = part->next;
    }
}

void printPartitionSize(PART* part)
{
	int i = 0;
	PART* temp = part;
	while(temp)
	{
		i++;
		temp = temp->next;
	}
	printf("\nPartition size: %d\n", i);
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
/////////////////////////////////////////////////////////////

///////// IMPLEMENTATION OF QUEUE DATA STRUCTURE ////////
/////////////////////////////////////////////////////////
struct QueueOfBlocks
{
    BLOCK* block;
    struct QueueOfBlocks* next;
};
typedef struct QueueOfBlocks QBLOCK;

void enqueue_block(QBLOCK** head, QBLOCK** tail, BLOCK* P)
{
	QBLOCK* new_element = (QBLOCK*)malloc(sizeof(QBLOCK));
	
	new_element->block = P;
	new_element->next = NULL;	
	
	if(*tail!=NULL) (*tail)->next = new_element;
	if(*head==NULL) (*head) = new_element;
	(*tail) = new_element;
}

void dequeue_block(QBLOCK** head, QBLOCK** tail)
{
	QBLOCK* top = *head;
	QBLOCK* next_qblock = NULL;	
	if(*head!=*tail) next_qblock = (*head)->next;

	free(top);
	if(*head==*tail) *tail = next_qblock;
	*head = next_qblock; 
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


