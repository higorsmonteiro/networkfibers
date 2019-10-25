/*	In this module, I define all the functions that are, direct or indirectly, useful
	for the main fibration routines needed to a correct implementation for the partition
	refinement algorithm. Also, here I define my own functions to construct undirected 
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

extern Graph* createGraph(int N)
{
	Graph* graph = (Graph*)malloc(sizeof(Graph));
	graph->size = N;

	// Array of NodeAdj's //
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

extern int** defineNetwork(int** edges, Graph* graph, char* filename)
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
	
    // Defines the network structure.
	addEdges(edges, graph, regulator, nlink);
    
    return edges;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/////////////////////////// GRAPH OPERATIONS ////////////////////////////////
extern int* GET_INNEIGH(Graph* graph, int node)
{
	int n_in = 0;	
	NodeAdj* NODE;
	for(NODE=graph->array[node].head_in; NODE!=NULL; NODE=NODE->next) n_in++;

	int n_index = 0;
	int* in_neighbors = (int*)malloc(n_in*sizeof(int));
	for(NODE=graph->array[node].head_in; NODE!=NULL; NODE=NODE->next)
		in_neighbors[n_index++] = NODE->neighbor;
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
	return in_neighbors; 
}

extern int GETNin(Graph* graph, int node)
{
	int n_in = 0;	
	NodeAdj* NODE;
	for(NODE=graph->array[node].head_in; NODE!=NULL; NODE=NODE->next) n_in++;
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

extern int CHECKLINK(Graph* graph, int pointing_node, int pointed_node)
{
	NodeAdj* Node;
	for(Node=graph->array[pointed_node].head_out; Node!=NULL; Node=Node->next)
		if(Node->neighbor==pointed_node) return 1;
	return 0;
}
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

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