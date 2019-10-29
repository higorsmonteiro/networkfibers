#ifndef FSTRUCTS_H
#define FSTRUCTS_H

#include <stdio.h>

/////////////////////////////////////////////////////////////////////
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
	char gene_name[60];	
	struct NodeAdj* head_in;
    struct NodeAdj* head_out;
};
typedef struct adjList AdjList;

// Structure to define the adjacency list of a node
struct NodeAdj
{
	int neighbor;
    int type_link;
	struct NodeAdj* next;
};
typedef struct NodeAdj NodeAdj;
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

////////////////////// STACK DATA STRUCTURE ///////////////////////
struct stack
{
	int extn;
	int node_ID;
	struct stack *next;
};
typedef struct stack Stack;
/////////////////////////////////////////////////////////////////////

////////////////// DOUBLE LINKED LIST DATA STRUCTURE //////////////////
struct NODELIST
{
    int data;
    struct NODELIST* prev;
    struct NODELIST* next;
};
typedef struct NODELIST NODELIST;

struct BLOCK
{
    int size;
	int index;
	int pos;
	int neg;
	int dual;
    NODELIST* head;
};
typedef struct BLOCK BLOCK;

struct Partition
{
	BLOCK* block;
	struct Partition* prev;
	struct Partition* next;
	int number_regulators;	
	double fundamental_number;
	NODELIST* regulators;
};
typedef struct Partition PART;
////////////////////////////////////////////////////////////////////////

struct QueueOfBlocks
{
    BLOCK* block;
    struct QueueOfBlocks* next;
};
typedef struct QueueOfBlocks QBLOCK;
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

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

void printBlock(BLOCK* P)
{
    NODELIST* List = P->head;
    while(List)
    {
        printf("%d ", List->data);
        List = List->next;
    }
	printf("\n");
}

void printBlockGene(BLOCK* P, Graph* graph)
{
    NODELIST* List = P->head;
	for(List=P->head; List!=NULL; List=List->next)
		printf("%s, ", graph->array[List->data].gene_name);
	printf("\n");
}

void printBlockSize(BLOCK* P)
{
	printf("Size: %d\n", P->size);
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

void printGenesPartition(PART* head, Graph* graph)
{
	if(head==NULL) printf("EMPTY\n");	
	PART* temp = head;
	BLOCK* aux_block = NULL;
	while(temp)
	{
		aux_block = temp->block;
		printf("Block %d with size %d: ", aux_block->index, aux_block->size);
		printBlockGene(aux_block, graph);
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
	printf("\n");
}

void printPartitionSize(PART* part)
{
	int i = 0;
	PART* temp;
	for(temp=part; temp!=NULL; temp=temp->next) i++;
	printf("Partition size: %d\n", i);
}

void printPartitionReg(PART* part)
{
	int valid_fibers = 0;
	PART* temp;
	for(temp=part; temp!=NULL; temp=temp->next)
		if(temp->block->size>1) printf("Number of external regulators of block %d: %d\n", temp->block->index, temp->number_regulators);
}

void ShowMainInfo(PART* partition)
{
	PART* current_part;
	for(current_part=partition; current_part!=NULL; current_part=current_part->next)
	{
		printf("Fiber %d: Size %d - Fundamental Class %lf - Subclass %d\n", current_part->block->index, current_part->block->size, current_part->fundamental_number, current_part->number_regulators);
	}
}
#endif
