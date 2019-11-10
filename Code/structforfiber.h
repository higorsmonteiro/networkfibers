#ifndef FSTRUCTS_H
#define FSTRUCTS_H

#include <stdio.h>
#include <stdlib.h>
#include "utilsforfiber.h"

/////////////////////////////////////////////////////////////////////
/// Structures to create and define a directed, unweighted graph ///
////////////////////////////////////////////////////////////////////
// Define a graph containing 'size' nodes
struct Graph
{
	// each node has its adjancency list
	int size;
	int num_component;
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
	int node_ID;
	struct stack *next;
};
typedef struct stack STACK;
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

struct storetype
{
	int node;
	int type;
};
typedef struct storetype STORETYPE;
/////////////////////////////////////////////////////////////////////////

struct QueueOfBlocks
{
    BLOCK* block;
    struct QueueOfBlocks* next;
};
typedef struct QueueOfBlocks QBLOCK;
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int cmp(const void *a, const void *b)
{
    STORETYPE *a1 = (STORETYPE *)a;
    STORETYPE *a2 = (STORETYPE *)b;
    if ((*a1).node > (*a2).node)
        return -1;
    else if ((*a1).node < (*a2).node)
        return 1;
    else
        return 0;
}

STORETYPE* GETIN_ADJTYPE(Graph* graph, int node)
{
	int n_in = 0;	
	NodeAdj* NODE;
	for(NODE=graph->array[node].head_in; NODE!=NULL; NODE=NODE->next) n_in++;

	int n_index = 0;
	STORETYPE* in_neighbors = (STORETYPE*)malloc(n_in*sizeof(STORETYPE));
	for(NODE=graph->array[node].head_in; NODE!=NULL; NODE=NODE->next)
	{
		in_neighbors[n_index].node = NODE->neighbor;
		in_neighbors[n_index].type = NODE->type_link;
		n_index++;
	}
	
	qsort(in_neighbors, n_in, sizeof(in_neighbors[0]), cmp);
	return in_neighbors; 
}

////////////////// VISUALIZATION UTILITIES ///////////////////
extern void PrintInNeighbors(Graph* graph, int node)
{
	NodeAdj* Inode = graph->array[node].head_in;
	printf("Node %d receives from: ", node);
	while(Inode)
	{
		printf("%d(type:%d) ", Inode->neighbor, Inode->type_link);
		Inode = Inode->next;
	}
	printf("\n");
}

extern void PrintOutNeighbors(Graph* graph, int node)
{
	NodeAdj* Inode = graph->array[node].head_out;
	printf("Node %d points to: ", node);
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
    NodeAdj* inode;
	for(j=0; j<graph->size; j++)
	{
		printf("NODE %d\nin:", j);
		for(inode=graph->array[j].head_in; inode!=NULL; inode=inode->next)
			printf("%d ", inode->neighbor);
		printf("\nout:");
        for(inode=graph->array[j].head_out; inode!=NULL; inode=inode->next)
			printf("%d ", inode->neighbor);
        printf("\n");
	}
}

extern void printGraphInFibers(Graph* graph, PART* partition, int* nodefiber)
{
	int i, n;
	STORETYPE* store;
	PART* current_part;
	NODELIST* nodelist;
	for(current_part=partition; current_part!=NULL; current_part=current_part->next)
	{
		if(current_part->block->size==1) continue;
		printf("FIBER %d WITH SIZE %d -> | %.4lf, %d >:\n", current_part->block->index, current_part->block->size, current_part->fundamental_number, current_part->number_regulators);
		for(nodelist=current_part->block->head; nodelist!=NULL; nodelist=nodelist->next)
		{
			n = GETNin(graph, nodelist->data);
			store = GETIN_ADJTYPE(graph, nodelist->data);
			printf("NODE %d receives from fiber(node,type): ", nodelist->data);
			for(i=0; i<n; i++) printf("%d(%d,%d) ", nodefiber[store[i].node], store[i].node, store[i].type);
			printf("\n");
		}
		printf("\n");
	}
}

extern void printGeneGraphInFibers(Graph* graph, PART* partition, int* nodefiber)
{
	int i, n;
	NodeAdj* Node;
	STORETYPE* temp;
	PART* current_part;
	NODELIST* nodelist;
	for(current_part=partition; current_part!=NULL; current_part=current_part->next)
	{
		if(current_part->block->size==1) continue;
		printf("FIBER %d WITH SIZE %d -> | %.4lf, %d >:\n", current_part->block->index, current_part->block->size, current_part->fundamental_number, current_part->number_regulators);
		for(nodelist=current_part->block->head; nodelist!=NULL; nodelist=nodelist->next)
		{
			printf("NODE %s receives from fiber(node,type): ", graph->array[nodelist->data].gene_name);
			n = GETNin(graph, nodelist->data);
			temp = GETIN_ADJTYPE(graph, nodelist->data);
			for(i=0; i<n; i++)
			{
				if(temp[i].type==0)
					printf("%d(%s,positive) ", nodefiber[temp[i].node], graph->array[temp[i].node].gene_name);
				else if(temp[i].type==1)
					printf("%d(%s,negative) ", nodefiber[temp[i].node], graph->array[temp[i].node].gene_name);
				else printf("%d(%s,dual) ", nodefiber[temp[i].node], graph->array[temp[i].node].gene_name);
			}
			printf("\n");
		}
		printf("\n");
	}
}
// ################################################# //

void printBlock(BLOCK* P)
{
    NODELIST* list;
	for(list=P->head; list!=NULL; list=list->next)	printf("%d ", list->data);
	printf("\n");
}

void printBlockGene(BLOCK* P, Graph* graph)
{
    NODELIST* List;
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

void printPartitionSize(PART* part)
{
	int i = 0;
	PART* temp;
	for(temp=part; temp!=NULL; temp=temp->next) i++;
	printf("Partition size: %d\n", i);
}

void printList(NODELIST* list)
{
	NODELIST* current;
	for(current=list; current!=NULL; current=current->next)
		printf("%d ", current->data);
	printf("\n");
}

void ShowClassification(PART* partition, int mode)
{
	PART* current_part;
	for(current_part=partition; current_part!=NULL; current_part=current_part->next)
	{
		if(current_part->block->size>1)
		{
			printf("Fiber %d: Size %d - n %lf - l %d\n", current_part->block->index, current_part->block->size, current_part->fundamental_number, current_part->number_regulators);
			//printBlock(current_part->block);
		}
		else if(mode==1)
		{
			printf("Fiber %d: Size %d - n %lf - l %d\n", current_part->block->index, current_part->block->size, current_part->fundamental_number, current_part->number_regulators);
		}
	}
}

void ShowClassification1(PART* partition, int mode)
{
	PART* current_part;
	for(current_part=partition; current_part!=NULL; current_part=current_part->next)
	{
		if(current_part->block->size>1)
		{
			printf("%d,%d,%lf,%d\n", current_part->block->index, current_part->block->size, current_part->fundamental_number, current_part->number_regulators);
			//printBlock(current_part->block);
		}
		else if(mode==1)
		{
			printf("Fiber %d: Size %d - n %lf - l %d\n", current_part->block->index, current_part->block->size, current_part->fundamental_number, current_part->number_regulators);
		}
	}
}

void ShowInfo(PART* partition, int mode)
{
	PART* current_part;
	for(current_part=partition; current_part!=NULL; current_part=current_part->next)
	{
		if(current_part->block->size>1)
		{
			printf("Fiber %d: Size %d - n %lf - l %d\n", current_part->block->index, current_part->block->size, current_part->fundamental_number, current_part->number_regulators);
			printBlock(current_part->block);
		}
		else if(mode==1)
		{
			printf("Fiber %d: Size %d - n %lf - l %d\n", current_part->block->index, current_part->block->size, current_part->fundamental_number, current_part->number_regulators);
			printBlock(current_part->block);
		}
	}
}



#endif
