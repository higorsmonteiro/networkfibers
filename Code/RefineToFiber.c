/*  The code given here has the aim to determine the minimal balanced color classes to identify
    the fibration building blocks of a directed network. The theory concerning the graph fibration
    morphism and its application on biological networks are detailed, respectively, mainly in the
    paper of Boldi and Vigna (2001) and in the paper of Morone et. al. (2019, to be published). Concerning
    the algorithm reproduced in this code, one should refer to the work of Paige and Tarjan (1987).

    The resulted fibers represent, for each fiber, the set of nodes that are isomorphic under a graph
    fibration morphism, i. e., all nodes belonging to the same fiber has the same information input-tree.

    -----------------------------------------------------------------------------------------------------

    The code has one command line argument, requiring the common identifier for two files: one single-line file 
	containing a integer number representing the number of nodes in the network, and a 'ARGedgelist.dat' file containing 
	all the directed links between nodes (3 columns -> Pointing Node/ Pointed Node/ Type of regulation).

    The result is stored in two arrays: 'nodecolor' and 'edgecolor'. Each array store the color of each
    component (node or link) and, thus, gives all the information needed, together with the files, for the
    base graph construction.

    Author: Higor da S. Monteiro - Departament of Physics/Universidade Federal do Ceará (UFC)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_eigen.h>
#include "utilsforfiber.h"   // Personal module containing helpful functions to handle with network data and other common operations.

void UPGRADE_PARTITION(PART** new_blocks, PART** old_blocks, PART** partition)
{
	BLOCK* temp_block = NULL;	
	PART* current_old = *old_blocks;
	PART* current_new = *new_blocks;
	PART* current_part;
	// Erase in 'partition' all the blocks contained in 'old_blocks'.
	while(current_old)
	{
		DeleteBlockInPartition(partition, current_old->block);
		current_old = current_old->next;
	}
	// push all the new blocks into the partition.	
	while(current_new)
	{
		push_block(partition, current_new->block);
		current_new = current_new->next;
	}
}

void EnqueueAllNotLargest(PART** new_blocks, QBLOCK** qhead, QBLOCK** qtail)
{
	int max_index = -1;
	int size = -1;
	PART* temp = *new_blocks;
	BLOCK* tempblock;	
	// Get the index of the largest block.	
	while(temp)
	{
		tempblock = temp->block;
		if(tempblock->size>size) { size = tempblock->size; max_index = tempblock->index; }
		temp = temp->next;	
	}

	temp = *new_blocks;
	while(temp)
	{
		tempblock = temp->block;
		if(tempblock->index!=max_index) enqueue_block(qhead, qtail, tempblock);
		temp = temp->next;
	}
}

/*	given a node and its index, adds it on a block that has the same index
	or create a new block in case there isn't any block with the given index. */
void Push_On_Block(int node, int nodeindex, PART** subpart2)
{
/*	First thing, we need to loop over the blocks of the given partition
	and check if there is a block with the same index. */
	PART* current_part = *subpart2;
	DoublyLinkNode* nodelist = NULL;
	
	while(current_part)
	{
		if(current_part->block->index==nodeindex)
		{
			push_doublylist(&(current_part->block->head), node);
			current_part->block->size += 1;
			return;
		}
		current_part = current_part->next;
	}

/*	If there isn't a block with the same index, then we create a new one. */
	BLOCK* new_block = (BLOCK*)malloc(sizeof(BLOCK));
	new_block->index = nodeindex;
	new_block->size = 0;
	new_block->head = NULL;
	push_doublylist(&(new_block->head), node);
	push_block(subpart2, new_block);
	return;
}

/*	Given the 'Set' block, now we select all the blocks in the 'partition' that have at 
	least one connection coming from 'Set'. We store these blocks into 'subpart'. */
void GET_NONSTABLE_BLOCKS(PART** partition, PART** subpart, Graph* graph, BLOCK* Set)
{
	BLOCK* tempblock;
	PART* temp_part = *partition;	
	DoublyLinkNode* tempdoublylist;	
	while(temp_part)
	{
		tempblock = temp_part->block;		// Get current block of the original partition.
		tempdoublylist = tempblock->head;	// Get the list of node of the current block.
		while(tempdoublylist)
		{
			int node = tempdoublylist->data;
			int n_in = intersection_edges(graph, node, Set);	// Get number of incoming links from 'Set' to 'node'.
			if(n_in>0) { push_block(subpart, tempblock); break; }
			tempdoublylist = tempdoublylist->next;
		}
		temp_part = temp_part->next;
	}
}

void BLOCKS_PARTITIONING(PART** subpart1, PART** subpart2, Graph* graph, BLOCK* Set)
{
	int node, n_in;	
	BLOCK* temp_block = NULL;
	PART* temp_part = *subpart1;	
	DoublyLinkNode* nodelist = NULL;		
	
	while(temp_part) // for each block to be splitted.
	{
		temp_block = temp_part->block;		
		nodelist = temp_block->head;	
		while(nodelist)	// for each node in the current block.
		{		
			node = nodelist->data;			
			n_in = intersection_edges(graph, node, Set);
			Push_On_Block(node, n_in, subpart2);	// *error* put 'node' in its proper new block.
			nodelist = nodelist->next;
		}
		temp_part = temp_part->next;
	}
}	

void S_SPLIT(PART** partition, BLOCK* Set, Graph* graph, QBLOCK** qhead, QBLOCK** qtail)
{	
	PART* subpart1 = NULL;
	PART* subpart2 = NULL;
/*	Given the 'Set' block, now we select all the blocks in the partition that have at 
	least one connection coming from 'Set'. */
	GET_NONSTABLE_BLOCKS(partition, &subpart1, graph, Set);
/* 'subpart1' contains all the blocks that have nonzero pointed nodes from 'Set'. */
	
	BLOCKS_PARTITIONING(&subpart1, &subpart2, graph, Set);

/*	After the process above, 'subpart2' contains the resulted splitted blocks. And now we
	delele all 'subpart1' blocks from partition and then we push all those 'subpart2' 
	blocks to the 'partition'. Before deleting 'subpart2', we enqueue all the blocks
	in the queue of refining blocks, except the largest one. */
	//printAllPartition(*partition);	
	UPGRADE_PARTITION(&subpart2, &subpart1, partition);
	//printAllPartition(*partition);	
	EnqueueAllNotLargest(&subpart2, qhead, qtail);
	
	// obs: 'subpart1' and 'subpart2' are erased by scope.
}

void ENQUEUE_ALONE_NODES(PART** null_partition, QBLOCK** qhead, QBLOCK** qtail)
{
	PART* current_part = *null_partition;
	while(current_part)
	{
		enqueue_block(qhead, qtail, current_part->block);
		current_part = current_part->next;
	}
}

void InitiateBlock(PART** Null_Part, int node)
{
	BLOCK* new_block = (BLOCK*)malloc(sizeof(BLOCK));
	new_block->index = -1;
	new_block->size = 0;
	new_block->head = NULL;
	add_to_block(&new_block, node);
	push_block(Null_Part, new_block);
}

void PREPROCESSING(BLOCK** P, BLOCK** NonP, PART** Null_Part, Graph* graph, int N)
{
	(*P)->size = 0;
	(*NonP)->size = 0;
	(*P)->index = 0;
	(*NonP)->index = -1;
	(*P)->head = NULL;
	(*NonP)->head = NULL;    

	int n_in, i;
    for(i=0; i<N; i++)
    {
		n_in = GETNin(graph, i);
        if(n_in>0) add_to_block(P, i);
		else add_to_block(NonP, i);
    }

	DoublyLinkNode* nodelist = (*NonP)->head;
	while(nodelist)
	{
		InitiateBlock(Null_Part, nodelist->data);
		nodelist = nodelist->next;
	}
}

int main(int argv, char** argc) 
{ 
	int N, M;                         	// Number of nodes and edges of the network.
	char netsize[100] = "../Data/";     // File containing (one line) the number of nodes in the network.
	char net_edges[100] = "../Data/";   // File containing all the directed links in the network.
	strcat(netsize, argc[1]);
	strcat(netsize, "Ngenes.dat");
	strcat(net_edges, argc[1]);
	strcat(net_edges, "edgelist.dat");
                                                          
    // Defines the size of the network.
	FILE* UTIL = fopen(netsize, "r");
	if(UTIL==NULL) printf("ERROR IN FILE READING\n");
	fscanf(UTIL, "%d\n", &N);
	fclose(UTIL);
	M = nlines_file(net_edges, 3);	// Number of edges.
    ///////////////////////////////////////////////////

    // Creates the graph structure for N nodes
	Graph* graph = createGraph(N);

    // Properly defines the network structure with the given 'edgelist.dat' file.
	int** edges;
	int* regulator;
	edges = defineNetwork(edges, regulator, graph, net_edges);
	///////////////////////////////////////////////////////////////////////////////////////

	/////////////////////// COARSEST REFINEMENT PARTITION ALGORITHM ////////////////////////
	int i, j, k;  		
           
	// Put aside from the algorithm (index -1) all nodes that don't receive any information.    
	BLOCK* P = (BLOCK*)malloc(sizeof(BLOCK));
	BLOCK* NonP = (BLOCK*)malloc(sizeof(BLOCK));
	
	PART* null_partition = NULL;
	PREPROCESSING(&P, &NonP, &null_partition, graph, N);
	// Define the initial partition as one block containing all operating nodes. 
		
	PART* partition = NULL;    
	push_block(&partition, P);	// Push initial block to the partition.
	
	// Initialize the queue of blocks with the initial block above.
	QBLOCK* qhead = NULL;
	QBLOCK* qtail = NULL;
	enqueue_block(&qhead, &qtail, P);
	ENQUEUE_ALONE_NODES(&null_partition, &qhead, &qtail);

	// Until L is empty, we procedure the splitting process.	
	int time = 0;
	printf("Time %d\n", time);
		
	BLOCK* CurrentSet;	
	//printAllPartition(partition);	
	while(qhead)
	{
		time++;		
		CurrentSet = peek_block(&qhead);
		dequeue_block(&qhead, &qtail);		
		S_SPLIT(&partition, CurrentSet, graph, &qhead, &qtail);
		printf("Time %d\n", time);
		printAllPartition(partition);
	}

	int size = GetPartitionSize(partition);
	int presize = GetPartitionSize(null_partition);
	printf("Number of fiber: %d\n", size+presize);
	//printAllPartition(partition);
    ////////////////////////////////////////////////////////////////////////////////////

}
