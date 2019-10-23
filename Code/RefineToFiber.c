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
#define BIGGERFLAG -96

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
	new_block->size = 1;
	new_block->head = NULL;
	push_doublylist(&(new_block->head), node);
	push_block(subpart2, new_block);
	return;
}

/*	Given the 'Set' block, now we select all the blocks in the 'partition' that have at 
	least one connection coming from 'Set'. We store these blocks into 'subpart'. */
void GET_NONSTABLE_BLOCKS(PART** partition, PART** subpart, int* n_fromSet, Graph* graph, BLOCK* Set)
{
	int check = 0;
	int blocksize, node, nextnode;
	PART* current_part = *partition;	
	DoublyLinkNode* nodelist;	
	while(current_part)
	{
		blocksize = current_part->block->size;
		nodelist = current_part->block->head;	// Get the list of node of the current block.
		while(nodelist)
		{
			if(nodelist->next!=NULL)
			{
				node = nodelist->data;
				nextnode = nodelist->next->data;
				if(n_fromSet[node]!=n_fromSet[nextnode]) { push_block(subpart, current_part->block); break; }
			}
			//if(n_fromSet[nodelist->data]>0) { push_block(subpart, current_part->block); break; }
			nodelist = nodelist->next;
		}
		current_part = current_part->next;
	}
}

/*	Given a block and the number of edges coming from 'Set', it splits
	this block and then adds to 'subpart2'. */
void SPLIT_BLOCK(BLOCK* block, int* n_fromSet, PART** subpart2)
{
	PART* splitted_blocks = NULL;
	PART* current_part;
	DoublyLinkNode* nodelist = block->head;
	
	while(nodelist)
	{
		Push_On_Block(nodelist->data, n_fromSet[nodelist->data], &splitted_blocks);
		nodelist = nodelist->next;
	}

/*	'splitted' has all the splitted blocks from 'block'. Now we find the largest
	block on it, and assigned its index as '-10' to the enqueue operation. */
	int index = 0;
	current_part = splitted_blocks;
	while(current_part)
	{
		current_part->block->index = index;
		current_part = current_part->next;
		index++;
	}

	int max_index = -1;
	int size = -1;
	current_part = splitted_blocks;
	while(current_part)
	{
		if((current_part->block->size)>size) { size = current_part->block->size; max_index = current_part->block->index; }
		current_part = current_part->next;	
	}

	current_part = splitted_blocks;
	while(current_part)
	{
		if(current_part->block->index==max_index) current_part->block->index = BIGGERFLAG;
		current_part = current_part->next;
	}	
/*	///////////////////////////////////////////////////////////////////////////	*/

	current_part = splitted_blocks;
	while(current_part)
	{
		push_block(subpart2, current_part->block);
		current_part = current_part->next;
	}
}

void BLOCKS_PARTITIONING(PART** subpart1, PART** subpart2, int* n_fromSet, Graph* graph, BLOCK* Set)
{
	PART* part_to_split = *subpart1;	
	DoublyLinkNode* nodelist = NULL;		
	
	while(part_to_split) // for each block to be splitted.
	{
		SPLIT_BLOCK(part_to_split->block, n_fromSet, subpart2);
		part_to_split = part_to_split->next;
	}
}

void S_SPLIT(PART** partition, BLOCK* Set, Graph* graph, QBLOCK** qhead, QBLOCK** qtail, int* onset)
{	
	int N, i;
	PART* subpart1 = NULL;
	PART* subpart2 = NULL;
	PART* current_part = NULL;

	N = 0;
	DoublyLinkNode* nodelist = Set->head;
	while (nodelist)
	{
		if(onset[nodelist->data]==0) N = 1;
		nodelist = nodelist->next;
	}
	if(N==0) return;

	//for(nodelist=Set->head; nodelist!=NULL; nodelist = nodelist->next)
	

/*	Precompute the number of edges coming from 'Set' for each node in the network. */
	int* n_fromSet = (int*)malloc((graph->size)*sizeof(int));
	for(i=0; i<(graph->size); i++) n_fromSet[i] = intersection_edges(graph, i, Set);

/*	Given the 'Set' block, now we select all the blocks in the partition that have at 
	least one connection coming from 'Set'. */
	GET_NONSTABLE_BLOCKS(partition, &subpart1, n_fromSet, graph, Set);

/*	Now 'subpart1' contains all the blocks that have nonzero pointed nodes from 'Set'. 
	Then, for each of these blocks, we split the ones that have different number of pointed
	links from 'Set' between their nodes. */
	//printf("Divide:\n");
	//printf("Number of blocks selected to be splitted -> ");
	//printPartitionSize(subpart1);
	//printAllPartition(subpart1);
	BLOCKS_PARTITIONING(&subpart1, &subpart2, n_fromSet, graph, Set);	// *error -> fixed* //
	//printf("Number of splitted blocks -> ");
	//printPartitionSize(subpart2);
	//printAllPartition(subpart2);
	

/*	After the process above, 'subpart2' contains the resulted splitted blocks. And now we
	delele all 'subpart1' blocks from partition and then we push all those 'subpart2' 
	blocks to the 'partition'. Before deleting 'subpart2', we enqueue all the blocks
	in the queue of refining blocks, except the largest one. */
	if(GetPartitionSize(subpart2)>GetPartitionSize(subpart1))
	{
		UPGRADE_PARTITION(&subpart2, &subpart1, partition);
	/*	Given a partition of blocks, insert all these blocks, except the largest ones, to
		the queue of refining blocks. */
		current_part = subpart2;
		while (current_part)
		{
			if(current_part->block->index!=(BIGGERFLAG)) enqueue_block(qhead, qtail, current_part->block);
			current_part = current_part->next;
		}
	}

	//printf("All Partition:\n");
	//printAllPartition(*partition);
	//printf("Queue After Divide:\n");
	//printQueue(*qhead);
	printPartitionSize(*partition);
	printQueueSize(*qhead);
	//free(subpart1);
	//free(subpart2);
	//free(n_fromSet);
}

void STABILITYCHECKER(PART** partition, BLOCK* Set, Graph* graph)
{	
	int N, i;
	PART* subpart1 = NULL;
	PART* subpart2 = NULL;
	PART* current_part = NULL;

	//for(nodelist=Set->head; nodelist!=NULL; nodelist = nodelist->next)
	

/*	Precompute the number of edges coming from 'Set' for each node in the network. */
	int* n_fromSet = (int*)malloc((graph->size)*sizeof(int));
	for(i=0; i<(graph->size); i++) n_fromSet[i] = intersection_edges(graph, i, Set);

/*	Given the 'Set' block, now we select all the blocks in the partition that have at 
	least one connection coming from 'Set'. */
	GET_NONSTABLE_BLOCKS(partition, &subpart1, n_fromSet, graph, Set);

/*	Now 'subpart1' contains all the blocks that have nonzero pointed nodes from 'Set'. 
	Then, for each of these blocks, we split the ones that have different number of pointed
	links from 'Set' between their nodes. */
	BLOCKS_PARTITIONING(&subpart1, &subpart2, n_fromSet, graph, Set);	// *error -> fixed* //
	

/*	After the process above, 'subpart2' contains the resulted splitted blocks. And now we
	delele all 'subpart1' blocks from partition and then we push all those 'subpart2' 
	blocks to the 'partition'. Before deleting 'subpart2', we enqueue all the blocks
	in the queue of refining blocks, except the largest one. */
	if(GetPartitionSize(subpart2)>GetPartitionSize(subpart1))
	{
		printf("Not stable\n");
	}
}

void ENQUEUE_SOLITAIRE(PART** null_partition, QBLOCK** qhead, QBLOCK** qtail)
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

int* UPGRADE_SET(BLOCK* block, int* onset)
{
	DoublyLinkNode* nodelist = block->head;
	while(nodelist)
	{
		onset[nodelist->data] = 1;
		nodelist = nodelist->next;
	}
	return onset;
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
		//add_to_block(P, i);
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
	int* onset = (int*)malloc(N*sizeof(int));
	for(i=0; i<N; i++) onset[i] = 0;  		
           
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
	ENQUEUE_SOLITAIRE(&null_partition, &qhead, &qtail);

	// Until L is empty, we procedure the splitting process.
	int time = 0;		
	BLOCK* CurrentSet;	
	while(qhead)
	{
		time++;		
		CurrentSet = dequeue_block(&qhead, &qtail);
		//printf("Time %d\n", time);
		//printf("Current Refinement Block: ");
		//printBlock(CurrentSet);
		//printf("Queue Before Divide:\n");
		//printQueue(qhead);
		S_SPLIT(&partition, CurrentSet, graph, &qhead, &qtail, onset);
		//if(time>1) onset = UPGRADE_SET(CurrentSet, onset);
	}
	int size = GetPartitionSize(partition);
	int presize = GetPartitionSize(null_partition);

	//PART* test = partition;
	//while(test)
	//{
	//	STABILITYCHECKER(&partition, test->block, graph);
	//	test = test->next;
	//}
	
	


	//printf("%d\n", presize);
	printf("Number of fiber: %d\n", size+presize);
	//printAllPartition(partition);
    ////////////////////////////////////////////////////////////////////////////////////

}
