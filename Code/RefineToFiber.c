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
#define SENTINEL -96

void CALCULATE_REGULATORS(PART** partition, Graph* graph)
{
	int* in_neighbors;	
	int i, current_node, boolean_in, boolean_out;
	
	PART* current_part;
	DoublyLinkNode* nodelist;
	for(current_part=(*partition); current_part!=NULL; current_part=current_part->next)
	{
		for(nodelist=current_part->block->head; nodelist!=NULL; nodelist=nodelist->next)
		{
			current_node = nodelist->data;
			int n_in = GETNin(graph, current_node);
			in_neighbors = GET_INNEIGH(graph, current_node);
			for(i=0; i<n_in; i++)
			{
				boolean_in = doublycheck_element(current_part->block->head, in_neighbors[i]);				
				boolean_out = doublycheck_element(current_part->regulators, in_neighbors[i]);
				if(boolean_out==0 && boolean_in==0)
				{ 
					push_doublylist(&(current_part->regulators), in_neighbors[i]); 
					current_part->number_regulators += 1;
				}
			}
		}
		free(in_neighbors);
	}
}

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
void Push_On_Block(int node, int nodeindex1, int nodeindex2, PART** subpart2)
{
/*	First thing, we need to loop over the blocks of the given partition
	and check if there is a block with the same index. */
	PART* current_part = *subpart2;
	DoublyLinkNode* nodelist = NULL;
	
	for(current_part=(*subpart2); current_part!=NULL; current_part=current_part->next)
	{
		if(current_part->block->aux1==nodeindex1 && current_part->block->aux2==nodeindex2)
		{
			push_doublylist(&(current_part->block->head), node);
			current_part->block->size += 1;
			return;
		}
	}
/*	If there isn't a block with the same index, then we create a new one. */
	BLOCK* new_block = (BLOCK*)malloc(sizeof(BLOCK));
	new_block->index = -1;
	new_block->aux1 = nodeindex1;
	new_block->aux2 = nodeindex2;
	new_block->size = 1;
	new_block->head = NULL;
	push_doublylist(&(new_block->head), node);
	push_block(subpart2, new_block);
	return;
}

/*	Given the 'Set' block, now we select all the blocks in the 'partition' that aren't 
	stable with respect to 'Set'. We store these blocks into 'subpart' partition structure. */
void GET_NONSTABLE_BLOCKS(PART** partition, PART** subpart, int* pos_fromSet, int* neg_fromSet, Graph* graph, BLOCK* Set)
{
	int node, nextnode;
	int npos_node, npos_nextnode, nneg_node, nneg_nextnode;
	
	PART* current_part;
	DoublyLinkNode* nodelist;
	for(current_part=(*partition); current_part!=NULL; current_part=current_part->next)
	{
		nodelist = current_part->block->head;
		for(nodelist=current_part->block->head; nodelist!=NULL; nodelist=nodelist->next)
		{
			if(nodelist->next!=NULL)
			{
				node = nodelist->data;
				nextnode = nodelist->next->data;
				npos_node = pos_fromSet[node];				// Number of 'Set' positive edges received by 'node'
				npos_nextnode = pos_fromSet[nextnode];		// Number of 'Set' positive edges received by 'nextnode'
				nneg_node = neg_fromSet[node];				// Number of 'Set' negatives edges received by 'node'
				nneg_nextnode = neg_fromSet[nextnode];		// Number of 'Set' negatives edges received by 'nextnode'
				if((npos_node!=npos_nextnode) || (nneg_node!=nneg_nextnode)) { push_block(subpart, current_part->block); break; }
			}
		}
	}
}

/*	Given a block and the number of edges coming from 'Set', it splits
	this block and then adds to 'subpart2'. */
void SPLIT_BLOCK(BLOCK* block, int* pos_fromSet, int* neg_fromSet, PART** subpart2)
{
	PART* current_part;
	PART* splitted_blocks = NULL;

	DoublyLinkNode* nodelist;
	for(nodelist=block->head; nodelist!=NULL; nodelist=nodelist->next)
		Push_On_Block(nodelist->data, pos_fromSet[nodelist->data], neg_fromSet[nodelist->data], &splitted_blocks);

	int index = 0;
	current_part = splitted_blocks;
	for(current_part=splitted_blocks; current_part!=NULL; current_part=current_part->next)
	{
		current_part->block->index = index;
		index++;
	}
/*	'splitted' has all the splitted blocks from 'block'. Now we find the largest
	block on it, and assigned its index as 'SENTINEL' to a proper enqueue operation. */
	int max_index = -1;
	int size = -1;
	for(current_part=splitted_blocks; current_part!=NULL; current_part=current_part->next)
		if((current_part->block->size)>size) { size = current_part->block->size; max_index = current_part->block->index; }

	for(current_part=splitted_blocks; current_part!=NULL; current_part=current_part->next)
		if(current_part->block->index==max_index) current_part->block->index = SENTINEL;
/*	////////////////////////////////////////////////////////////////////////////////	*/

	for(current_part=splitted_blocks; current_part!=NULL; current_part=current_part->next)
		push_block(subpart2, current_part->block);
}

void BLOCKS_PARTITIONING(PART** subpart1, PART** subpart2, int* pos_fromSet, int* neg_fromSet, Graph* graph, BLOCK* Set)
{
	PART* part_to_split = *subpart1;	
	DoublyLinkNode* nodelist = NULL;		
	
	for(part_to_split=(*subpart1); part_to_split!=NULL; part_to_split = part_to_split->next)
		SPLIT_BLOCK(part_to_split->block, pos_fromSet, neg_fromSet, subpart2);
}

void S_SPLIT(PART** partition, BLOCK* Set, Graph* graph, QBLOCK** qhead, QBLOCK** qtail, int* onset)
{	
	int N, i;
	PART* subpart1 = NULL;
	PART* subpart2 = NULL;
	PART* current_part = NULL;

	//N = 0;
	//DoublyLinkNode* nodelist;
	//for(nodelist=Set->head; nodelist!=NULL; nodelist = nodelist->next) if(onset[nodelist->data]==1) N=1;
	//if(N==1) return;
	
/*	Precompute the number of edges coming from 'Set' for each node in the network. */
	int* pos_fromSet = (int*)malloc((graph->size)*sizeof(int));
	int* neg_fromSet = (int*)malloc((graph->size)*sizeof(int));
	for(i=0; i<(graph->size); i++)
	{
		pos_fromSet[i] = edgesfromSet(graph, i, Set, 0);
		neg_fromSet[i] = edgesfromSet(graph, i, Set, 1);
		//pos_fromSet[i] += edgesfromSet(graph, i, Set, 2);
		//neg_fromSet[i] += edgesfromSet(graph, i, Set, 2);
	} 
/*	Given the 'Set' block, now we select all the blocks in the partition that have at 
	least one connection coming from 'Set'. */
	GET_NONSTABLE_BLOCKS(partition, &subpart1, pos_fromSet, neg_fromSet, graph, Set);

/*	Now 'subpart1' contains all the blocks that have nonzero pointed nodes from 'Set'. 
	Then, for each of these blocks, we split the ones that have different number of pointed
	links from 'Set' between their nodes. */
	BLOCKS_PARTITIONING(&subpart1, &subpart2, pos_fromSet, neg_fromSet, graph, Set);

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
		for(current_part=subpart2; current_part!=NULL; current_part=current_part->next)
			if(current_part->block->index!=(SENTINEL)) enqueue_block(qhead, qtail, current_part->block);
	}
}

int STABILITYCHECKER(PART** partition, BLOCK* Set, Graph* graph)
{	
	int i;
	PART* subpart1 = NULL;
/*	Precompute the number of edges coming from 'Set' for each node in the network. */
	int* pos_fromSet = (int*)malloc((graph->size)*sizeof(int));
	int* neg_fromSet = (int*)malloc((graph->size)*sizeof(int));
	for(i=0; i<(graph->size); i++)
	{
		pos_fromSet[i] = edgesfromSet(graph, i, Set, 0);
		neg_fromSet[i] = edgesfromSet(graph, i, Set, 1);
		//pos_fromSet[i] += edgesfromSet(graph, i, Set, 2);
		//neg_fromSet[i] += edgesfromSet(graph, i, Set, 2);
	} 
/*	Given the 'Set' block, now we select all the blocks in the partition that have at 
	least one connection coming from 'Set'. */
	GET_NONSTABLE_BLOCKS(partition, &subpart1, pos_fromSet, neg_fromSet, graph, Set);
	
	free(pos_fromSet);
	free(neg_fromSet);
	if(GetPartitionSize(subpart1)>0) return -1;
	else return 1;
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

int UPGRADE_SET(BLOCK* block, int* onset)
{
	DoublyLinkNode* nodelist;
	for(nodelist=block->head; nodelist!=NULL; nodelist=nodelist->next)	onset[nodelist->data] = 1;
	//return onset;
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
	// Define one block for each solitaire node and put in the 'Null_Part'.
	DoublyLinkNode* nodelist = (*NonP)->head;
	for(nodelist=(*NonP)->head; nodelist!=NULL; nodelist=nodelist->next)	InitiateBlock(Null_Part, nodelist->data);
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
    ///////////////////////////////////////////////////

    // Creates the graph structure for N nodes
	Graph* graph = createGraph(N);

    // Properly defines the network structure with the given 'edgelist.dat' file.
	int** edges;
	edges = defineNetwork(edges, graph, net_edges);
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
	int stability;
	int time = 0;		
	BLOCK* CurrentSet;	
	while(qhead)
	{
		time++;		
		CurrentSet = dequeue_block(&qhead, &qtail);
		S_SPLIT(&partition, CurrentSet, graph, &qhead, &qtail, onset);
		//stability = STABILITYCHECKER(&partition, CurrentSet, graph);
		//f(stability==1) continue;
		//printf("%d\n", stability);
		//printPartitionSize(partition);
		//if(time>1) UPGRADE_SET(CurrentSet, onset);
	}
	int size = GetPartitionSize(partition);
	int presize = GetPartitionSize(null_partition);
	// Partition contains all the fibers, including the trivial ones.

	/////////////////////// FIBER STATISTICS //////////////////////////
	PART* current_part;
	int index = 0;					// Proper block indexing.
	int total_nodes = 0;			// Number of nodes in fibers.
	int non_trivial_fibers = 0;		// Number of non-trivial fibers (size larger than one).
	for(current_part=partition; current_part!=NULL; current_part=current_part->next)
	{
		current_part->regulators = NULL;	// Initialize the list of external regulators for each block.		
		current_part->number_regulators = 0;		
		current_part->block->index = index++;		
		//if(current_part->block->size>1) { non_trivial_fibers++; total_nodes += current_part->block->size; }
	}
	
	CALCULATE_REGULATORS(&partition, graph);
	CALCULATE_REGULATORS(&null_partition, graph);
	//printPartitionReg(partition);
	//printPartitionReg(null_partition);

	printAllPartition(partition);
	

	//printPartitionSize(null_partition);
	//printf("%d %d\n", non_trivial_fibers, total_nodes);
	printf("Number of fiber: %d\n", size+presize);
	//printAllPartition(partition);
    ////////////////////////////////////////////////////////////////////////////////////

}
