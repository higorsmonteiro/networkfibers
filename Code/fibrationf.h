#ifndef FIBRATIONF_H
#define FIBRATIONF_H

/*  
    The code here represents the central module for the implementation of the coarsest
    refinement graph partitioning algorithm, containing the functions used in the main
    code. The necessary comments are given at the beginning of each function (to do).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
///////////////////
#include <gsl/gsl_eigen.h>
#include "utilsforfiber.h"
#include "structforfiber.h"
#define SENTINEL -96

extern double GET_EIGMAX(NODELIST* scc_nodes, Graph* graph)
{
	int i, j, k;
	int nn = GetListSize(scc_nodes);
	int* temp_index = (int*)malloc(nn*sizeof(int));
    double* temp_adjmatrix = (double*)malloc(nn*nn*sizeof(double));

	int index = 0;
	NODELIST* nodelist;
	for(nodelist=scc_nodes; nodelist!=NULL; nodelist=nodelist->next)
		temp_index[index++] = nodelist->data;

	for(j=0; j<nn; j++)
    {
		for(k=0; k<nn; k++)
		{
			int check = CHECKLINK(graph, temp_index[k], temp_index[j]); // A_{jk} = 1 if j -> k
			if(check==1) temp_adjmatrix[k*nn + j] = 1.0;
			else temp_adjmatrix[k*nn + j] = 0.0;
		}
	}

	/////////// GSL PACKAGE ROUTINES TO EIGENSYSTEMS PROBLEMS ////////////
    gsl_matrix_view m = gsl_matrix_view_array(temp_adjmatrix, nn, nn);
    // 'eval' will gonna stores all the 'nn' eigenvalues of the fiber adjacency matrix.
    gsl_vector_complex *eval = gsl_vector_complex_alloc (nn);

    gsl_eigen_nonsymm_workspace* w = gsl_eigen_nonsymm_alloc(nn);
    gsl_eigen_nonsymm(&m.matrix, eval, w);
    gsl_eigen_nonsymm_free(w);
    //////////////////////////////////////////////////////////////////////

    double temp;
    double eigmax = gsl_vector_get(eval, 0);
    for(j=1; j<nn; j++) { temp = gsl_vector_get(eval, j); if(temp>eigmax) eigmax = temp; }
	free(temp_index);
	free(temp_adjmatrix);
	return eigmax;
}

/*	Calculates the fundamental class number. For that, for each fiber block, we select
	all the nodes belonging to the block and all the external nodes that regulates that 
	same fiber.	After that, construct the adjacency matrix for that set and obtain its
	largest eigenvalue. */
extern void DEF_BRANCH_RATIO(PART** partition, Graph* graph)
{
	PART* current_part;
	NODELIST* nodelist;
	NODELIST* scc_nodes;
	// Loop over all the fiber blocks.
	for(current_part=(*partition); current_part!=NULL; current_part=current_part->next)
	{
		scc_nodes = NULL;
		// For each node in the fiber gets all nodes in its SCC.
		for(nodelist=current_part->block->head; nodelist!=NULL; nodelist=nodelist->next)
			KOSAJARU(&scc_nodes, nodelist->data, graph);
		double eigmax = GET_EIGMAX(scc_nodes, graph);
		current_part->fundamental_number = eigmax;
	}
}

extern void CALCULATE_REGULATORS(PART** partition, Graph* graph)
{
	int* in_neighbors;	
	int i, current_node, boolean_in, boolean_out;
	
	PART* current_part;
	NODELIST* nodelist;
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
					current_part->number_regulators++;
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
void Push_On_Block(int node, int nodeindex1, int nodeindex2, int nodeindex3, PART** subpart2)
{
/*	First thing, we need to loop over the blocks of the given partition
	and check if there is a block with the same index. */
	PART* current_part = *subpart2;
	NODELIST* nodelist = NULL;
	
	for(current_part=(*subpart2); current_part!=NULL; current_part=current_part->next)
	{
		if(current_part->block->pos==nodeindex1 && current_part->block->neg==nodeindex2 && current_part->block->dual==nodeindex3)
		{
			push_doublylist(&(current_part->block->head), node);
			current_part->block->size += 1;
			return;
		}
	}
/*	If there isn't a block with the same index, then we create a new one. */
	BLOCK* new_block = (BLOCK*)malloc(sizeof(BLOCK));
	new_block->index = -1;
	new_block->pos = nodeindex1;
	new_block->neg = nodeindex2;
	new_block->dual = nodeindex3;
	new_block->size = 1;
	new_block->head = NULL;
	push_doublylist(&(new_block->head), node);
	push_block(subpart2, new_block);
	return;
}

/*	Given the 'Set' block, now we select all the blocks in the 'partition' that aren't 
	stable with respect to 'Set'. We store these blocks into 'subpart' partition structure. */
extern void GET_NONSTABLE_BLOCKS1(PART** partition, PART** subpart, int* pos_fromSet, int* neg_fromSet, int* dual_fromSet, Graph* graph, BLOCK* Set)
{
	int node, nextnode;
	int npos_node, npos_nextnode;
	int nneg_node, nneg_nextnode;
	int ndual_node, ndual_nextnode;
	
	PART* current_part;
	NODELIST* nodelist;
	for(current_part=(*partition); current_part!=NULL; current_part=current_part->next)
	{
		nodelist = current_part->block->head;
		if(nodelist!=NULL)
		{
			node = nodelist->data;
			npos_node = pos_fromSet[node];
			nneg_node = neg_fromSet[node];
			ndual_node = dual_fromSet[node];
		} 		
		for(nodelist=current_part->block->head; nodelist!=NULL; nodelist=nodelist->next)
		{
			nextnode = nodelist->data;
			npos_nextnode = pos_fromSet[nextnode];		// Number of 'Set' positive edges received by 'nextnode'
			nneg_nextnode = neg_fromSet[nextnode];		// Number of 'Set' negatives edges received by 'nextnode'
			ndual_nextnode = dual_fromSet[nextnode];		// Number of 'Set' dual edges received by 'nextnode'
			if(((npos_node!=npos_nextnode) || (nneg_node!=nneg_nextnode)) || (ndual_node!=ndual_nextnode)) 
			{
				push_block(subpart, current_part->block); 
				break; 
			}
		}
	}
}

/*	Given a block and the number of edges coming from 'Set', it splits
	this block and then adds to 'subpart2'. */
void SPLIT_BLOCK(BLOCK* block, int* pos_fromSet, int* neg_fromSet, int* dual_fromSet, PART** subpart2)
{
	PART* current_part;
	PART* splitted = NULL;

	NODELIST* nodelist;
	for(nodelist=block->head; nodelist!=NULL; nodelist=nodelist->next)
		Push_On_Block(nodelist->data, pos_fromSet[nodelist->data], neg_fromSet[nodelist->data], dual_fromSet[nodelist->data], &splitted);

	int index = 0;
	current_part = splitted;
	for(current_part=splitted; current_part!=NULL; current_part=current_part->next)
	{
		current_part->block->index = index;
		index++;
	}
/*	'splitted' has all the splitted blocks from 'block'. Now we find the largest
	block on it, and assigned its index as 'SENTINEL' to a proper enqueue operation. */
	int max_index = -1;
	int size = -1;
	for(current_part=splitted; current_part!=NULL; current_part=current_part->next)
		if((current_part->block->size)>size) { size = current_part->block->size; max_index = current_part->block->index; }

	for(current_part=splitted; current_part!=NULL; current_part=current_part->next)
		if(current_part->block->index==max_index) current_part->block->index = SENTINEL;
/*	////////////////////////////////////////////////////////////////////////////////	*/

	for(current_part=splitted; current_part!=NULL; current_part=current_part->next)
		push_block(subpart2, current_part->block);
}

void BLOCKS_PARTITIONING(PART** subpart1, PART** subpart2, int* pos_fromSet, int* neg_fromSet, int* dual_fromSet, Graph* graph, BLOCK* Set)
{
	PART* part_to_split = *subpart1;	
	NODELIST* nodelist = NULL;		
	for(part_to_split=(*subpart1); part_to_split!=NULL; part_to_split = part_to_split->next)
		SPLIT_BLOCK(part_to_split->block, pos_fromSet, neg_fromSet, dual_fromSet, subpart2);
}

extern void S_SPLIT(PART** partition, BLOCK* Set, Graph* graph, QBLOCK** qhead, QBLOCK** qtail)
{	
	int N, i;
	PART* subpart1 = NULL;
	PART* subpart2 = NULL;
	PART* current_part = NULL;
	
/*	Precompute the number of edges coming from 'Set' for each node in the network. */
	int* pos_fromSet = (int*)malloc((graph->size)*sizeof(int));
	int* neg_fromSet = (int*)malloc((graph->size)*sizeof(int));
	int* dual_fromSet = (int*)malloc((graph->size)*sizeof(int));
	for(i=0; i<(graph->size); i++)
	{
		pos_fromSet[i] = edgesfromSet(graph, i, Set, 0);
		neg_fromSet[i] = edgesfromSet(graph, i, Set, 1);
		dual_fromSet[i] = edgesfromSet(graph, i, Set, 2);
	} 
/*	Given the 'Set' block, now we select all the blocks in the partition that have at 
	least one connection coming from 'Set'. */
	GET_NONSTABLE_BLOCKS1(partition, &subpart1, pos_fromSet, neg_fromSet, dual_fromSet, graph, Set);

/*	Now 'subpart1' contains all the blocks that have nonzero pointed nodes from 'Set'. 
	Then, for each of these blocks, we split the ones that have different number of pointed
	links from 'Set' between their nodes. */
	BLOCKS_PARTITIONING(&subpart1, &subpart2, pos_fromSet, neg_fromSet, dual_fromSet, graph, Set);

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
	free(pos_fromSet);
	free(neg_fromSet);
	free(dual_fromSet);
}

// Checks if 'partition' is input-tree stable with respect to 'Set'. //
int STABILITYCHECKER(PART** partition, BLOCK* Set, Graph* graph)
{	
/*	Precompute the number of edges coming from 'Set' for each node in the network. */
	int i;
	int* pos_fromSet = (int*)malloc((graph->size)*sizeof(int));
	int* neg_fromSet = (int*)malloc((graph->size)*sizeof(int));
	int* dual_fromSet = (int*)malloc((graph->size)*sizeof(int));	
	for(i=0; i<(graph->size); i++)
	{
		pos_fromSet[i] = edgesfromSet(graph, i, Set, 0);
		neg_fromSet[i] = edgesfromSet(graph, i, Set, 1);
		dual_fromSet[i] = edgesfromSet(graph, i, Set, 2);
	} 
	PART* subpart1 = NULL;
	GET_NONSTABLE_BLOCKS1(partition, &subpart1, pos_fromSet, neg_fromSet, dual_fromSet, graph, Set);
	free(pos_fromSet);
	free(neg_fromSet);
	free(dual_fromSet);

	if(GetPartitionSize(subpart1)>0) return -1;
	else return 1;
}

///////////////////////////////////////////////////////////////////////////////////////
////////////////// PREPROCESSING FUNCTIONS FOR REFINEMENT ALGORITHM ///////////////////
extern void ENQUEUE_BLOCKS(PART** partition, QBLOCK** qhead, QBLOCK** qtail)
{
	PART* current_part;
	for(current_part=(*partition); current_part!=NULL; current_part=current_part->next)
		enqueue_block(qhead, qtail, current_part->block);
}

extern void InitiateBlock(PART** Null_Part, int node)
{
	BLOCK* new_block = (BLOCK*)malloc(sizeof(BLOCK));
	new_block->index = -1;
	new_block->size = 0;
	new_block->head = NULL;
	add_to_block(&new_block, node);
	push_block(Null_Part, new_block);
}

//extern void PREPROCESSING(BLOCK** P, BLOCK** NonP, PART** Null_Part, Graph* graph, int N)
//{
//	(*P)->size = 0;
//	(*NonP)->size = 0;
//	(*P)->index = 0;
//	(*NonP)->index = -1;
//	(*P)->head = NULL;
//	(*NonP)->head = NULL;    
//
//	int n_in, i;
//    for(i=0; i<N; i++)
//    {
//		n_in = GETNin(graph, i);
//        if(n_in>0) add_to_block(P, i);
//		else add_to_block(NonP, i);
//    }
//	// Define one block for each solitaire node and put in the 'Null_Part'.
//	NODELIST* nodelist = (*NonP)->head;
//	for(nodelist=(*NonP)->head; nodelist!=NULL; nodelist=nodelist->next)
//        InitiateBlock(Null_Part, nodelist->data);
//}

extern void push_block_by_index(int node, int index, PART** partition)
{
/*	First thing, we need to loop over the blocks of the given partition
	and check if there is a block with the same index. */
	PART* current_part = *partition;
	NODELIST* nodelist = NULL;
	
	for(current_part=(*partition); current_part!=NULL; current_part=current_part->next)
	{
		if(current_part->block->index==index)
		{
			push_doublylist(&(current_part->block->head), node);
			current_part->block->size += 1;
			return;
		}
	}
/*	If there isn't a block with the same index, then we create a new one. */
	BLOCK* new_block = (BLOCK*)malloc(sizeof(BLOCK));
	new_block->index = index;
	new_block->size = 1;
	new_block->head = NULL;
	push_doublylist(&(new_block->head), node);
	push_block(partition, new_block);
	return;
}

extern void PREPROCESSING(PART** partition, PART** null_partition, PART** null_partition1, int* components, Graph* graph)
{    
	int N = graph->size;
	int root, n_in, n_out, i;
	int* roots = (int*)malloc(N*sizeof(int));
	//for(i=0; i<N; i++)
	//{
	//	n_in = GETNin(graph, i);
	//	if(n_in==0) roots[i] = -1;
	//	else
	//	{
	//		root = findroot(i, components);
	//		roots[i] = root;
	//	}
	//}
	for(i=0; i<N; i++)
	{
		n_out = GETNout(graph, i);
		n_in = GETNin(graph, i);
		int solitaire = IDENTIFY_SOLITAIRE(graph, i);
		if(solitaire==0) roots[i] = -1;
		else if(solitaire==1) roots[i] = -2;
		else
		{
			root = findroot(i, components);
			roots[i] = root;
		}
	}
	BLOCK* P = (BLOCK*)malloc(sizeof(BLOCK));
	BLOCK* NonP = (BLOCK*)malloc(sizeof(BLOCK));
	P->size = 0;
	P->index = -1;
	P->head = NULL;
	(NonP)->size = 0;
	(NonP)->index = -1;
	(NonP)->head = NULL;   
	for(i=0; i<N; i++)
	{
		if(roots[i]>=0) push_block_by_index(i, roots[i], partition);
		else if(roots[i]==-1) add_to_block(&NonP, i);
		else
		{
			add_to_block(&P, i);
			roots[i] = findroot(i, components);
			push_block_by_index(i, roots[i], partition);
		}
	}

	// Define one block for each solitaire node and put in the 'Null_Part'.
	NODELIST* nodelist = (NonP)->head;
	for(nodelist=(NonP)->head; nodelist!=NULL; nodelist=nodelist->next)
		InitiateBlock(null_partition, nodelist->data);
	for(nodelist=P->head; nodelist!=NULL; nodelist=nodelist->next)
		InitiateBlock(null_partition1, nodelist->data);        
}
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

#endif