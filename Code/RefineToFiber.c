/*	The code given here has the aim to determine the coarsest refinement partition sets of a graph to identify 
	the fibration building blocks of a directed gene regulatory network. The theory concerning the graph fibration 
	morphism and its application on biological networks are detailed, respectively, mainly in the paper of Boldi 
	and Vigna (2001) and in the paper of Morone et. al. (2019, to be published). Concerning the algorithm reproduced 
	in this code, one should refer to the work of Paige and Tarjan (1987).

	The resulted fibers represent, for each fiber, the set of nodes that are isomorphic under a graph fibration morphism, 
	i. e., all nodes belonging to the same fiber has identical input-set or input-trees, representing the classes of nodes 
	that receive equivalent information from other fibers on the network. Moreover, after the proper identification of the 
	building blocks, we classify each one based on specific topological features, represented by two parameter: |n,l>.

	---------------------------------------------------------------------------------------------------------------

	The code receives ONE command line argument, that is the common identifier for two required files: one single-line 
	file containing a integer number "%d\n" representing the number of nodes in the network ('ARG1Ngenes.dat'), and a edgelist 
	file ('ARG1edgelist.dat') containing all the directed links between nodes (3 columns: "%d\t%d\t%s\n" -> Pointing Node/ Pointed Node/ 
	Type of regulation). For gene regulatory networks, the type of the regulation can be positive, negative or dual.

	The result is stored in two arrays: 'nodecolor' and 'edgecolor'. Each array store the color of each component (node or 
	link) and, thus, gives all the information needed, together with the files, for the base graph construction.

	Author: Higor da S. Monteiro
	Email: higor.monteiro@fisica.ufc.br
	Complex System Lab - Departament of Physics/Universidade Federal do Ceará (UFC)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_eigen.h>
// Separated personal constructed modules for graph data and fibration specific functions.
#include "utilsforfiber.h"
#include "structforfiber.h"
#define SENTINEL -96

void CALCULATE_FUNDAMENTAL(PART** partition, Graph* graph)
{
	int i, j, k;
	int neigh, index, nn;
	int nregulators;
	int n_nodes;

	int* temp_index;
	double* temp_adjmatrix;

	int h = 0;
	PART* current_part;
	NODELIST* nodelist;
	for(current_part=(*partition); current_part!=NULL; current_part=current_part->next)
	{
		n_nodes = current_part->block->size;
		nregulators = current_part->number_regulators;
		nn = n_nodes+nregulators;
		printf("%d\n", nn);

		index = 0;
		temp_index = (int*)malloc(nn*sizeof(int));
    	temp_adjmatrix = (double*)malloc(nn*nn*sizeof(double));
		for(nodelist=current_part->regulators; nodelist!=NULL; nodelist=nodelist->next)
		{
			temp_index[index++] = nodelist->data;
			//index++;
		}
		for(nodelist=current_part->block->head; nodelist!=NULL; nodelist=nodelist->next)
		{
			temp_index[index++] = nodelist->data;
			//index++;
		}
		//Now we construct the current fiber adjacency matrix to calculate its eigenvalues.
		for(j=0; j<nn; j++)
    	{
			for(k=0; k<nn; k++)
			{
				neigh = CHECKLINK(graph, temp_index[j], temp_index[k]); // A_{jk} = 1 if j -> k
				if(neigh==1) temp_adjmatrix[j*nn + k] = 1.0;
				else temp_adjmatrix[j*nn + k] = 0.0;
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
		current_part->fundamental_number = eigmax;
		free(temp_index);
		free(temp_adjmatrix);
	}
}

void CALCULATE_REGULATORS(PART** partition, Graph* graph)
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
void GET_NONSTABLE_BLOCKS(PART** partition, PART** subpart, int* pos_fromSet, int* neg_fromSet, int* dual_fromSet, Graph* graph, BLOCK* Set)
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
				ndual_node = dual_fromSet[node];				// Number of 'Set' dual edges received by 'node'
				ndual_nextnode = dual_fromSet[nextnode];		// Number of 'Set' dual edges received by 'nextnode'
				if((npos_node!=npos_nextnode) || (nneg_node!=nneg_nextnode) || (ndual_node!=ndual_nextnode)) 
				{ push_block(subpart, current_part->block); break; }
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

void S_SPLIT(PART** partition, BLOCK* Set, Graph* graph, QBLOCK** qhead, QBLOCK** qtail)
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
		//neg_fromSet[i] += edgesfromSet(graph, i, Set, 2);
	} 
/*	Given the 'Set' block, now we select all the blocks in the partition that have at 
	least one connection coming from 'Set'. */
	GET_NONSTABLE_BLOCKS(partition, &subpart1, pos_fromSet, neg_fromSet, dual_fromSet, graph, Set);

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
}

int STABILITYCHECKER(PART** partition, BLOCK* Set, Graph* graph)
{	
	int i;
	PART* subpart1 = NULL;
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
	GET_NONSTABLE_BLOCKS(partition, &subpart1, pos_fromSet, neg_fromSet, dual_fromSet, graph, Set);
	
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
	// Define one block for each solitaire node and put in the 'Null_Part'.
	NODELIST* nodelist = (*NonP)->head;
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

    // Creates the network for N nodes and defines its structure with the given 'edgelist.dat' file.
	int** edges;
	Graph* graph = createGraph(N);
	edges = defineNetwork(edges, graph, net_edges);
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
	ENQUEUE_SOLITAIRE(&null_partition, &qhead, &qtail);

	// Until L is empty, we procedure the splitting process.
	BLOCK* CurrentSet;	
	while(qhead)
	{
		CurrentSet = dequeue_block(&qhead, &qtail);
		S_SPLIT(&partition, CurrentSet, graph, &qhead, &qtail);
	}
	int size = GetPartitionSize(partition) + GetPartitionSize(null_partition);
	int nfibers = GetFiberNumber1(partition, null_partition);
	// Partition contains all the fibers, including the trivial ones.

	/////////////////////////////// FIBER STATISTICS ////////////////////////////////////
	// Proper block unique indexation and external regulators initialization.
	PART* current_part;
	int index = 0;					
	int total_nodes = 0;			// Number of nodes in fibers.
	int non_trivial_fibers = 0;		// Number of non-trivial fibers (size larger than one).
	for(current_part=partition; current_part!=NULL; current_part=current_part->next)
	{
		current_part->regulators = NULL;	// Initialize the list of external regulators for each block.		
		current_part->number_regulators = 0;
		current_part->fundamental_number = 0.0;		
		current_part->block->index = index++;		
	}
	////////////////////////////////////////////////////////////////////////
	
	// Defines number of external regulators and set list of external regulators for each block.
	CALCULATE_REGULATORS(&partition, graph);
	CALCULATE_REGULATORS(&null_partition, graph);
	// Calculates fundamental class number for each fiber block.
	CALCULATE_FUNDAMENTAL(&partition, graph);
	////////////////////////////////////////////////////////////////////////////////////////////
	//printPartitionReg(partition);
	//printPartitionReg(null_partition);
	ShowMainInfo(partition);

	//printAllPartition(partition);
	

	//printPartitionSize(null_partition);
	//printf("%d %d\n", non_trivial_fibers, total_nodes);
	//printf("Number of fiber: %d\n", nfibers);
	//printAllPartition(partition);
    ////////////////////////////////////////////////////////////////////////////////////

}
