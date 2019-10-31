/*	The code given here has the aim to determine the coarsest refinement partition sets of a graph to identify 
	the fibration building blocks of a directed gene regulatory network. The theory concerning the graph fibration 
	morphism and its application on biological networks are detailed, respectively, mainly in the paper of Boldi 
	and Vigna (2001) and in the paper of Morone et. al. (2019, to be published). Concerning the algorithm reproduced 
	in this code, one should refer to the work of Paige and Tarjan (1987).

	The resulted fibers represent, for each fiber, the set of nodes that are isomorphic under a graph fibration morphism, 
	i. e., all nodes belonging to the same fiber has identical input-set or input-trees, representing the classes of nodes 
	that receive equivalent information from other fibers on the network. Moreover, after the proper identification of the 
	building blocks, we classify each one based on specific topological features, represented by two parameter: |n,l>.

	----------------------------------------------------------------------------------------------------------------------------

	The code receives TWO commands line arguments. The first one is the string identifier for the edgelist file 
	('ARG1edgelist.dat') containing all the directed links between nodes (3 columns: "%d\t%d\t%s\n" -> Pointing Node/ 
	Pointed Node/ Type of regulation). For gene regulatory networks, the type of the regulation can be 'positive', 'negative' 
	or 'dual'. The second argument is a flag used to signal the code to properly get the gene names of each node number. For 
	that, it is necessary an auxiliary file called 'ARG1genename.dat' containing two columns (formatted as "%s\t%d\n" -> Gene 
	name/ Gene ID number). Thus, if there is a gene name file, the code will properly link all the node numbers with their 
	corresponding name if 'ARG2' is passed as '-y', otherwise just the node numbers is stored for each node.

	The result is stored in the 'partition' and 'null_partition' structures, together with the 'graph' structure. To check 
	which data each one of this structures stores the user can refer to the 'structforfiber.h' module. In general, a partition 
	stores all the fiber blocks and each block stores the list of node that belongs to it. The values of n and l are stored in 
	'partition' and 'null_partition' for each block.

	------------------------------------------------------------------------------------------------------------------------------

	Author: Higor da S. Monteiro
	Email: higor.monteiro@fisica.ufc.br
	Complex System Lab (Prof. José Soares de Andrade Jr.)
	Departament of Physics/Universidade Federal do Ceará (UFC) - Fortaleza, Ceará.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// Used fot the calculation of the eigenvalues for each fibration block adjacency matrix.
#include <gsl/gsl_eigen.h>
// Separated personal constructed modules for graph data and graph fibration specific functions.
#include "utilsforfiber.h"
#include "structforfiber.h"
////////////////////////////////////////////////////////////////////////////////////////////////
#define SENTINEL -96


/*	Calculates the fundamental class number. For that, for each fiber block, we select
	all the nodes belonging to the block and all the external nodes that regulates that 
	same fiber.	After that, construct the adjacency matrix for that set and obtain its
	largest eigenvalue. */
void CALCULATE_FUNDAMENTAL(PART** partition, Graph* graph)
{
	int i, j, k;
	int check, index, nn;
	int all_nodes;	
	int number_fnodes;
	int number_regulators;

	int* temp_index;			// Temporary indexation for each node in the circuit set.
	double* temp_adjmatrix;		// n*n array to define the adjacency matrix.

	PART* current_part;
	NODELIST* nodelist;
	NODELIST* scc_nodes;
	// Loop over all the fiber blocks.
	for(current_part=(*partition); current_part!=NULL; current_part=current_part->next)
	{
		//temp_index = (int*)malloc(all_nodes*sizeof(int));
    	//temp_adjmatrix = (double*)malloc(all_nodes*all_nodes*sizeof(double));

		scc_nodes = NULL;
		// Get the correct nodes.
		for(nodelist=current_part->block->head; nodelist!=NULL; nodelist=nodelist->next)
		{
			push_doublylist(&scc_nodes, nodelist->data);
			KosajaruSCC(&scc_nodes, nodelist->data, graph);
		}
			
		//for(nodelist=current_part->regulators; nodelist!=NULL; nodelist=nodelist->next)
		//	temp_index[index++] = nodelist->data;
		
		//Now we construct the current circuit fiber adjacency matrix to calculate its eigenvalues.
		for(j=0; j<all_nodes; j++)
    	{
			for(k=0; k<all_nodes; k++)
			{
				check = CHECKLINK(graph, temp_index[j], temp_index[k]); // A_{jk} = 1 if j -> k
				if(check==1) temp_adjmatrix[j*all_nodes + k] = 1.0;
				else temp_adjmatrix[j*all_nodes + k] = 0.0;
			}
		}
		
		/////////// GSL PACKAGE ROUTINES TO EIGENSYSTEMS PROBLEMS ////////////
        gsl_matrix_view m = gsl_matrix_view_array(temp_adjmatrix, all_nodes, all_nodes);
        // 'eval' will gonna stores all the 'nn' eigenvalues of the fiber adjacency matrix.
        gsl_vector_complex *eval = gsl_vector_complex_alloc (all_nodes);

        gsl_eigen_nonsymm_workspace* w = gsl_eigen_nonsymm_alloc(all_nodes);
        gsl_eigen_nonsymm(&m.matrix, eval, w);
        gsl_eigen_nonsymm_free(w);
        //////////////////////////////////////////////////////////////////////

        double temp;
        double eigmax = gsl_vector_get(eval, 0);
        for(j=1; j<all_nodes; j++) { temp = gsl_vector_get(eval, j); if(temp>eigmax) eigmax = temp; }
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

void GET_NONSTABLE_BLOCKS1(PART** partition, PART** subpart, int* pos_fromSet, int* neg_fromSet, int* dual_fromSet, Graph* graph, BLOCK* Set)
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
			//node = nodelist->data;
			nextnode = nodelist->data;
			//npos_node = pos_fromSet[node];				// Number of 'Set' positive edges received by 'node'
			npos_nextnode = pos_fromSet[nextnode];		// Number of 'Set' positive edges received by 'nextnode'
			//nneg_node = neg_fromSet[node];				// Number of 'Set' negatives edges received by 'node'
			nneg_nextnode = neg_fromSet[nextnode];		// Number of 'Set' negatives edges received by 'nextnode'
			//ndual_node = dual_fromSet[node];				// Number of 'Set' dual edges received by 'node'
			ndual_nextnode = dual_fromSet[nextnode];		// Number of 'Set' dual edges received by 'nextnode'
			if(((npos_node!=npos_nextnode) || (nneg_node!=nneg_nextnode)) || (ndual_node!=ndual_nextnode)) 
			{
				push_block(subpart, current_part->block); 
				break; 
			}
		}
	}
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
				{
					push_block(subpart, current_part->block); 
					break; 
				}
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
	//free(pos_fromSet);
	//free(neg_fromSet);
	//free(dual_fromSet);
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
	GET_NONSTABLE_BLOCKS1(partition, &subpart1, pos_fromSet, neg_fromSet, dual_fromSet, graph, Set);
	
	free(pos_fromSet);
	free(neg_fromSet);
	free(dual_fromSet);
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
	char net_edges[100] = "../Data/";   // File containing all the directed links in the network.
	char nodename[100] = "../Data/";	// File containing all the nodes name.
	strcat(net_edges, argc[1]);
	strcat(net_edges, "edgelist.dat");
	strcat(nodename, argc[1]);
	strcat(nodename, "nameID.dat");

	// From the edgelist get the number of nodes in the network.
	N = GetNodeNumber(net_edges);

	//// Check if it is necessary to check for the file containing names for the nodes ////
	int nodename_bool;
	if(strcmp(argc[2], "-y")==0) nodename_bool = 1;
	else nodename_bool = 0;
	///////////////////////////////////////////////////////////////////////////////////////

    // Creates the network for N nodes and defines its structure with the given 'edgelist.dat' file.
	int** edges;
	Graph* graph = createGraph(N, nodename, nodename_bool);
	edges = defineNetwork(edges, graph, net_edges);
	//printf("%d\n", graph->size);
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
	int nontrivial_fibers = GetFiberNumber1(partition, null_partition);
	// 'partition' contains all the fibers, except the solitaire ones.

	/////////////////////////////// FIBER STATISTICS ////////////////////////////////////
	// Proper block unique indexation and external regulators initialization.
	PART* current_part;
	int index = 0;					
	for(current_part=partition; current_part!=NULL; current_part=current_part->next)
	{
		current_part->regulators = NULL;	// Initialize the list of external regulators for each block.		
		current_part->number_regulators = 0;
		current_part->fundamental_number = 0.0;		
		current_part->block->index = index++;		
	}

	// Defines number of external regulators and set list of external regulators for each block.
	CALCULATE_REGULATORS(&partition, graph);
	CALCULATE_REGULATORS(&null_partition, graph);
	// Calculates fundamental class number for each fiber block.
	CALCULATE_FUNDAMENTAL(&partition, graph);
	////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////
	
	////// 'nodefibers' directly relates nodes with their fiber index ///////
	int* nodefibers = (int*)malloc(N*sizeof(int));
	for(i=0; i<N; i++) nodefibers[i] = -1;

	NODELIST* nodelist;	
	int total_nodes = 0;	// Number of nodes inside non-trivial fibers.
	for(current_part=partition; current_part!=NULL; current_part=current_part->next)
	{
		if(current_part->block->size>1) { total_nodes+=current_part->block->size; }
		for(nodelist=current_part->block->head; nodelist!=NULL; nodelist=nodelist->next)
			nodefibers[nodelist->data] = current_part->block->index;
	}
	//////////////////////////////////////////////////////////////////////////

	ShowMainInfo(partition);
	
	int n;
	STORETYPE* tempx;
	int* temp;
	printf("%d\n", nontrivial_fibers);
	//for(current_part=partition; current_part!=NULL; current_part=current_part->next)
	//{
	//	if(current_part->block->size==1) continue;
	//	printf("INSIDE BLOCK %d:\n", current_part->block->index);
	//	for(nodelist=current_part->block->head; nodelist!=NULL; nodelist=nodelist->next)
	//	{
	//		printf("NODE %d receives from fiber(node,type): ", nodelist->data);
	//		n = GETNin(graph, nodelist->data);
	//		tempx = GET_INTYPENEIGH(graph, nodelist->data);
	//		for(i=0; i<n; i++) printf("%d(%d,%d) ", nodefibers[tempx[i].node], tempx[i].node, tempx[i].type);
	//		printf("\n");
	//	}
	//	printf("\n");
	//}
	//int ll = atoi(argc[3]);
	//PrintInNeighbors(graph, ll, N);
	//PrintOutNeighbors(graph, ll, N);
    ////////////////////////////////////////////////////////////////////////////////////

}
