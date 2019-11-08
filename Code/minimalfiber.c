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
	that, it is necessary an auxiliary file called 'ARG1nameID.dat' containing two columns (formatted as "%s\t%d\n" -> Gene 
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
//#include <gsl/gsl_eigen.h>
// Separated personal constructed modules for graph data and graph fibration specific functions.
#include "fibrationf.h"
#include "utilsforfiber.h"
#include "structforfiber.h"
////////////////////////////////////////////////////////////////////////////////////////////////

void main(int argv, char** argc) 
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

    // Creates the network for N nodes and defines its structure with the given edgelist file.
	int** edges;
	int* components = (int*)malloc(N*sizeof(int));
	Graph* graph = createGraph(N, nodename, nodename_bool);
	edges = defineNetwork(edges, components, graph, net_edges);
	///////////////////////////////////////////////////////////////////////////////////////

	/////////////////////// COARSEST REFINEMENT PARTITIONING ALGORITHM ////////////////////////
	int i, j, k;
	
	// Define the initial partition as one block containing all operating nodes. 
	PART* partition = NULL;    
	PART* null_partition = NULL;
	PART* null_partition2 = NULL;
	PREPROCESSING(&partition, &null_partition, &null_partition2, components, graph);

	// Initialize the queue of blocks with the initial block above.
	QBLOCK* qhead = NULL;
	QBLOCK* qtail = NULL;
	ENQUEUE_BLOCKS(&partition, &qhead, &qtail);
	ENQUEUE_BLOCKS(&null_partition, &qhead, &qtail);
	ENQUEUE_BLOCKS(&null_partition2, &qhead, &qtail);
	//printAllPartition(null_partition2);

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
	// Calculates branch ratio number for each fiber block.
	DEF_FUNDAMENTAL(&partition, graph);
	//DEF_BRANCH_RATIO(&partition, graph);
	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	
	////// 'nodefibers' directly relates nodes with their fiber index ///////
	NODELIST* nodelist;	
	int total_nodes = 0;	// Number of nodes inside non-trivial fibers.
	int* nodefibers = (int*)malloc(N*sizeof(int));
	for(i=0; i<N; i++) nodefibers[i] = -1;

	for(current_part=partition; current_part!=NULL; current_part=current_part->next)
	{
		if(current_part->block->size>1) { total_nodes+=current_part->block->size; }
		for(nodelist=current_part->block->head; nodelist!=NULL; nodelist=nodelist->next)
			nodefibers[nodelist->data] = current_part->block->index;
	}
	////////////////////////////////////////////////////////////////////////////////////

	// Uncomment line below to get the fiber input details.
	/* With gene names */ //printGeneGraphInFibers(graph, partition, nodefibers);
	/* Without gene names */ //printGraphInFibers(graph, partition, nodefibers);
	
	/* Classification Info */ //ShowClassification1(partition, 0);
	/* Fiber blocks and classification info */ ShowInfo(partition, 0);

	/*	Show the number of non-trivial fibers	*/
	//printf("%d %d\n", nontrivial_fibers, total_nodes);
	
	//int node = atoi(argc[3]);
	//PrintInNeighbors(graph, node);
	//PrintOutNeighbors(graph, node);
	//int n = GETNin(graph, node);

	//printf("%d\n", n);
	/////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////

}
