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

void Upgrade_partition(PART** new_blocks, PART** old_blocks, PART** partition)
{
	BLOCK* tempblock = NULL;	
	PART* temp_old = *old_blocks;
	PART* temp_new = *new_blocks;
	while(temp_old) { deletePart(old_blocks, *old_blocks); }
	
}

void push_on_block(int node, int n_in, PART** phead)
{
	BLOCK* tempblock = NULL;
	PART* parttemp = *phead;
	while(parttemp)
	{
		tempblock = parttemp->block;
		if(tempblock->index==n_in) { add_to_block(&tempblock, node); return; }
	}
	// If the function reaches this line, it didnt find any block. So we create a new one.
	BLOCK* newblock = (BLOCK*)malloc(sizeof(BLOCK));
	newblock->index = n_in;
	newblock->size = 0;
	add_to_block(&newblock, node);
	push_block(phead, newblock);
}

int intersection_edges(Graph* graph, int node, BLOCK* Set)
{
	int i, check;	
	NodeAdj* Node = graph->array[node].head_in;

	int n_in = 0;
	while(Node)
	{
		check = doublycheck_element(Set->head, Node->neighbor);
		if(check==1) n_in++;
		Node = Node->next;
	}
	return n_in;	
}

void S_SPLIT(PART** partition, BLOCK* Set, Graph* graph, int N, QBLOCK** head, QBLOCK** tail)
{	
	PART* subpart1 = NULL;
	PART* subpart2 = NULL;
	BLOCK* tempblock = NULL;	
	PART* temp_part = *partition;
	DoublyLinkNode* tempdoublylist;

	/*	Given the 'Set' block, now we select all the blocks in the partition that have at 
		least one inward connection coming from 'Set'. */
	while(temp_part)
	{
		tempblock = temp_part->block;		// Get current block of the original partition.
		tempdoublylist = tempblock->head;	// Get the list of node of the current block.
		while(tempdoublylist)
		{
			int node = tempdoublylist->data;
			int n_in = intersection_edges(graph, node, Set);	// Get number of incoming links of from 'Set' to 'node'.
			if(n_in>0) { tempblock->index = n_in; push_block(&subpart1, tempblock); break; }
			tempdoublylist = tempdoublylist->next;
		}
		temp_part = temp_part->next;
	}
	// 'subpart' contains all the blocks that have nonzero pointed nodes from 'Set'.
	
	// New blocks creation.	
	temp_part = subpart1;
	while(temp_part) // for each block to be splitted.
	{
		tempblock = temp_part->block;		
		tempdoublylist = tempblock->head;
		while(tempdoublylist)	// for each node in the current block.
		{
			int node = tempdoublylist->data;			
			int n_in = intersection_edges(graph, node, Set);
			push_on_block(node, n_in, &subpart2);	// put 'node' in its proper new block.
		}
		temp_part = temp_part->next;
	}
	
	/*	After the process above, 'subpart2' contains the new splitted blocks. And now we
		delele all 'subpart1' blocks from partition and then we push all those 'subpart2' 
		blocks to the 'partition'. Before deleting 'subpart2', we enqueue all the blocks
		in the queue of refining blocks, except the largest one. */
	 
	Upgrade_partition(&subpart2, &subpart1, partition);
	//EnqueueAllNotLargest(&subpart2, qhead, qtail);
}

void PREPROCESSING(BLOCK** P, BLOCK** NonP, Graph* graph, int N)
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
    
    // INITIAL STATE: All nodes and links have the same color '0' //        
	// Put aside from the algorithm (color -1) all nodes that don't receive any information.    
	BLOCK* P = (BLOCK*)malloc(sizeof(BLOCK));
	BLOCK* NonP = (BLOCK*)malloc(sizeof(BLOCK));	
	PREPROCESSING(&P, &NonP, graph, N);

	// Define the initial partition as one block containing all operating nodes. 
	PART* partition = NULL;    
	push_block(&partition, P);	// Push initial block to the partition.

	// Initialize the queue of blocks with the initial block above.
	QBLOCK* qhead = NULL;
	QBLOCK* qtail = NULL;
	enqueue_block(&qhead, &qtail, P);

	// Until L is empty, we procedure the splitting process.
	//while(qhead)
	//{
		//BLOCK* Set = peek_block(&qhead);
		//S_SPLIT(&partition, Set, graph, N, &qhead, &qtail);
	//}
    ////////////////////////////////////////////////////////////////////////////////////

}
