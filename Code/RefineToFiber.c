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

int intersection_edges(Graph* graph, int node, BLOCK* Set)
{
	int i, check;	
	NodeAdj* Node = graph->array[node].head_in;
	BLOCK* SetHead = Set->head;

	int n_in = 0;
	while(Node)
	{
		check = doublycheck_element(SetHead, Node->neighbor);
		if(check==1) n_in++;
		Node = Node->next;
	}
	return n_in;	
}

void SPLIT(PART** partition, BLOCK* Set, Graph* graph, int N, QBLOCK** head, QBLOCK** tail)
{	
	PART* subpart1 = NULL;
	PART* subpart2 = NULL;
	PART* temp = partition;
	BLOCK* tempblock = NULL;	
	DoublyLinkNode* tempdoublylist;
	int* nedges = (int*)malloc(N*sizeof(int));

	/*	Given the 'Set' block, now we select all the blocks in the partition that have at 
		least one inward connection coming from 'Set'. */
	
	// 'subpart' will contains all the blocks that have nonzero	pointed nodes from 'Set'.
	while(temp)
	{
		tempdoublylist = temp->block->head;
		while(tempdoublylist)
		{
			int node = tempdoublylist->data;
			int n_in = intersection_edges(graph, node, Set);
			nedges[node] = n_in;
			if(n_in>0) { push_block(&subpart1, temp->block); break; }
			tempdoublylist = tempdoublylist->next;
		}
		temp = temp->next;
	}
	
	temp = subpart1;
	while(temp) // for each block to be splitted.
	{
		tempdoublylist = temp->block->head;
		while(tempdoublylist)	// for each node in block.
		{
			int node = tempdoublylist->data;			
			int n_in = intersection_edges(graph, node, Set);
			push_on(node, n_in, &subpart2);
		}
	}
	
	/*	After the process above, 'subpart2' contains the splitted blocks. And now we
		delele all 'subpart1' blocks from partition and then we push all those 'subpart2' 
		blocks to the 'partition'. Before deleting 'subpart2', we enqueue all the blocks
		in the queue of refining blocks, except the largest one. */
	 
	Upgrade_partition();
	EnqueueAllNotLargest();
}

BLOCK* DEFINEBLOCKS(BLOCK* P, int* nodecolor, int N)
{
    int i;
    for(i=0; i<N; i++) if(nodecolor[i]>-1) { push_doublylist(&P->head, i); P->size++; }
    return P;
}

void PREPROCESSING(Graph* graph, int* nodecolor, int N)
{
    int n_in, i;
    for(i=0; i<N; i++)
    {
        n_in = GETNin(graph, i);
        if(n_in==0) nodecolor[i] = -1;
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
	int ncolors = 1;
    int* nodecolor = (int*)malloc(N*sizeof(int));
    for(i=0; i<N; i++) nodecolor[i] = 0;        
	
	// Put aside from the algorithm (color -1) all nodes that don't receive any information.    
	PREPROCESSING(graph, nodecolor, N);

	// Define the initial partition as one block containing all operating nodes. 
	PART* partition = NULL;    
	BLOCK* P = (BLOCK*)malloc(sizeof(BLOCK));
    
	P->size = 0;
	P->index = 0;
    P->head = NULL;
    P = DEFINEBLOCKS(P, nodecolor, N);
    push_block(&partition, P);	// Push initial block to the partition.

	// Initialize the queue of blocks with the initial block above.
	QBLOCK* head = NULL;
	QBLOCK* tail = NULL;
	enqueue_block(&head, &tail, P);

	// Until L is empty, we procedure the splitting process.
	while(head)
	{
		BLOCK* Set = peek_block(&head);
		SPLIT(&partition, Set, graph, N, &head, &tail);
	}
    
    ////////////////////////////////////////////////////////////////////////////////////

}
