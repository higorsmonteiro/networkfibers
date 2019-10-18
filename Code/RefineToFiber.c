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
#include "utilsforfiber.h"   	// Personal module containing helpful functions to handle with network data and other common operations.
#define SEQSIZE 500     	// It's possible that the size of the string sequences should be larger for larger networks.

struct QueueOfBlocks
{
    int block_index;
    BLOCK* block;
    struct QueueOfBlocks* next;
};
typedef struct QueueOfBlocks QBLOCK;

BLOCK* DEFINEBLOCKS(BLOCK* P, int* nodecolor, int N)
{
    int i;
    for(i=0; i<N; i++) if(nodecolor[i]>-1) push_doublylist(&P->head, i);
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
    int N;                              // Number of nodes in the network.
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
    //////////////////////////////////

    // Creates the graph structure for N nodes (Not used in the balanced coloring algorithm)
    Graph* graph = createGraph(N);

    // Properly defines the network structure with the given 'edgelist.dat' file.
    int** edges;
    int* regulator;
    edges = defineNetwork(edges, regulator, graph, net_edges);

	/////////////////////// MINIMAL BALANCED COLORING ALGORITHM ////////////////////////
    int i, j, k;
    int M = nlines_file(net_edges, 3);              // Number of edges.  		
    
    // INITIAL STATE: All nodes and links have the same color '0' //
	int ncolors = 1;
    int* nodecolor = (int*)malloc(N*sizeof(int));
    for(i=0; i<N; i++) nodecolor[i] = 0;        
    PREPROCESSING(graph, nodecolor, N); // Nodes that do not receive any information are put aside (color -1) from the algorithm.
	
    BLOCK* P = (BLOCK*)malloc(sizeof(BLOCK));
    P->size = 0;
    P->head = NULL;

    P = DEFINEBLOCKS(P, nodecolor, N);
    printBlock(P);
    
    ////////////////////////////////////////////////////////////////////////////////////

}