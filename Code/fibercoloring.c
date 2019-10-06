/*	Code for phase synchronization evaluation on complex brain networks.
	To study the possible coherence states of the brain we apply the kuramoto
	model on functional brain networks (constructed by correlation between blood oxygen
	activity on fMRI measurements) to analyze the synchronization behavior between
	nodes and modules within the brain.
	Author: Higor da S. Monteiro - Universidade Federal do Ceará
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "graphdef.h"

int main() 
{ 
    int N;
    
	char genes[100] = "../Data/Ecoli/ngenes.dat";
    char genes_edges[100] = "../Data/Ecoli/edgelist.dat";
    FILE* UTIL = fopen(genes, "r");
    fscanf(UTIL, "%d\n", &N);

	// Creates the graph structure for N nodes.
    Graph* graph = createGraph(N);

	int nlinks = nlines_file(genes_edges, 2);	  			// Effective number of voxels within the defined modules.
	
    int** edges;
    int* regulator;
    defineNetwork(edges, regulator, graph, genes_edges);

    ///// Minimal balanced coloring algorithm /////
    int i, j;
    int ncolors = 1;

    int* colorset = (int*)malloc(N*sizeof(int));
    for(i=0; i<N; i++) colorset[i] = 0; // Initial state (all nodes have color '0').



    
}