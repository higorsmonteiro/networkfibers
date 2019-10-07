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

int set_group(Graph* graph, int** table, int** edges, int ncolors, int* edgecolor, int* nodecolor, int  M)
{
    int i, j;
    int node1, node2, colorindex;
    for(i=0; i<M; i++)
    {
        colorindex = edgecolor[i];
        node2 = edges[i][1];
        table[node2][colorindex]++;
    }

    for(j=0; j<ncolors; j++)
    {
        
    }


}

int main() 
{ 
    int N;
    
	char genes[100] = "../Data/Ecoli/ngenes.dat";
    char genes_edges[100] = "../Data/Ecoli/edgelist.dat";
    FILE* UTIL = fopen(genes, "r");
    fscanf(UTIL, "%d\n", &N);

	// Creates the graph structure for N nodes.
    Graph* graph = createGraph(N);

    int** edges;
    int* regulator;
    edges = defineNetwork(edges, regulator, graph, genes_edges);

    ///// Minimal balanced coloring algorithm /////
	int M = nlines_file(genes_edges, 3);	  			// Effective number of voxels within the defined modules.
    
    int i, j;
    int ncolors = 1;

    int* nodecolor = (int*)malloc(N*sizeof(int));
    int* edgecolor = (int*)malloc(M*sizeof(int));
    for(i=0; i<N; i++) nodecolor[i] = 0;        // Initial state (all nodes have color '0').
    for(j=0; j<M; j++) edgecolor[j] = 0;        // Initial state (all links have color '0').

    int** table = arrint2d(N, ncolors);
    set_groups(graph, nodecolor, edgecolor, table);



    
}