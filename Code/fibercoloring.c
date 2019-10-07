/*	Code to determine the minimal balanced coloring for a proper identification
    of the fiber building blocks of a given network. The theory and description
    of the algorithm is detailed in the papers of F. Morone et. al. (2019) and Cardon 
    and Crochemore (1982).

    The resulted fibers classes distribution represents the set of nodes that are 
    isomorphics under a graph fibration transformation.

	Author: Higor da S. Monteiro - Universidade Federal do Ceará
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "graphdef.h"

struct sequences
{
    int node;
    char strseq[200];
};
typedef struct sequences STRSEQ;


// After the nodes and edges in the network has their colors set, we upgrade the
// values given in the table structure.
void upgrade_table(int M, int** table, int** edges, int* edgecolor, int* nodecolor)
{
    int j;
    int node_pointed, ind_color;
    for(j=0; j<M; j++)
    {
        ind_color = edgecolor[j];

        node_pointed = edges[j][1];
        table[node_pointed][ind_color] += 1;
    }
}

int check_states(int ncolors, int** table, int* nodecolor, int* edgecolor, int N, int M)
{
    int pluscolor;
    int i, j, nn, nodetemp;

    char* strvar;
    STRSEQ* seqlist;
    Stack* node_in_class = NULL;

    for(i=0; i<ncolors; i++)    // For each color class C_i.
    {
        nn = 0;          // Number of nodes inside the color class 'i'.
        pluscolor = 0;   // If the current color class is going to be partitioned, the 'pluscolor' stores the additional number of colors.
       
        // Stores in a stack all the nodes belonging to C_i.
        for(j=0; j<N; j++)   if(nodecolor[j]==i) { node_in_class = push(node_in_class, j); nn++; }

        seqlist = (STRSEQ*)malloc(nn*sizeof(STRSEQ));
        strvar = (char*)malloc((ncolors+1)*sizeof(char));
        strcpy(strvar, "");

        // Get node at the top of the 'node_in_class' stack until stack is EMPTY.
        while(node_in_class)     
        {
            nodetemp = node_in_class->node_ID;
           
            //  For each color class C_j.
            for(j=0; j<ncolors; j++)
            {
                char temp[10];
                sprintf(temp, "%d", table[nodetemp][j]);
                strcat(strvar, temp);
            }
            seqlist[nn-1].node = nodetemp;
            strcpy(seqlist[nn-1].strseq, strvar);

            node_in_class = pop(node_in_class);
            nn--;
        }

        // Now, 'seqlist' is a list of structures cointaining the node in the color class and its
        // string sequence of the number of in-links for all colors.

        // Sort lexicographically that list. The blocks that have the same sequence will be neighbors
        // so that we can find the number of partitions inside the class color C_i.

        


        free(seqlist);
        free(strvar);
    }
}

int** fibration(Graph* graph, int** edges, int* edgecolor, int* nodecolor, int** table, int ncolors, int N, int M)
{
    int i,j;
    table = arrint2d(N, ncolors);

    // Upgrade the initial state //
    upgrade_table(M, table, edges, edgecolor, nodecolor);

    int temp_ncolor = ncolors;
    do
    {
        temp_ncolor = check_states(ncolors, table, nodecolor, edgecolor, N, M);
    } while(temp_ncolor!=ncolors);


}

int main() 
{ 
    int N;      // Number of nodes in the network.
    
	char genes[100] = "../Data/Ecoli/ngenes.dat";           // File containing (one line) the number of nodes in the network.
    char genes_edges[100] = "../Data/Ecoli/edgelist.dat";   // File containing all the links in the network.
    
    // Defines the size of the network.
    FILE* UTIL = fopen(genes, "r");
    fscanf(UTIL, "%d\n", &N);

	// Creates the graph structure for N nodes.
    Graph* graph = createGraph(N);

    // Properly defines the network structure with the given 'edgelist.dat' file.
    int** edges;
    int* regulator;
    edges = defineNetwork(edges, regulator, graph, genes_edges);
    /////////////////////////////////////////////////////////////////////////////

    ///////////// Minimal balanced coloring algorithm //////////////
	int M = nlines_file(genes_edges, 3);	     // Number of edges.  		
    
    int i, j;

    // INITIAL STATE: ALL NODES AND LINKS HAVE THE SAME COLOR '0' //
    int ncolors = 1;
    int* nodecolor = (int*)malloc(N*sizeof(int));
    int* edgecolor = (int*)malloc(M*sizeof(int));
    for(i=0; i<N; i++) nodecolor[i] = 0;        
    for(j=0; j<M; j++) edgecolor[j] = 0;
    ////////////////////////////////////////////////////////////////        

    int** table;
    table = fibration(graph, edges, nodecolor, edgecolor, table, ncolors, N, M);



    
}