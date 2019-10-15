/*  The code given here has the aim to determine the minimal balanced color classes to identify
    the fibration building blocks of a directed network. The theory concerning the graph fibration
    morphism and its application on biological networks are detailed, respectively, mainly in the
    paper of Boldi and Vigna (2001) and in the paper of Morone et. al. (2019, to be published). Concerning
    the algorithm reproduced in this code, one should refer to the work of Cardon and Crochemore (1982).

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
#include "graphdef.h"   	// Personal module containing helpful functions to handle with network data and other common operations.
#define SEQSIZE 500     	// It's possible that the size of the string sequences should be larger for larger networks.

int nfibers;

struct aux
{
	int color;
	int number;
}
typedef struct aux AUX;

struct strseq
{
    int node;
    Node* strseq;       
};
typedef struct strseq STRSEQ;

// Comparison function to be used in string sorting operation.
int cmp(const void *a, const void *b)
{
    STRSEQ *a1 = (STRSEQ *)a;
    STRSEQ *a2 = (STRSEQ *)b;

    return strcmp((*a1).strseq, (*a2).strseq);
}

STRSEQ* INITSEQ(STRSEQ* seqlist, int N)
{
	int i;
	for(i=0; i<N; i++) seqlist[i].strseq = NULL;
	return seqlist;
}
/////////////////////////////////////////////////////////////// 

void OPERATION(Graph* graph, Stack* color_index, int ncolor, int* nodecolor, int M, int N)
{
	int i,j;
	int n_in, current_color;
	neigh_vector* inlist;
	
	Node* List;
	Node* ListColor;
	
	STRSEQ* seqlist = (STRSEQ*)malloc(N*sizeof(STRSEQ));
	seqlist = INITSEQ(seqlist, N);
	
	// Determination of 'List', 'ListColor' and the sequences from 'seqlist'.	
	while(color_index)	// For each color class.
	{
		current_color = color_index->node_ID;
		color_index = pop(color_index);	
		
		for(i=0; i<N; i++)
		{
			if(nodecolor[i]==current_color)
			{
				inlist = GET_inNEIGH(graph, nodecolor, i, N);	// Get all the nodes that have outgoing links to the current node 'i'.
				n_in = GETNin(graph, i);						// Number of inbound links
				for(j=0; j<n_in; j++)
				{
					seqlist[inlist[j]]->strseq = push_data(seqlist[inlist[j]]->strseq, current_color);	// add 'current_color' to sequence of 'j'.					
					if(check_value(List, inlist[j].node)==0) List = push_data(List, inlist[j].node)
				}
			} 
		}
	}
	// At the end of the process above, the linked list 'List' stores all the nodes that have non-empty classes.
	// From that, the linked list 'ListColor' will stores the color of each node in 'List'.
	Node* temp = List;
	while(temp)
	{		
		int w = temp->data;
		ListColor = push_data(ListColor, nodecolor[w]);
	}
	// Sort in ascending order the sequences list for each node. *
	for(j=0; j<N; j++) SORTLIST(seqlist[j]->strseq, List);

	temp = ListColor;
	int* nclass = (int*)malloc(ncolor*sizeof(int));	// index 'i' stores the number of nodes in the class 'i'.
	for(i=0; i<ncolor; i++) nclass[i] = 0;			// Initial setting.
	while(temp)
	{
		nclass[temp->data]++;
		temp = temp->next;
	}

	while(List)
	{
		int vertex = List->data;
		List = pop_head(List);
		int i = nodecolor[vertex];
		if(nclass[i]!=0)
		{
			T[i].strseq = seqlist[vertex].strseq;
			R[i] = push(R[i], i) // wat?
			if(nclass[i]!=sizecolor[i]) // add a new index to R[i] //
			nclass[i] = 0;
		}	
	}
	
	
}

void SETFIBERS(Graph* graph, int** edges, int* edgecolor, int* nodecolor, int** table, int ncolors, int N, int M)
{
	Stack* color_index = (Stack*)malloc(ncolors*sizeof(Stack));
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
    int* edgecolor = (int*)malloc(M*sizeof(int));
    for(i=0; i<N; i++) nodecolor[i] = 0;        
    for(j=0; j<M; j++) edgecolor[j] = 0;

	SETFIBERS(graph, edges, edgecolor, nodecolor, table, ncolors, N, M);
	////////////////////////////////////////////////////////////////////////////////////

}