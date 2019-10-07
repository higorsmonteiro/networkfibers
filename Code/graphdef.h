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

////////////////////////////////////////////////////////////////////
/// Structures to create and define a directed, unweighted graph ///
////////////////////////////////////////////////////////////////////

// Define a graph containing N nodes
struct Graph
{
	// each node has its adjancency list
	int size;
	struct adjList* array;
};
typedef struct Graph Graph;

struct adjList
{
	struct NodeAdj* head_in;
    struct NodeAdj* head_out;
};
typedef struct adjList AdjList;

// Structure to define the adjacency list of a node
struct NodeAdj
{
	int neighbor;
    int regulator;
	struct NodeAdj* next_in;
    struct NodeAdj* next_out;
};
typedef struct NodeAdj NodeAdj;

NodeAdj* createNode(int neigh)
{
	NodeAdj* newnode = (NodeAdj*)malloc(sizeof(NodeAdj));
	newnode->neighbor = neigh;
	newnode->next_in = NULL;
    newnode->next_out = NULL;
    newnode->regulator = -2;
	return newnode;
}

extern Graph* createGraph(int N)
{
	Graph* graph = (Graph*)malloc(sizeof(Graph));
	graph->size = N;

	// Array of NodeAdj's
	graph->array = (AdjList*)malloc(N*sizeof(AdjList));

	int j;
	for(j=0; j<N; j++) { graph->array[j].head_in = NULL; graph->array[j].head_out = NULL; }
	return graph;
}

void addEdges(int** edges, Graph* graph, int nE)
{
	int index, j, node1, node2, reg;
	
	for(j=0; j<nE; j++)
	{
		// 'node1' -> 'node2' directed link.
        node1 = edges[j][0];
		node2 = edges[j][1];

		NodeAdj* newnode1 = createNode(node2);
		NodeAdj* newnode2 = createNode(node1);

		// The created nodes points to the head of graph and
		// next become the head itselfs.
        newnode1->next_out = graph->array[node1].head_out;
        graph->array[node1].head_out = newnode1;

        newnode2->next_in = graph->array[node2].head_in;
        graph->array[node2].head_in = newnode2;
	}
}

extern int** defineNetwork(int** edges, int* regulator, Graph* graph, char* filename)
{
	FILE *EDGE_FILE = fopen(filename, "r");
	if(EDGE_FILE==NULL) printf("ERROR in file reading");

	int i, j;
	char k[5];
	int r = 1;
	int nlink = 0; // number of links.
	while(r) // Calculates the number of lines in the file
	{
		r = fscanf(EDGE_FILE, "%d\t%d\t%s\n", &i, &j, &k);
		if(r==EOF) break;
		nlink++;
	}
	rewind(EDGE_FILE);

	edges = (int**)malloc(nlink*sizeof(int*));
    regulator = (int*)malloc(nlink*sizeof(int));
	for(j=0; j<nlink; j++)
	{
		edges[j] = (int*)malloc(2*sizeof(int));
		r = fscanf(EDGE_FILE, "%d\t%d\t%s\n", &edges[j][0], &edges[j][1], &k);
        if(strcmp("+", k)==0) regulator[j] = 0;
        else if(strcmp("-", k)==0) regulator[j] = 1;
        else if(strcmp("b", k)==0) regulator[j] = 2;
        else regulator[j] = -1;
        printf("%d\n", regulator[j]);
	}
	fclose(EDGE_FILE);
	
    // Defines the network structure.
	addEdges(edges, graph, nlink);
    
    return edges;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// Reads the number of effective voxels and total number of voxels of the subject map.
extern int nlines_file(char* filename, int ncolumns)
{
	FILE *MODFILE = fopen(filename, "r");
	if(MODFILE==NULL) printf("ERROR in file reading");

	char ch[5];
    int i, j, k, m, x, y, z;
	int ne = 0;
	int r = 1;
	while(r) // Calculates the number of lines in the file
	{
		if(ncolumns==3) r = fscanf(MODFILE, "%d\t%d\t%s\n", &i, &j, &ch);
		else if(ncolumns==6) r = fscanf(MODFILE, "%d\t%d\t%d\t%d\t%d\t%d\n", &i, &j, &k, &x, &y, &z);
		if(r==EOF) break;
		ne++;
	}
	fclose(MODFILE);
	return ne;
}

extern void printGraph(Graph* graph)
{
	int j;
	for(j=0; j<graph->size; j++)
	{
		printf("Node %d, in:\n", j);
        NodeAdj* SeeAux = graph->array[j].head_in;
        NodeAdj* SeeAux1 = graph->array[j].head_out;
		while(SeeAux)
		{
			printf("<- %d", SeeAux->neighbor);
			SeeAux = SeeAux->next_in; 
		}
		printf("\nout:\n");
        while(SeeAux1)
        {
            printf("-> %d", SeeAux1->neighbor);
			SeeAux1 = SeeAux1->next_out;
        }
        printf("\n");
	}
}

int** arrint2d(int lines, int columns)
{
    int i, j;
    int** arr = (int**)malloc(lines*sizeof(int*));
    for(i=0; i<lines; i++)
    {
        arr[i] = (int*)malloc(columns*sizeof(int));
        for(j=0; j<columns; j++)
            arr[i][j] = 0;
    }
    return arr;
}