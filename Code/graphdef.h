/*	This module contains a compilation of utility functions of some
	type of situations. Specially, here we define the functions necessary
	to construct directed and undirected networks. 

	Author: Higor da S. Monteiro - Universidade Federal do Ceará
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

////////////////////////////////////////////////////////////////////
/// Structures to create and define a directed, unweighted graph ///
////////////////////////////////////////////////////////////////////

// Define a graph containing 'size' nodes
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
	struct NodeAdj* next;
    //struct NodeAdj* next_out;
};
typedef struct NodeAdj NodeAdj;

NodeAdj* createNode(int neigh, int type_link)
{
	NodeAdj* newnode = (NodeAdj*)malloc(sizeof(NodeAdj));
	newnode->neighbor = neigh;
	newnode->next = NULL;
    newnode->regulator = type_link;
    //newnode->next_out = NULL;
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

void addEdges(int** edges, Graph* graph, int* regulator, int nE)
{
	int index, j, node1, node2, reg;
	
	for(j=0; j<nE; j++)
	{
		// 'node1' -> 'node2' directed link.
        node1 = edges[j][0];
		node2 = edges[j][1];
		reg = regulator[j];

		NodeAdj* newnode1 = createNode(node1, reg);
		NodeAdj* newnode2 = createNode(node2, reg);

		newnode2->next = graph->array[node1].head_out;
		newnode1->next = graph->array[node2].head_in;
		graph->array[node1].head_out = newnode2;
		graph->array[node2].head_in = newnode1;
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
        else if(strcmp("+-", k)==0) regulator[j] = 2;
        else regulator[j] = -1;
	}
	fclose(EDGE_FILE);
	
    // Defines the network structure.
	addEdges(edges, graph, regulator, nlink);
    
    return edges;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//////////////////////////// GRAPH OPERATIONS ////////////////////////////////

struct neigh_vector
{
	int neigh;
	int color;
};
typedef struct neigh_vector neigh_vector;

// To sort in ascendent order //
int cmpcolor(const void *a, const void *b)
{
    neigh_vector *a1 = (neigh_vector *)a;
    neigh_vector *a2 = (neigh_vector *)b;
    if ((*a1).color > (*a2).color)
        return 1;
    else if ((*a1).color < (*a2).color)
        return -1;
    else
        return 0;
}

extern neigh_vector* GET_inNEIGH(Graph* graph, int* nodecolor, int node, int N)
{
	int n_in = 0;	
	NodeAdj* NODE = graph->array[node].head_in;
	while(NODE)
	{
		n_in++;
		NODE = NODE->next;
	}
	neigh_vector* neighbors = (neigh_vector*)malloc((n_in)*sizeof(neigh_vector));

	int n_index = 0;
	NODE = graph->array[node].head_in;
	while(NODE)
	{
		neighbors[n_index].neigh = NODE->neighbor;
		neighbors[n_index].color = nodecolor[NODE->neighbor];
		n_index++;
		NODE = NODE->next;
	}
	qsort(neighbors, n_in, sizeof(neighbors[0]), cmpcolor);
	return neighbors;
}

extern neigh_vector* GET_outNEIGH(Graph* graph, int* nodecolor, int node, int N)
{
	int n_out = 0;	
	NodeAdj* NODE = graph->array[node].head_out;
	while(NODE)
	{
		n_out++;
		NODE = NODE->next;
	}

	neigh_vector* neighbors = (neigh_vector*)malloc(n_out*sizeof(neigh_vector));

	n_out = 0;
	NODE = graph->array[node].head_out;
	while(NODE)
	{
		neighbors[n_out].neigh = NODE->neighbor;
		neighbors[n_out].color = nodecolor[NODE->neighbor];
		n_out++;
		NODE = NODE->next;
	}
	qsort(neighbors, n_out, sizeof(neighbors[0]), cmpcolor);
	
	return neighbors;
}

extern int GETNin(Graph* graph, int node)
{
	int n_in = 0;	
	NodeAdj* NODE = graph->array[node].head_in;
	while(NODE)
	{
		n_in++;
		NODE = NODE->next;
	}
	return n_in;
}

extern char* GETSEQUENCE(neigh_vector* vec, int n, const int STRSIZE)
{
	int i;	
	int n_color = 1;
	for(i=0; i<(n-1); i++)
		if(vec[i].color!=vec[i+1].color) n_color++;

	int* seq = (int*)malloc(n_color*sizeof(int));
	for(i=0; i<n_color; i++) seq[i] = 0;

	int index = 0;	
	int* which_color = (int*)malloc(n_color*sizeof(int));
	
	which_color[index] = vec[0].color;
	for(i=0; i<(n-1); i++)
	{
		seq[index]++;
		if(vec[i].color!=vec[i+1].color) { which_color[++index] = vec[i+1].color; }
	}
	seq[index]++;

	char* strvar1 = (char*)malloc(STRSIZE*sizeof(char));
	strcpy(strvar1, "");

	for(i=0; i<n_color; i++)
	{
		char temp[10];
		char ctemp[10];
		sprintf(temp, "%d", seq[i]);
		sprintf(ctemp, "%d", which_color[i]);
		strcat(strvar1, temp);
		strcat(strvar1, ctemp);
	}
	//printf("%s\n", strvar1);
	return strvar1;
}

extern void GetInNeighbors(Graph* graph, int node, int N)
{
	NodeAdj* Inode = graph->array[node].head_in;
	printf("Node %d: ", node);
	while(Inode)
	{
		printf("%d ", Inode->neighbor);
		Inode = Inode->next;
	}
	printf("\n");
}

extern void GetOutNeighbors(Graph* graph, int node, int N)
{
	NodeAdj* Inode = graph->array[node].head_out;
	printf("Node %d: ", node);
	while(Inode)
	{
		printf("%d ", Inode->neighbor);
		Inode = Inode->next;
	}
	printf("\n");
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
/////////////////// STACK DATA STRUCTURE ////////////////////
struct strstack
{
	char color[10];
	struct strstack *next;
};
typedef struct strstack STRSTACK;

struct stack
{
	int extn;
	int node_ID;
	struct stack *next;
};
typedef struct stack Stack;
////////////////////////////////
////// STACK OPERATIONS ///////

// To insert an element in the stack.
extern STRSTACK *pushstr(STRSTACK* top, char* color)
{
	STRSTACK *ptr;
	ptr = (STRSTACK*)malloc(sizeof(STRSTACK));
	strcpy(ptr->color, color);
	if(top==NULL)
	{
		ptr->next = NULL;
		top = ptr;
	}
	else
	{
		ptr->next = top;
		top = ptr;
	}
	return top;
}

// To delete an element in the stack.
extern STRSTACK *popstr(STRSTACK* top)
{
	STRSTACK *ptr;
	ptr = top;
	if(top!=NULL)
	{
		top = top->next;
		free(ptr);
	}
	return top;
}

// To insert an element in the stack.
extern Stack *push(Stack* top, int node)
{
	Stack *ptr;
	ptr = (Stack*)malloc(sizeof(Stack));
	ptr->node_ID = node;
	if(top==NULL)
	{
		ptr->next = NULL;
		top = ptr;
	}
	else
	{
		ptr->next = top;
		top = ptr;
	}
	return top;
}

// To delete an element in the stack.
extern Stack *pop(Stack* top)
{
	Stack *ptr;
	ptr = top;
	if(top!=NULL)
	{
		top = top->next;
		free(ptr);
	}
	return top;
}
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

////////////// LINKED LIST DATA STRUCTURE ////////////////

struct Node
{
	int data;
	struct Node* next;
}
typedef struct Node Node;

extern Node* push_data(Node* head, int data)
{
	if(head==NULL)
	{
		head = (Node*)malloc(sizeof(Node));
		head->data = data;
		head->next = NULL;
	}
	else
	{
		Node* temp = head;		
		while(temp->next) { temp = temp->next; }
		temp->next = (Node*)malloc(sizeof(Node));
		temp = temp->next;
		temp->data = data;
		temp->next = NULL;
	}
	return head;
}

extern Node* push_to_head(Node* head, int data)
{
	Node* newnode = (Node*)malloc(sizeof(Node));
	newnode->data = data;
	newnode->next = head;
	head = newnode;
	return head;
}

extern Node* pop_head(Node* head)
{
	if(head==NULL) return NULL;	
	Node* newhead = head->next;
	free(head);
	head = newhead;
	return head;
}

extern Node* pop(Node* head)
{
	if(head==NULL) return NULL;	
	if(head->next==NULL) free(head);
	else
	{
		/* get to the second to last node in the list */
    	Node* current = head;
    	while (current->next->next != NULL) current = current->next;

    	free(current->next);
    	current->next = NULL;
	}
}

extern Node* remove_by_index(Node* head, int index)
{
	int i = 0;
    Node* current = head;
    Node* temp_node = NULL;

    if(index==0)	return pop_head(head);

    for (i = 0; i < n-1; i++) 
	{
        if (current->next == NULL)	return head;
        current = current->next;
    }

    temp_node = current->next;
    current->next = temp_node->next;
    free(temp_node);
	return head;
}

extern int check_value(Node* head, int value)
{
    Node* current = head;

	if(current->data==value) return 1;
	
	while(current)
	{
		if(current->next==NULL) return 0;		
		if(current->data==value) return 1;	
		else
			current = current->next;
	}
	return 0;
}

//////////////////////////////////////////////////////////

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
		printf("NODE %d\nin:", j);
        NodeAdj* SeeAux = graph->array[j].head_in;
        NodeAdj* SeeAux1 = graph->array[j].head_out;
		while(SeeAux)
		{
			printf("<-%d", SeeAux->neighbor);
			SeeAux = SeeAux->next; 
		}
		printf("\nout:");
        while(SeeAux1)
        {
            printf("->%d", SeeAux1->neighbor);
			SeeAux1 = SeeAux1->next;
        }
        printf("\n");
	}
}

extern int** arrint2d(int lines, int columns)
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

extern void free2d(int** arr, int rows, int columns)
{
	int i;
	for(i=0; i<rows; i++) free(arr[i]);
	free(arr);
}

extern double largest(double* arr, int N)
{
	int i;
	double max = arr[0];
	if(N>1) for(i=1; i<N; i++) if(arr[i]>max) max = arr[i];
	printf("max %lf\n", max);
	return max;
}

extern int INCHECK(int* arr, int size, int value)
{
	int i;
	for(i=0; i<size; i++) if(arr[i]==value) return 1;
	return 0;
}