/*  The code given here has the aim to determine the minimal balanced color classes for a 
    proper identification of the fibration building blocks of a directed network. The theory
    concerning the graph fibration morphism and its application on biological networks are 
    detailed, respectively, mainly in the paper of Boldi and Vigna (2001) and in the paper of Morone
    et. al. (2019, to be published). Concerning the algorithm reproduced in this code, one should refer
    to the work of Cardon and Crochemore (1982).

    The resulted fibers represent, for each fiber, the set of nodes that are isomorphics under a graph
    fibration morphism, i. e., all nodes belonging to the same fiber has the same information input-tree.

    -----------------------------------------------------------------------------------------------------

    The code hasn't any command line arguments, requiring only two files: one single-line file cointaing a
    integer number representing the number of nodes in the network, and a 'edgelist.dat' file containing all
    the directed links between nodes (3 columns -> Pointing Node/ Pointed Node/ Type of regulation).

    The result is stored in two arrays: 'nodecolor' and 'edgecolor'. Each array store the color of each
    component (node or link) and, thus, gives all the information needed, together with the files, for the
    base graph construction.

    Author: Higor da S. Monteiro - Universidade Federal do Ceará
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include "graphdef.h" // Personal module containing helpful functions to handle with network data and other common operations.

#define SEQSIZE 200     // It's possible that the size of that string should be larger for larger networks.
int nfibers;

////////////////////////////////////////////////////////////////
struct strseq
{
    int node;
    char strseq[SEQSIZE];       
};
typedef struct strseq STRSEQ;

// Comparison function to be used in string sorting operation.
int cmp(const void *a, const void *b)
{
    STRSEQ *a1 = (STRSEQ *)a;
    STRSEQ *a2 = (STRSEQ *)b;

    return strcmp((*a1).strseq, (*a2).strseq);
}
/////////////////////////////////////////////////////////////////

//////////////////////// GRAPH OPERATIONS //////////////////////////
int COUNT_EXTREGULATORS(Graph* graph, int node, int* nodecolor)
{
    int l = 0;
    NodeAdj* Node = graph->array[node].head_in;
    while(Node)
    {
        if(nodecolor[Node->neighbor]==nodecolor[node]) l++;
        Node = Node->next_in;
    }
    return l;
}

int CHECKLINK(Graph* graph, int node1, int node2)
{
    NodeAdj* Node = graph->array[node1].head_out;
    while(Node)
    {
        if(Node->neighbor==node2) return 1;
        Node = Node->next_out;
    }
    return 0;
}
////////////////////////////////////////////////////////////////////

/*  After that all nodes and edges in the network have their colors set correctly, we upgrade 
    the values given in the table structure. */
void UPGRADETABLE(int M, int** table, int** edges, int* edgecolor, int* nodecolor)
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


/*  After the upgrade of the color table structure is done, we have to check if the minimal balance condition
    is satisfied. For this task, we loop for each color class, stacking all the nodes belonging to the 
    current color class. Then, we store and sort lexicographically the strings sequence to each node, in
    order to check if all the sequences are equal in the current color. If this is not true for at least
    one class, then the minimal balance condition is not satisfied.                                         */
int CHECKBALANCE(int ncolors, int** edges, int** table, int* nodecolor, int* edgecolor, int N, int M)
{
    int i,j;
    int nodetemp, nn, nsize;

    STRSEQ* seqlist;
    Stack* node_in_class = NULL;
    
    int newncolors = ncolors;
    for(i=0; i<ncolors; i++)    // For each color class C_i.
    {
        nn = 0;          // Number of nodes inside the color class 'i'.
        // Stores in a stack all the nodes belonging to C_i.
        for(j=0; j<N; j++) if(nodecolor[j]==i) { node_in_class = push(node_in_class, j); nn++; }

        nsize = nn;
        seqlist = (STRSEQ*)malloc(nn*sizeof(STRSEQ));
    
        // Get node at the top of the 'node_in_class' stack until stack is EMPTY.
        while(node_in_class)     
        {
            char strvar[SEQSIZE];
            strcpy(strvar, "");
            nodetemp = node_in_class->node_ID;
           
            //  For each color class C_j.
            for(j=0; j<ncolors; j++)
            {
                char temp[10];
                sprintf(temp, "%d", table[nodetemp][j]);
                strcat(strvar, temp);
            }
            seqlist[nsize-1].node = nodetemp;
            strcpy(seqlist[nsize-1].strseq, strvar);

            node_in_class = pop(node_in_class);
            nsize--;
        }

        // Sort lexicographically that list. The blocks that have the same sequence will be neighbors
        // so that we can find the number of partitions inside the class color C_i.

        qsort(seqlist, nn, sizeof(seqlist[0]), cmp);
        for(j=0; j<(nn-1); j++) if(strcmp(seqlist[j].strseq, seqlist[j+1].strseq)!=0) return 1;
            
        free(seqlist);
    }
    return 0;
}

/*  We define a string sequence for each node in a current color class. That sequence is all the
    number of inwards links of a node for each color placed side by side in a string. This way, when we
    sort lexicographically a list of the sequences for all nodes in class C_i, if there is different
    string sequences in that class, then equal sequences will be located side by side and it is possible
    to define new color classes inside C_i.                                                                 */
int CHECKSTATE(int ncolors, int** edges, int** table, int* nodecolor, int* edgecolor, int N, int M)
{
    int pluscolor;
    int newncolors;
    int i, j, nn, nodetemp, nsize;

    STRSEQ* seqlist;
    Stack* node_in_class = NULL;

    newncolors = ncolors;
    for(i=0; i<ncolors; i++)    // For each color class C_i.
    {
        nn = 0;          // Number of nodes inside the color class 'i'.
        pluscolor = 0;   // If the current color class is going to be partitioned, the 'pluscolor' stores the additional number of colors.
       
        // Stores in a stack all the nodes belonging to C_i.
        for(j=0; j<N; j++)   if(nodecolor[j]==i) { node_in_class = push(node_in_class, j); nn++; }

        seqlist = (STRSEQ*)malloc(nn*sizeof(STRSEQ));
    
        // Get node at the top of the 'node_in_class' stack until stack is EMPTY.
        nsize = nn;
        while(node_in_class)     
        {
            char strvar[SEQSIZE];
            strcpy(strvar, "");
            nodetemp = node_in_class->node_ID;
           
            //  For each color class C_j.
            for(j=0; j<ncolors; j++)
            {
                char temp[10];
                sprintf(temp, "%d", table[nodetemp][j]);
                strcat(strvar, temp);
            }
            seqlist[nsize-1].node = nodetemp;
            strcpy(seqlist[nsize-1].strseq, strvar);

            node_in_class = pop(node_in_class);
            nsize--;
        }

        // Now 'seqlist' is a list of structures cointaining the node in the color class and its
        // string sequence of the number of in-links for all colors.

        // Sort lexicographically that list. The blocks that have the same sequence will be neighbors
        // so that we can find the number of partitions inside the class color C_i.

        qsort(seqlist, nn, sizeof(seqlist[0]), cmp);

        for(j=0; j<(nn-1); j++)
        {
            if(strcmp(seqlist[j].strseq, seqlist[j+1].strseq)!=0) pluscolor++;
            
            nodecolor[seqlist[j+1].node] = (newncolors-1) + pluscolor;
        }
        free(seqlist);

        newncolors += pluscolor;
    }

    // Upgrades the color of each directed link.
    for(j=0; j<M; j++)  edgecolor[j] = nodecolor[edges[j][0]];

    return newncolors;
}

int** SETFIBERS(Graph* graph, int** edges, int* edgecolor, int* nodecolor, int** table, int ncolors, int N, int M)
{
    int i,j;
    int ncolortemp, minbalance;

    // Table for the initial state.
    table = arrint2d(N, ncolors);   
    UPGRADETABLE(M, table, edges, edgecolor, nodecolor);
    minbalance = CHECKBALANCE(ncolors, edges, table, nodecolor, edgecolor, N, M);

    while(minbalance==1)
    {
        ncolortemp = CHECKSTATE(ncolors, edges, table, nodecolor, edgecolor, N, M);

        free2d(table, N, ncolors);
        ncolors = ncolortemp;
        table = arrint2d(N, ncolors);
        UPGRADETABLE(M, table, edges, edgecolor, nodecolor);

        minbalance = CHECKBALANCE(ncolors, edges, table, nodecolor, edgecolor, N, M);
    }

    printf("%d\n", ncolors);
    for(i=0; i<N; i++) printf("%d %d \n", i, nodecolor[i]);
    printf("\n");
    nfibers = ncolors;
}

int main() 
{ 
    int N;                                                  // Number of nodes in the network.
    
    char netsize[100] = "../Data/Ecoli/ngenes.dat";         // File containing (one line) the number of nodes in the network.
    char net_edges[100] = "../Data/Ecoli/edgelist.dat";   // File containing all the directed links in the network.
    
    // Defines the size of the network.
    FILE* UTIL = fopen(netsize, "r");
    fscanf(UTIL, "%d\n", &N);
    fclose(UTIL);
    //////////////////////////////////

    // Creates the graph structure for N nodes (Not used in the balanced coloring algorithm)
    Graph* graph = createGraph(N);

    // Properly defines the network structure with the given 'edgelist.dat' file.
    int** edges;
    int* regulator;
    edges = defineNetwork(edges, regulator, graph, net_edges);
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////// Minimal balanced coloring algorithm ////////////////////////
	int M = nlines_file(net_edges, 3);              // Number of edges.  		
    
    int i, j, k;

    // INITIAL STATE: ALL NODES AND LINKS HAVE THE SAME COLOR '0' //
    int ncolors = 1;
    int* nodecolor = (int*)malloc(N*sizeof(int));
    int* edgecolor = (int*)malloc(M*sizeof(int));
    for(i=0; i<N; i++) nodecolor[i] = 0;        
    for(j=0; j<M; j++) edgecolor[j] = 0;

    int** table = SETFIBERS(graph, edges, edgecolor, nodecolor, table, ncolors, N, M);
    ////////////////////////////////////////////////////////////////////////////////////       

    //////////////////////////////// FIBER STATISTICS //////////////////////////////////
    int* nloop = (int*)malloc(nfibers*sizeof(int));
    int* lextreg = (int*)malloc(nfibers*sizeof(int));

    // For each fiber C_i, we define its adjacency matrix in order to calculate the vector <n,l|.
    int* temp_index;
    double* temp_adjmatrix;
    int nexternal, nn, neigh;
    Stack* node_in_class;
    for(i=0; i<nfibers; i++)
    {
        nn = 0;
        nexternal = 0;
        node_in_class = NULL;
        for(j=0; j<N; j++) if(nodecolor[j]==i) { node_in_class = push(node_in_class, j); nn++; }

        temp_index = (int*)malloc(nn*sizeof(int));
        temp_adjmatrix = (double*)malloc(nn*nn*sizeof(double));

        for(j=0; j<nn; j++) { temp_index[j] = node_in_class->node_ID; node_in_class = pop(node_in_class); }

        int index = 0;
        for(j=0; j<nn; j++)
        {
            int aux = COUNT_EXTREGULATORS(graph, temp_index[j], nodecolor);
            nexternal += aux;
            for(k=0; k<nn; k++)
            {
                neigh = CHECKLINK(graph, temp_index[j], temp_index[k]);
                if(neigh==1) temp_adjmatrix[j*nn + k] = 1.0;
                else temp_adjmatrix[j*nn + k] = 0.0;
            }
        }

        gsl_matrix_view m = gsl_matrix_view_array(temp_adjmatrix, nn, nn);
        gsl_vector_complex *eval = gsl_vector_complex_alloc (nn);

        gsl_eigen_nonsymm_workspace* w = gsl_eigen_nonsymm_alloc(nn);
        gsl_eigen_nonsymm(&m.matrix, eval, w);

        gsl_eigen_nonsymm_free(w);
        printf("eigenvalues of fiber %d:", i);
        for(j=0; j<nn; j++) { double eval_j = gsl_vector_get(eval, j); printf("%lf ", eval_j); }
        printf("\n");
        printf("%d\n", nexternal);

        lextreg[i] = nexternal;
    } 

}