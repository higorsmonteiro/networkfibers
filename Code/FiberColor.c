/*  The code given here has the aim to determine the minimal balanced color classes to identify
    the fibration building blocks of a directed network. The theory concerning the graph fibration
    morphism and its application on biological networks are detailed, respectively, mainly in the
    paper of Boldi and Vigna (2001) and in the paper of Morone et. al. (2019, to be published). Concerning
    the algorithm reproduced in this code, one should refer to the work of Cardon and Crochemore (1982).

    The resulted fibers represent, for each fiber, the set of nodes that are isomorphic under a graph
    fibration morphism, i. e., all nodes belonging to the same fiber has the same information input-tree.

    -----------------------------------------------------------------------------------------------------

    The code hasn't any command line arguments, requiring only two files: one single-line file containing a
    integer number representing the number of nodes in the network, and a 'edgelist.dat' file containing all
    the directed links between nodes (3 columns -> Pointing Node/ Pointed Node/ Type of regulation).

    The result is stored in two arrays: 'nodecolor' and 'edgecolor'. Each array store the color of each
    component (node or link) and, thus, gives all the information needed, together with the files, for the
    base graph construction.

    Author: Higor da S. Monteiro - Departament of Physics - Universidade Federal do Ceará (UFC)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_eigen.h>
#include "graphdef.h"   // Personal module containing helpful functions to handle with network data and other common operations.
#define SEQSIZE 2000     // It's possible that the size of the string sequences should be larger for larger networks.

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
int INSERT_EXTREGULATORS(Stack* extregul, Stack* block, Graph* graph, int nn, int* nodecolor, int curcolor)
{
    int i;
    int n = nn;
    Stack* Temp = block;
    int* temp = (int*)malloc(n*sizeof(int));
    // 'temp' stores all the nodes of a color class.
    for(i=0; i<n; i++) { temp[i] = Temp->node_ID; Temp = Temp->next; }

    
    // We store in 'Temp' all the external nodes that point to any node in the current color class.
    int extn = 0;
    extregul = NULL;
    for(i=0; i<n; i++)
    {
        NodeAdj* Node = graph->array[temp[i]].head_in;
        while(Node)
        {
            if(nodecolor[Node->neighbor]!=curcolor) { extregul = push(extregul, Node->neighbor); extn++; }
            Node = Node->next_in;
        }
    }

    // 'extregul' may have repeated elements.
    int* EXT = (int*)malloc(extn*sizeof(int));
    int* aux = (int*)malloc(extn*sizeof(int));
    for(i=0; i<extn; i++) aux[i] = -1;
    for(i=0; i<extn; i++) { EXT[i] = extregul->node_ID; extregul = pop(extregul); }

    int nn_add = 0;
    for(i=0; i<extn; i++)
    {
        int check = INCHECK(aux, extn, EXT[i]);
        if(check==0)
        {
            aux[nn_add++] = EXT[i];
            block = push(block, EXT[i]);
        }
    }
    return nn_add;
}

int COUNT_EXTREGULATORS(Graph* graph, int node, int* nodecolor)
{
    int l = 0;
    NodeAdj* Node = graph->array[node].head_in;
    while(Node)
    {
        if(nodecolor[Node->neighbor]!=nodecolor[node]) l++;
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

/*  After all nodes and edges in the network have their colors set correctly, we upgrade 
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
    int added_color, newncolor;
    int i, j, nn, nodetemp, nsize;

    STRSEQ* seqlist;
    Stack* node_in_class = NULL;    // Stack data structure to store the nodes of a specific color class.  

    newncolor = ncolors;
    for(i=0; i<ncolors; i++)    // For each color class C_i.
    {
        nn = 0;             // Number of nodes inside the color class 'i'.
        added_color = 0;   // If the current color class is going to be partitioned, it stores the additional number of colors.
       
        // Stores in a stack all the nodes belonging to C_i.
        for(j=0; j<N; j++)   if(nodecolor[j]==i) { node_in_class = push(node_in_class, j); nn++; }

        // Get node at the top of the 'node_in_class' stack until stack is EMPTY.
        nsize = nn;
        seqlist = (STRSEQ*)malloc(nn*sizeof(STRSEQ));
        while(node_in_class)     
        {
            char strvar[SEQSIZE];
            strcpy(strvar, "");                 // Initially an empty sequence.
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
            if(strcmp(seqlist[j].strseq, seqlist[j+1].strseq)!=0) added_color++;
            if(added_color>0) nodecolor[seqlist[j+1].node] = (newncolor-1) + added_color;
        }
        newncolor += added_color;
        free(seqlist);
    }

    // Upgrades the color of each directed link.
    for(j=0; j<M; j++)  edgecolor[j] = nodecolor[edges[j][0]];
    return newncolor;
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
        printf("Fibers: %d\n", ncolors);
        //for(j=0; j<ncolors; j++)
        //{
        //    printf("Fiber number %d:", j);
        //    for(i=0; i<N; i++)
        //        if(nodecolor[i]==j) printf("%d ", i);
        //    printf("\n");
        //}
        //printf("\n\n");
    }
    nfibers = ncolors;
}

int main() 
{ 
    int N;                                                  // Number of nodes in the network.
    char netsize[100] = "../Data/ECOLINgenes.dat";          // File containing (one line) the number of nodes in the network.
    char net_edges[100] = "../Data/ECOLIedgelist.dat";      // File containing all the directed links in the network.
    
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
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////// MINIMAL BALANCED COLORING ALGORITHM ////////////////////////
    int i, j, k;
    int M = nlines_file(net_edges, 3);              // Number of edges.  		
    

    // INITIAL STATE: All nodes and links have the same color '0' //
    int ncolors = 1;
    int* nodecolor = (int*)malloc(N*sizeof(int));
    int* edgecolor = (int*)malloc(M*sizeof(int));
    for(i=0; i<N; i++) nodecolor[i] = 0;        
    for(j=0; j<M; j++) edgecolor[j] = 0;

    int** table = SETFIBERS(graph, edges, edgecolor, nodecolor, table, ncolors, N, M);
    ////////////////////////////////////////////////////////////////////////////////////       

    //////////////////////////////// FIBER STATISTICS //////////////////////////////////
    int* lextreg = (int*)malloc(nfibers*sizeof(int));
    double* nloop = (double*)malloc(nfibers*sizeof(double));

    // For each fiber C_i, we define its adjacency matrix to calculate the fiber vector <n,l|.
    int* temp_index;
    double* temp_adjmatrix;

    Stack* blocknodes;
    Stack* repeated_extregul;
    int nn_add, nexternal, nn, neigh;
    for(i=0; i<nfibers; i++)
    {
        nn = 0;
        nexternal = 0;

        // Stack all the nodes belonging to the current color class and its external regulators. 
        printf("%d\n", nexternal);
        blocknodes = NULL;
        repeated_extregul = NULL;
        for(j=0; j<N; j++)
        {
            if(nodecolor[j]==i)
            {
                nn++;
                blocknodes = push(blocknodes, j);
            }
        }
        nexternal = INSERT_EXTREGULATORS(repeated_extregul, blocknodes, graph, nn, nodecolor, i);
        nn += nexternal;

        temp_index = (int*)malloc(nn*sizeof(int));
        temp_adjmatrix = (double*)malloc(nn*nn*sizeof(double));
        for(j=0; j<nn; j++) { temp_index[j] = blocknodes->node_ID; blocknodes = pop(blocknodes); }

        // Now we construct the current fiber adjacency matrix to calculate its eigenvalues.
        for(j=0; j<nn; j++)
        {
            for(k=0; k<nn; k++)
            {
                neigh = CHECKLINK(graph, temp_index[j], temp_index[k]); // A_{jk} = 1 if j -> k
                if(neigh==1) temp_adjmatrix[j*nn + k] = 1.0;
                else temp_adjmatrix[j*nn + k] = 0.0;
            }
        }

        /////////// GSL PACKAGE ROUTINES TO EIGENSYSTEMS PROBLEMS ////////////
        gsl_matrix_view m = gsl_matrix_view_array(temp_adjmatrix, nn, nn);
        // 'eval' will gonna stores all the 'nn' eigenvalues of the fiber adjacency matrix.
        gsl_vector_complex *eval = gsl_vector_complex_alloc (nn);

        gsl_eigen_nonsymm_workspace* w = gsl_eigen_nonsymm_alloc(nn);
        gsl_eigen_nonsymm(&m.matrix, eval, w);
        gsl_eigen_nonsymm_free(w);
        //////////////////////////////////////////////////////////////////////

        double temp;
        double eigmax = gsl_vector_get(eval, 0);
        for(j=1; j<nn; j++) { temp = gsl_vector_get(eval, j); if(temp>eigmax) eigmax = temp; }

        nloop[i] = eigmax;
        lextreg[i] = nexternal;

        free(temp_index);
        free(temp_adjmatrix);
    }

    for(i=0; i<nfibers; i++) printf("Fiber %d: n = %lf, l = %d \n", i, nloop[i], lextreg[i]); 

}