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
};
typedef struct aux AUX;

struct strseq
{
    int node;
    Node* strseq;       
};
typedef struct strseq STRSEQ;

STRSEQ* INITSEQ(STRSEQ* seqlist, int N)
{
	int i;
	for(i=0; i<N; i++) seqlist[i].strseq = NULL;
	return seqlist;
}

int main(int argv, char** argc) 
{ 
    int N, i;                              // Number of nodes in the network.
                                                          
	Node* List = NULL;
	STRSEQ* seqlist = (STRSEQ*)malloc(10*sizeof(STRSEQ));
	seqlist = INITSEQ(seqlist, 10);

	int arr[] = {5, 3, 2, 9, 4, 1, 8, 6, 7, 0};

	for(i=0; i<10; i++)
	{
		List = push_data(List, arr[i]);
		seqlist[i].node = i;
	}
	
	Node* temp = List;
	while(temp)
	{
		printf("%d ", temp->data);
		temp = temp->next;
	}

	SelectionSort(List);
	temp = List;	
	while(temp)
	{
		printf("%d ", temp->data);
		temp = temp->next;
	}

}