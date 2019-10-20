#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
};
typedef struct NodeAdj NodeAdj;

struct neigh_vector
{
	int neigh;
	int color;
};
typedef struct neigh_vector neigh_vector;

struct stack
{
	int extn;
	int node_ID;
	struct stack *next;
};
typedef struct stack Stack;

struct DoublyLinkNode
{
    int data;
    struct DoublyLinkNode* prev;
    struct DoublyLinkNode* next;
};
typedef struct DoublyLinkNode DoublyLinkNode;

struct Block
{
    int size;
    DoublyLinkNode* head;   // The doubly linked list of the block containing the elements on it.
};
typedef struct Block BLOCK;

struct QueueOfBlocks
{
    BLOCK* block;
    struct QueueOfBlocks* next;
};
typedef struct QueueOfBlocks QBLOCK;
